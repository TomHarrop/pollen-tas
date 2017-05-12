#!/usr/bin/env Rscript

library(data.table)

# parse taxonomy files 
GenerateQiimeTaxonomyLine <- function(
    x,
    allowed_taxa = allowed_taxa,
    ranks = ranks) {
    # remove blank entries
    all_names <- as.character(x[x != ''])
    
    # subset allowed_taxa by taxon names in x
    matched_taxa <- allowed_taxa[x]
    
    # make sure every rank in ranks is represented
    all_ranks <- merge(data.table(rank = ranks),
                       matched_taxa,
                       by = "rank",
                       all.x = TRUE)
    all_ranks[, rank := factor(rank, levels = ranks)]
    setorder(all_ranks, rank)
    
    # replace NA with blanks
    all_ranks[is.na(name_txt), name_txt := ""]
    
    # use qiime format for taxa names (prefix with node class letter)
    all_ranks[, name_qiime := paste(substr(rank, 1, 1), name_txt, sep = "_")]
    return(all_ranks[, paste0(name_qiime, collapse = "; ")])
}

# load data
rutils::GenerateMessage("Loading tax_names")
tax_names <- readRDS("data/ncbi/names.dmp.Rds")
rutils::GenerateMessage("Loading tax_nodes")
tax_nodes <- readRDS("data/ncbi/nodes.dmp.Rds")

# generate list of acceptable names
rutils::GenerateMessage("Generating allowed_taxa")
ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")
tax_ids_to_keep <- tax_nodes[rank %in% ranks, unique(tax_id)]
filtered_tax_names <- tax_names[tax_id %in% tax_ids_to_keep &
                                    name_class == "scientific name"]
allowed_taxa <- merge(filtered_tax_names[, .(tax_id, name_txt)],
                      tax_nodes[, .(tax_id, rank)],
                      all.x = TRUE,
                      all.y = FALSE)
setkey(allowed_taxa, name_txt)

# taxonomy from genbank rbcL search results
rutils::GenerateMessage("Loading GenBank results")
acc_tax <- fread("output/taxonomy/genbank_results.txt")

# generate a taxonomy line for each accession
rutils::GenerateMessage("Searching for taxa in allowed_taxa")
taxo <- acc_tax[, .(
    taxonomy = GenerateQiimeTaxonomyLine(
        taxon, allowed_taxa, ranks)),
    by = accession]

# add species from acc_tax
acc_spec <- unique(acc_tax[, .(accession, species)])
taxo_spec <- merge(taxo,
                   acc_spec,
                   all.x = TRUE,
                   all.y = FALSE,
                   by = "accession")
taxo_spec[, species := gsub("[^[:alnum:]]+", ".", species)]
taxo_spec[, taxonomy := paste(taxonomy, species, sep = "; s_")]

# write output
rutils::GenerateMessage("Writing output")
fwrite(taxo_spec,
       file = "output/taxonomy/taxonomy.txt",
       quote = FALSE,
       sep = "\t",
       col.names = FALSE)

