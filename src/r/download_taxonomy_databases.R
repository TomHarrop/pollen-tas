#!/usr/bin/env Rscript

library(data.table)

# download taxdump (complete NCBI taxonomy)
rutils::GenerateMessage("Downloading taxdump")
taxdmp.url <- "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"
temp2 <- tempfile(fileext = ".zip")
temp3 <- tempdir()
download.file(taxdmp.url, temp2)

# read nodes.dmp from inside taxdmp
rutils::GenerateMessage("fread-ing nodes.dmp")
nodes.dmp.file <- unzip(temp2, files = "nodes.dmp", exdir = temp3)
nodes.dmp.raw <- fread(nodes.dmp.file, sep = "|", quote = "\t", header = FALSE)

# tidy up nodes.dmp
rutils::GenerateMessage("Munging nodes.dmp")
nodes.dmp.raw[, V14 := NULL]
nodes.names <- c("tax_id", "parent_tax_id", "rank", "embl_code", "division_id",
                 "inherited_div_flag", "genetic_code_id", "inherited_GC",
                 "mitochondrial_genetic_code_id", "inherited_MGC_flag",
                 "GenBank_hidden_flag", "hidden_subtree_root_flag", "comments")
setnames(nodes.dmp.raw, names(nodes.dmp.raw), nodes.names)
nodes.dmp.raw[, tax_id := as.numeric(gsub("\t", "", tax_id, fixed = TRUE))]
setkey(nodes.dmp.raw, tax_id)

# read names.dmp from inside taxdmp
rutils::GenerateMessage("fread-ing names.dmp")
names.dmp.file <- unzip(temp2, files = "names.dmp", exdir = temp3)
names.dmp.raw <- fread(names.dmp.file, sep = "|", quote = "\t", header = FALSE)

# tidy up names.dmp
rutils::GenerateMessage("Munging names.dmp")
names.dmp.raw[, V5 := NULL]
names.names <- c("tax_id", "name_txt", "unique_name", "name_class")
setnames(names.dmp.raw, names(names.dmp.raw), names.names)
names.dmp.raw[, tax_id := as.numeric(gsub("\t", "", tax_id, fixed = TRUE))]
setkey(names.dmp.raw, tax_id)

# read merged.dmp from taxdmp
rutils::GenerateMessage("fread-ing merged.dmp")
merged.dmp.file <- unzip(temp2, files = "merged.dmp", exdir = temp3)
merged.dmp.raw <- fread(merged.dmp.file, sep = "|", quote = "\t", header = FALSE)

# tidy up merged.dmp
rutils::GenerateMessage("Munging merged.dmp")
merged.dmp.raw[, V3 := NULL]
merged.names <- c("old_tax_id", "new_tax_id")
setnames(merged.dmp.raw, names(merged.dmp.raw), merged.names)
merged.dmp.raw[, old_tax_id :=
                   as.numeric(gsub("\t", "", old_tax_id, fixed = TRUE))]
merged.dmp.raw[, new_tax_id := as.numeric(new_tax_id)]

# save output
outdir <- "data/ncbi"
rutils::GenerateMessage("Writing output")
rutils::PrintF("outdir: %s\n", outdir)
if (!dir.exists(outdir)) {
    dir.create(outdir)
}

saveRDS(nodes.dmp.raw, paste0(outdir, "/nodes.dmp.Rds"))
saveRDS(names.dmp.raw, paste0(outdir, "/names.dmp.Rds"))

# log metadata
rutils::GenerateMessage("Logging metadata")
log.file <- paste0(outdir, "/SessionInfo.txt")
s.inf <- rutils::GitSessionInfo()
writeLines(s.inf, log.file)

quit("no", 0)
