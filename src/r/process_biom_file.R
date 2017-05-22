#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

# read data
raw_otu_data <- fread("output/qiime/open_references/otu_table_mc2_w_tax.txt")

# generate per-sample amplicons per taxa
otu_data_long <- melt(raw_otu_data,
                      id.vars = c("#OTU ID", "taxonomy"),
                      variable.name = "sample",
                      value.name = "amplicons")

amplicons_per_taxa <- otu_data_long[, .(
    amplicons = sum(amplicons)
), by = .(sample, taxonomy)]
amplicons_per_taxa[, percent_sample_total := 100 * amplicons / sum(amplicons),
                   by = sample]
amplicons_per_taxa[, total_reads_per_sample := sum(amplicons), by = sample]

# try to calculate psedo-size factors: median of count / average count by
# sample.
amp_wide <- dcast(amplicons_per_taxa,
                  taxonomy ~ sample,
                  value.var = "amplicons")
count.matrix <- as.matrix(data.frame(amp_wide, row.names = "taxonomy"))
cm <- count.matrix[rowSums(count.matrix) > 10, ]
sjmat <- cm / rowMeans(cm)
sf <- apply(sjmat, 2, function(x)
    median(x[is.finite(x) & x != 0], na.rm = TRUE))
# ggplot2::ggplot(data.frame(sample = names(sf), sf = sf),
#                 aes(x = sample, y = sf)) +
#     ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
#     ggplot2::geom_col()

# write output
fwrite(amplicons_per_taxa, "output/otu_tables/amplicons_per_taxon.csv")

# top 10 overall taxa
ranked_taxa <- amplicons_per_taxa[, .(
    total_amplicons = sum(amplicons)),
    by = taxonomy]
setorder(ranked_taxa, -total_amplicons)
plot_taxa <- ranked_taxa[1:10, taxonomy]

# quick plot
pd <- amplicons_per_taxa[taxonomy %in% plot_taxa]
pd[, taxonomy := factor(taxonomy, levels = rev(plot_taxa))]
so_tab <- pd[taxonomy == plot_taxa[1]]
setorder(so_tab, -percent_sample_total)
so <- so_tab[, unique(as.character(sample))]
pd[, sample := factor(sample, levels = so)]

g <- ggplot(pd, aes(x = sample, y = percent_sample_total, fill = taxonomy)) +
    theme_minimal() + xlab(NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    geom_col()

ggsave(filename = "output/otu_tables/quick_plot_top_10_taxa.pdf",
       plot = g,
       width = 10,
       height = 7.5)

