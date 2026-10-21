#!/usr/bin/env Rscript
# Dotplot of one chromosome against its counterpart in another assembly, with
# axes in Mb. Small isolated alignments can be dropped with min_len to remove
# scattered noise.
#
# Usage:
#   Rscript R/dotplot_chromosome.R <coords> <ref_chrom> <qry_chrom> <ref_label> <qry_label> \
#           <outprefix> [min_len=10000] [tick_mb=5]

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("usage: dotplot_chromosome.R <coords> <ref_chrom> <qry_chrom> <ref_label> <qry_label> <outprefix> [min_len] [tick_mb]")
}
coords_file <- args[1]
ref_chrom   <- args[2]
qry_chrom   <- args[3]
ref_label   <- args[4]
qry_label   <- args[5]
outprefix   <- args[6]
min_len     <- if (length(args) >= 7) as.numeric(args[7]) else 10000
tick_mb     <- if (length(args) >= 8) as.numeric(args[8]) else 5

lines <- readLines(coords_file)
first <- which(grepl("^[0-9]", lines))[1]
coords <- read.table(coords_file, skip = first - 1, sep = "\t", header = FALSE,
                     fill = TRUE, strip.white = TRUE)
colnames(coords)[1:13] <- c("ref_start", "ref_end", "qry_start", "qry_end",
                            "len_ref", "len_qry", "pct_id", "chr_len_ref", "chr_len_qry",
                            "cov_ref", "cov_qry", "chr_ref", "chr_qry")

d <- coords %>%
  filter(chr_ref == ref_chrom, chr_qry == qry_chrom, len_ref >= min_len) %>%
  mutate(orientation = ifelse((ref_end - ref_start) * (qry_end - qry_start) > 0, "Forward", "Reverse"),
         x = (ref_start + ref_end) / 2 / 1e6,
         y = (qry_start + qry_end) / 2 / 1e6)
if (nrow(d) == 0) stop("No alignments between ", ref_chrom, " and ", qry_chrom)

xmax <- ceiling(max(d$chr_len_ref) / 1e6)
ymax <- ceiling(max(d$chr_len_qry) / 1e6)

p <- ggplot(d, aes(x = x, y = y, colour = orientation)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_colour_manual(values = c(Forward = "red", Reverse = "blue"), name = "Orientation") +
  scale_x_continuous(breaks = seq(0, xmax, tick_mb), limits = c(0, xmax), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, ymax, tick_mb), limits = c(0, ymax), expand = c(0, 0)) +
  labs(title = paste(ref_label, "vs", qry_label),
       subtitle = paste("Alignments:", nrow(d), "| mean identity:", round(mean(d$pct_id), 1), "%"),
       x = paste(ref_label, "position (Mb)"),
       y = paste(qry_label, "position (Mb)")) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(paste0(outprefix, ".png"), p, width = 8, height = 8, dpi = 300)
ggsave(paste0(outprefix, ".pdf"), p, width = 8, height = 8)

summary <- data.frame(pair = paste(ref_label, "vs", qry_label),
                      alignments = nrow(d),
                      forward = sum(d$orientation == "Forward"),
                      reverse = sum(d$orientation == "Reverse"),
                      mean_identity = round(mean(d$pct_id), 2),
                      ref_covered_bp = sum(abs(d$ref_end - d$ref_start)))
write.table(summary, paste0(outprefix, ".summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary)
