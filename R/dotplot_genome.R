#!/usr/bin/env Rscript
# Whole-genome dotplot from a show-coords -T -r -c -l table, with chromosomes
# laid end to end and labelled Chr1..ChrN, plus a chromosome correspondence
# table (which query chromosome gets most aligned sequence for each
# reference chromosome) to check that chromosome numbering agrees between
# assemblies.
#
# Usage:
#   Rscript R/dotplot_genome.R <coords> <ref_label> <qry_label> <outprefix> [chrom_regex]

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: dotplot_genome.R <coords> <ref_label> <qry_label> <outprefix> [chrom_regex]")
}
coords_file <- args[1]
ref_label   <- args[2]
qry_label   <- args[3]
outprefix   <- args[4]
chrom_regex <- if (length(args) >= 5) args[5] else NA

read_coords <- function(path) {
  lines <- readLines(path)
  first <- which(grepl("^[0-9]", lines))[1]
  d <- read.table(path, skip = first - 1, sep = "\t", header = FALSE,
                  fill = TRUE, strip.white = TRUE)
  colnames(d)[1:13] <- c("ref_start", "ref_end", "qry_start", "qry_end",
                         "len_ref", "len_qry", "pct_id", "chr_len_ref", "chr_len_qry",
                         "cov_ref", "cov_qry", "chr_ref", "chr_qry")
  d %>% filter(!is.na(ref_start), !is.na(qry_start))
}

# Order chromosomes by the trailing number in their name
chrom_number <- function(x) as.numeric(sub(".*?(\\d+)$", "\\1", x))

coords <- read_coords(coords_file)
if (!is.na(chrom_regex)) {
  coords <- coords %>% filter(grepl(chrom_regex, chr_ref), grepl(chrom_regex, chr_qry))
}

axis_table <- function(chr, len) {
  data.frame(chr = chr, len = len) %>%
    distinct() %>%
    arrange(chrom_number(chr)) %>%
    mutate(offset = lag(cumsum(as.numeric(len)), default = 0),
           mid = offset + len / 2,
           label = paste0("Chr", chrom_number(chr)))
}
ref_axis <- axis_table(coords$chr_ref, coords$chr_len_ref)
qry_axis <- axis_table(coords$chr_qry, coords$chr_len_qry)

coords <- coords %>%
  left_join(ref_axis %>% select(chr_ref = chr, ref_off = offset), by = "chr_ref") %>%
  left_join(qry_axis %>% select(chr_qry = chr, qry_off = offset), by = "chr_qry") %>%
  mutate(orientation = ifelse((ref_end - ref_start) * (qry_end - qry_start) > 0, "Forward", "Reverse"),
         x = (ref_start + ref_end) / 2 + ref_off,
         y = (qry_start + qry_end) / 2 + qry_off)

p <- ggplot(coords, aes(x = x, y = y, colour = orientation)) +
  geom_point(size = 0.9, alpha = 0.6) +
  scale_colour_manual(values = c(Forward = "red", Reverse = "blue"), name = "Match orientation") +
  geom_vline(xintercept = ref_axis$offset + ref_axis$len, colour = "grey70", linewidth = 0.4) +
  geom_hline(yintercept = qry_axis$offset + qry_axis$len, colour = "grey70", linewidth = 0.4) +
  scale_x_continuous(breaks = ref_axis$mid, labels = ref_axis$label, name = ref_label, expand = c(0, 0)) +
  scale_y_continuous(breaks = qry_axis$mid, labels = qry_axis$label, name = qry_label, expand = c(0, 0)) +
  labs(title = paste(ref_label, "vs", qry_label, "whole-genome alignment"),
       subtitle = paste("Forward:", sum(coords$orientation == "Forward"),
                        "| Reverse:", sum(coords$orientation == "Reverse"))) +
  coord_fixed() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = c(0.98, 0.02),
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(colour = "black"))

ggsave(paste0(outprefix, ".png"), p, width = 14, height = 12, dpi = 300)
ggsave(paste0(outprefix, ".pdf"), p, width = 14, height = 12)

# Chromosome correspondence: aligned bp per reference/query chromosome pair
corr <- coords %>%
  group_by(chr_ref, chr_qry) %>%
  summarise(aligned_bp = sum(abs(ref_end - ref_start)), alignments = n(),
            mean_identity = round(mean(pct_id), 2), .groups = "drop") %>%
  group_by(chr_ref) %>%
  mutate(share_of_ref = round(aligned_bp / sum(aligned_bp), 3)) %>%
  arrange(chrom_number(chr_ref), desc(aligned_bp)) %>%
  ungroup()
best <- corr %>% group_by(chr_ref) %>% slice_max(aligned_bp, n = 1) %>% ungroup()

write.table(corr, paste0(outprefix, ".correspondence_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(best, paste0(outprefix, ".correspondence_best.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote", paste0(outprefix, c(".png", ".pdf", ".correspondence_all.tsv", ".correspondence_best.tsv")), sep = "\n")
print(as.data.frame(best))
