library(ggplot2)
library(ggrepel)
library(dplyr)

setwd("rna-seq_workflow/scripts/")


# Import Files:
# degs = differentially expressed genes (results from a software like limma-voom)
degs = read.csv("../data_files/limma_voom/POLvsSAL.csv")

# N.B. important to distinguish 3 items used in code here,
# (these will be named differently according to which pipeline data is from)
## the adjusted p-value.
## the logFC. (this can be variants such as log2FC, , etc.)
## the non-adjusted (normal) p-values.


# Prepare data:
degs$updown <- ifelse(degs$adj.P.Val > 0.05, 'Not Sig',
                      ifelse(degs$logFC < 0, 'down', 'up'))


# Plot Volcano plot

## No labels:
p <- ggplot(degs, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = updown)) + 
  scale_color_manual(values = c("blue", "gray", "red")) +
  theme_bw() + 
  theme(
    # legend.position = "bottom",
    legend.position = c(0.1, 0.825),
    legend.background = element_rect(linetype="solid", colour ="black", size=0.5),
  )

p

# Save to both PNG and SVG formats:
ggsave("../sample_figures/volcano_plot_degs.png", p, width = 6, height = 4)
ggsave("../sample_figures/volcano_plot_degs.svg", p, width = 11, height = 8.5)


## With genes labeled:
pl <- ggplot(degs, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = updown)) + 
  scale_color_manual(values = c("blue", "gray", "red")) +
  theme_bw() + 
  theme(
    # legend.position = "bottom",
    legend.position = c(0.1, 0.825),
    legend.background = element_rect(linetype="solid", colour ="black", size=0.5),
  ) +
  geom_text_repel(data = subset(degs, adj.P.Val < 0.005 & (logFC >= 1 | logFC <= 1)),
                  aes(label = external_gene_name),
                  size = 2.5,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")
  )

pl

# Save to both PNG and SVG formats:
ggsave("../sample_figures/volcano_plot_degs_labeled.png", pl, width = 6, height = 4)
ggsave("../sample_figures/volcano_plot_degs_labeled.svg", pl, width = 11, height = 8.5)
