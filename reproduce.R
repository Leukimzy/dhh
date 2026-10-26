# 从仓库根目录运行：Rscript reproduce.R
# 三图所用 CSV 均来自论文补充表或公开原始数据，不从原图估读。

dir.create("results", showWarnings = FALSE)

# Fig.3h：小变异的频率组成 -----------------------------------------------
h <- read.csv("data/fig3h_five_classes.csv")
methods <- c("MC", "PGGB", "Dipcall", "HiFi-DeepVariant")
categories <- c("Common", "Low frequency", "Rare", "Singleton", "Missing")
bar_colors <- list(
  MC = c("#DF3475", "#EC719F", "#F4ADC6", "#FBE2EC", "white"),
  PGGB = c("#37B4CD", "#82CDDA", "#BCE4EA", "#E8F4F6", "white"),
  Dipcall = c("#FFB817", "#FFD17C", "#FFE3AC", "#FFF3D8", "white"),
  `HiFi-DeepVariant` = c("#6268A6", "#969AC8", "#C3C6DF", "#E9EAF4", "white")
)
png("results/fig3h.png", width = 1500, height = 650, res = 200)
par(mar = c(4.5, 8.5, 2.4, 1), bg = "white", xaxs = "i", yaxs = "i")
plot(0, 0, type = "n", xlim = c(0, 100), ylim = c(.5, 5.2),
     axes = FALSE, xlab = "", ylab = "")
for (i in 1:4) {
  one <- h[h$method == methods[i], ]
  one <- one[match(categories, one$category), ]
  stopifnot(abs(sum(one$percent) - 100) < 1e-6)
  left <- 0
  y <- 5 - i
  for (j in 1:5) {
    right <- left + one$percent[j]
    rect(left, y - .28, right, y + .28,
         col = bar_colors[[methods[i]]][j], border = "grey20")
    left <- right
  }
}
axis(1, at = c(0, 50, 100), las = 1)
axis(2, at = 4:1, labels = methods, las = 1, tick = FALSE)
title(xlab = "Fraction of small variants (%)")
text(c(7, 29, 53, 77, 93), 5.05, categories, cex = .75)
rect(c(4, 26, 50, 74, 90), 4.69, c(10, 32, 56, 80, 96), 4.84,
     col = c("grey65", "grey78", "grey85", "grey94", "white"), border = "grey20")
dev.off()

# Fig.3b：MC、PGGB 的小变异位点数量 --------------------------------------
d <- read.csv("data/fig3b_counts_wide.csv")
stopifnot(nrow(d) == 6)
stopifnot(all(d$type == rep(c("SNP", "MNP", "INDEL"), each = 2)))
stopifnot(all(d$method == rep(c("MC", "PGGB"), 3)))
two <- d$two / 1e6
three <- d$three / 1e6
more <- d$more / 1e6
method_col <- rep(c("#E33772", "#43B8CE"), 3)
png("results/fig3b.png", width = 1800, height = 900, res = 200)
par(mar = c(5, 5, 2, 1), bg = "white", las = 1, bty = "l")
mid <- barplot(two, col = method_col, border = "grey35", ylim = c(0, 25),
               space = c(.12, .07, .85, .07, .85, .07), axes = FALSE,
               ylab = "Number of small variant sites (x10^6)")
rect(mid - .5, two, mid + .5, two + three, col = "grey82", border = "grey35")
rect(mid - .5, two + three, mid + .5, two + three + more,
     col = "white", border = "grey35")
axis(2, at = c(0, 5, 10, 15, 20, 25), las = 1)
axis(1, at = c(mean(mid[1:2]), mean(mid[3:4]), mean(mid[5:6])),
     labels = c("SNP", "MNP", "INDEL"), tick = FALSE)
legend("topright", legend = c("MC", "PGGB", "3 alleles", ">3 alleles"),
       fill = c("#E33772", "#43B8CE", "grey82", "white"),
       border = "grey35", cex = .82, bty = "n", ncol = 2)
dev.off()

# Fig.6c：逐位点频率与六角形密度 -----------------------------------------
if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("hexbin", quietly = TRUE)) {
  stop("Fig.6c needs ggplot2 and hexbin: install.packages(c('ggplot2','hexbin'))")
}
freq <- read.csv("data/fig6c_autosomal_sites.csv")
freq$type <- factor(freq$type,
                    levels = c("SV deletions", "SV insertions", "SV other"))
plot6 <- ggplot2::ggplot(freq, ggplot2::aes(AF_assemblies, AF_PanGenie)) +
  ggplot2::stat_bin_hex(bins = 35) +
  ggplot2::facet_wrap(~type, nrow = 1) +
  ggplot2::scale_fill_gradientn(
    colours = c("white", "#fff4a1", "#ffc900", "#e13b16", "#160000"),
    trans = "log10", name = "Sites") +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(x = "AF assemblies", y = "AF PanGenie") +
  ggplot2::theme_classic()
ggplot2::ggsave("results/fig6c.png", plot = plot6,
                width = 11.25, height = 3.75, dpi = 200)
for (type_name in levels(freq$type)) {
  z <- freq[freq$type == type_name, ]
  cat(type_name, "sites", nrow(z), "Pearson r",
      round(cor(z$AF_assemblies, z$AF_PanGenie), 3), "\n")
}
