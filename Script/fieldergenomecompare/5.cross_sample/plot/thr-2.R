library(ggplot2)

#### part1 ----

turn_label_log <- function(x) {paste0("10e-", 6-x)}
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

toplot <- c("AABB_all", "D_all")
tmplist <- list()
for (j in 1:length(toplot)){
  CHR <- toplot[j]
  DF2 <- read.table(paste0("../log_scale_dist/", CHR, "_sample_by_count.txt"))
  
  DF3 <- colSums(DF2)
  cumu <- c()
  for (i in 1:100){
    cumu <- c(cumu, rep(i, DF3[i]))
  }
  mixmdl <- normalmixEM(x=cumu, k = 2)

  p1 <- data.frame(x = mixmdl$x) %>%
    ggplot() +
    geom_histogram(aes(x, ..density..), bins=50, colour = "black", fill = "grey90", 
                   fill = "white") +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                  colour = "red", lwd = 1.5) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                  colour = "blue", lwd = 1.5) +
    scale_x_continuous(limits=c(1, 40)) +
    ylab("Density") + xlab("Variant density (per bp)") +
    cowplot::theme_cowplot()
  ggsave(paste0(CHR, "_norm.pdf"), width = 5, height = 2)
}

