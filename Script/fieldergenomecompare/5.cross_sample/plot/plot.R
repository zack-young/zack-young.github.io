library(ggplot2)

#### part1 ----

my_confidence_intervals <- function(DF){
  sample_mean <- mean(DF)
  standard_dev <- sd(DF)
  mean_minus_margin_of_error <- sample_mean - standard_dev
  mean_plus_margin_of_error <- sample_mean + standard_dev
  c(mean_minus_margin_of_error, sample_mean, mean_plus_margin_of_error)
}

turn_label <- function(x) {x/10000}
turn_label_log <- function(x) {paste0("10e-", 6-x/10)}

#toplot <- c("AABB_all", "AABB_HH", "AABB_HT", "AABB_TT", "AA_all", "AA_HH", "AA_HT", "AA_TT", "BB_all", "BB_HH", "BB_HT", "BB_TT", "D_all", "D_HH", "D_HD", "D_DD")
toplot <- c("AB", "D")
tmplist <- list()
for (j in 1:length(toplot)){
  CHR <- toplot[j]
  DF2 <- read.table(paste0("/data/user/yangzz/mapping/fieldergenomecompare/1.diff_dev/201_final/", CHR, "_sample_by_count.txt"))

  tmp <- data.frame()
  for (i in seq(1,100)){
    t <- my_confidence_intervals(DF2[,i])
    tmp <- rbind(tmp, c(i, t))
  }
  colnames(tmp) <- c("density", "lower", "mean", "upper")
  tmp[tmp<0] <- 0
  p <- ggplot(tmp)+
    geom_ribbon(aes(density, ymin=lower,ymax=upper), fill="grey80", col=NA) +
    # geom_errorbar(aes(density, ymin = lower, ymax = upper)) +
    geom_line(aes(density, mean), size=1.2, col="red") +
    geom_vline(xintercept = 30) +
    scale_x_continuous(limits=c(0,40), labels = c(expression(10^-6), expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2))) +
    cowplot::theme_cowplot() +
    ylab(NULL) + xlab(NULL)
  # ggsave(paste0(CHR, ".logscale.pdf"), p, height = 2.5, width = 4)
  tmplist[[j]] <- p
}
pall <- cowplot::plot_grid(plotlist = tmplist)
ggsave("all.logscale.pdf", pall, height = 7, width = 10)
