library(mixtools)

turn_label_log <- function(x) {paste0("10e-", 6-x)}
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

#### AABB ----
DF1<- read.table("AABB_hist.txt")
wait <- log(DF1$V1+1,10)
wait[wait==0]=NA
wait<-wait[!is.na(wait)]
mixmdl <- normalmixEM(wait, k = 2)

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
  scale_x_continuous(labels = c(expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2)), limits=c(1, 4)) +
  ylab("Density") + xlab("Variant density (per bp)") +
  cowplot::theme_cowplot()
ggsave("AABB_all_density_distribution.pdf", width = 5, height = 2)


DF1<- read.table("DD_hist.txt")
wait <- log(DF1$V1+1,10)
wait[wait==0]=NA
wait<-wait[!is.na(wait)]
mixmdl <- normalmixEM(x=wait, k = 2)


p2 <- data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), bins=50, colour = "black", fill = "grey90", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
              args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
              colour = "purple", lwd = 1.5) +
  scale_x_continuous(labels = c(expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2)), limits=c(1, 4)) +
  ylab("Density") + xlab("Variant density (per bp)") +
  cowplot::theme_cowplot()
ggsave("DD_all_density_distribution.pdf", width = 5, height = 2)

cowplot::plot_grid(p1, p2, ncol=1, align = "v")
ggsave("all_density_distribution.pdf", width = 4.5, height = 3)
