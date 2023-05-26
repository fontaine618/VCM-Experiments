library(tidyverse)
library(magrittr)
theme_set(theme_minimal())

xs = seq(0, 1, 0.001)

K = 20

X = splines::bs(xs, df=K)

beta = c(rep(1, 3), rep(0, 2), rep(1, 5), rep(0, 5), rep(1, 5))


y = X %*% beta

df = data.frame(basis=X, spline=y, x=xs)

df_long = df %>% tidyr::pivot_longer(
  cols=-x
)

df_coef = data.frame(
  x=apply(X, 2, function(x) xs[which.max(x)]),
  y=beta
)


g0 = ggplot(
  data=df_long,
  mapping=aes(x=x, y=value, group=name, color=name)
) + geom_line() + 
  scale_color_manual(values=c("spline"="black")) + 
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  coord_cartesian(xlim=c(-0.01, 1.01))

g1 = ggplot(
  data=df_coef,
  mapping=aes(x=x, y=pmax(y, 0.01))
) + geom_bar(stat="identity") + coord_cartesian(xlim=c(-0.01, 1.01)) +
   ylab("Coef")


g = cowplot::plot_grid(g0, g1, nrow=2, ncol=1, align="hv", axis="lr", rel_heights=c(2, 1))
ggsave("./figures/splines.pdf", g, width=8, height=6)
