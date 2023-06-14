library(VCMM)
library(ggplot2)
theme_set(theme_minimal())

# plot masking
ts = seq(0, 1, 0.1)
h = 0.25
rho = 10
sig2 = 1


k = function(x, y) exp(-(outer(x, y, "-") / h)^2)/h
# k = function(x, y) 1*(abs(outer(x, y, "-")) < h) /h
# k = function(x, y) 0.75*pmax(1-(outer(x, y, "-")/h)^2, 0)/h

t0 = 0.4
wit = sqrt(k(ts, t0))

Vi = (diag(length(ts)) + rho ) *sig2
Pi = solve(Vi)

to_long = function(mat) reshape2::melt(mat, c("x", "y"), value.name = "z")

plot_mat = function(mat) ggplot() + 
  geom_tile(
    data=to_long(mat),
    mapping=aes(x=x, y=y, fill=z)
  ) + xlab("") + ylab("") + 
  labs(fill="Z") + 
  scale_fill_gradient2() + 
  scale_y_reverse() + 
  theme(
    legend.position="none",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) + coord_fixed() + 
  geom_text(
    data=to_long(mat),
    mapping=aes(x=x, y=y, label=round(z, 1)),
    size=2
    )

PWi = Pi * (wit %*% t(wit))
witinv = ifelse(wit < 1e-4, 0, 1/wit)
Vit = Vi * (witinv %*% t(witinv))
# Vit = MASS::ginv(PWi)
dVit = sqrt(diag(Vit))
Cit = Vit / (dVit %*% t(dVit))

gw = plot_mat(wit %*% t(wit)) + ggtitle(latex2exp::TeX("$\\textbf{W}_i(t)$"))
gp = plot_mat(Pi) + ggtitle(latex2exp::TeX("$\\widetilde{\\textbf{P}}_i$"))
gpw = plot_mat(PWi) + ggtitle(latex2exp::TeX("$\\widetilde{\\textbf{P}}_i(t)$"))
gvw = plot_mat(Cit) + ggtitle(latex2exp::TeX("$\\widetilde{\\textbf{V}}_i(t)$"))



odot = ggplot() + 
  annotate("text", label="\u2299", x=1, y=1, size=10) + 
  theme_void()
eq = ggplot() + 
  annotate("text", label="=", x=1, y=1, size=10) + 
  theme_void()

cowplot::plot_grid(gp, odot, gw, eq, gpw, nrow=1, rel_widths = c(3, 1, 3, 1, 3))
cowplot::plot_grid(gp, odot, gw, eq, gpw, gvw, nrow=1, rel_widths = c(3, 1, 3, 1, 3, 3))

round(Vit %*% PWi, 2)









GP = diag(length(ts)) * 0.
GW = GP
for(t0 in ts){
  wit = sqrt(k(ts, t0))
  GW = GW + (wit %*% t(wit))
  PWi = Pi * (wit %*% t(wit))
  GP = GP + PWi
}
GV = solve(GP)
Gsd = sqrt(diag(GV))
GC = diag(1/Gsd) %*% GV %*% diag(1/Gsd)

gGV = plot_mat(GV) + ggtitle("Global variance")
gGC = plot_mat(GC) + ggtitle("Global correlation")
gGP = plot_mat(GP) + ggtitle("Global precision")
gGW = plot_mat(GW) + ggtitle("Total weight")



cowplot::plot_grid(gGW, gGP, gGV, gGC, nrow=1)
GC
