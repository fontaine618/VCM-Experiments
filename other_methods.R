library(tidyverse)
library(magrittr)
# Generate data
set.seed(0)
pu = 1
n = 400 #*3
I = 60 #*3
noise_variance = 1.
re_sd = 0.5
f0 = function(t) 2/(1+exp((t-0.3)*10))
f1 = function(t) 2/(1+exp((0.8-t)*20))
s =  matrix(sample.int(I, n, T) - 1, n, 1)
t =  matrix(runif(n), n, 1)
t = round(t * 2^3) / 2^3
x = sample.int(2, I, T)-1
u = matrix(rnorm(n*pu), n, pu)
z = matrix(1., n, 1)
g = matrix(rnorm(I), I, 1)*re_sd
a = matrix(rnorm(pu), pu, 1)
t0 = seq(0, 1, 0.01)
b = matrix(
  c(f0(t0), f1(t0)),
  nrow=2, byrow=T
)
X = cbind(1, x[s+1])
bb = matrix(
  c(f0(t), f1(t)),
  ncol=2, byrow=F
)
m = rowSums(X*bb) + u %*% a + rowSums(z * matrix(g[s+1, ], ncol=1))
y = m + matrix(rnorm(n)*sqrt(noise_variance), n, 1)
data = data.frame(id=s, t=t, y=y, group=X[,2], u=u)


# refund FPCA # maybe this should be done in each group
fpca_df = data %>% select(id, t, y) %>% rename(.id=id, .index=t, .value=y)
fpca_results = refund::fpca.sc(
  ydata=fpca_df,
  nbasis=4
)
matplot(fpca_results$Yhat %>% t, type="l")

# commands
install.packages("spfda")
spfda::spfda()

install.packages("LocKer")
LocKer::LocKer()

install.packages("slasso")



xs = seq(-2, 2, 0.01)
ys = 1/xs * (1/(1+exp(-xs))-0.5)
