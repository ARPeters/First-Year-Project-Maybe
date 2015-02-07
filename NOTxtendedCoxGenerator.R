install.packages("psych.r")
install.packages("survsim.r")

?rnorm()

install.packages("msm.R")
require(msm)

require(survival)

# CREATING g() AND g^-1()
g.inv <- sqrt
g <- function(x) {
  x^2
}
# CREATING THE TIME SCALE AND TRANSFORMED TIME SCALE
t <- 0:199
t.diff <- (t[-1] - t[1:(length(t) - 1)])[-(length(t) - 1)]
g.inv.t <- g.inv(t)
g.inv.t.diff <- (g.inv(t[-1]) - g.inv(t[1:(length(t) - 1)]))[-(length(t) - 1)]

#CREATING THE BOUNDS OF TRUNCATION
t.max <- 150
t.min <- 10
g.inv.t.max <- g.inv(t.max)
g.inv.t.min <- g.inv(t.min)
#DATA GENERATING PROCESS FOR COVARIATE
B <- function(N, m, M) {
  runif(N, m, M)
}
#BETA
b <- 2
#NUMBER OF OBSERVATIONS
n <- 1000
#CREATING DATA VECTOR
z.list <- list()
for (i in 1:n) {
  z <- B(length(t), -0.5, 0.5)
  z.list[[i]] <- cbind(z, exp(b * z))
}
#GENERATING DATA USING ACCEPT-REJECT METHOD
k <- function(x, m, M, rates, t){
  ifelse(x <= m | x >= M, 0, dpexp(x, rates, t))
}
gen.y <- function(x) {
  x1 <- x[, 2]
  d <- ppexp(g.inv.t.max, x1, g.inv.t) - ppexp(g.inv.t.min, x1, g.inv.t)
  M <- 1 / d
  r <- 60
  repeat{
    y <- rpexp(r, x1, g.inv.t)
    u <- runif(r)
    t <- M * ((k(y, g.inv.t.min, g.inv.t.max, x1, g.inv.t) / d /
                 dpexp(y, x1, g.inv.t)))
    y <- y[u <= t][1]
    if (!is.na(y)) break
  }
  y
}
y <- sapply(z.list, gen.y)
g.y <- g(y)


#CREATING CENSORING INDICATOR
prop.cen <- 0.5
d <- sample(0:1, n, replace = TRUE, prob = c(prop.cen, 1 - prop.cen))
23#CREATING DATASET
data <- NULL
for (i in 1:n) {
  id.temp <- rep(i, ceiling(g.y[i]))
  time.temp <- c(1:ceiling(g.y[i]))
  time0.temp <- 0:ceiling(g.y[i] - 1)
  d.temp <- c(rep(0, length(time.temp) - 1), d[i])
  z.temp <- z.list[[i]][1:(ceiling(g.y[i])), 1]
  data.temp <- cbind(id.temp, time.temp, time0.temp, d.temp, z.temp)
  data <- rbind(data, data.temp)
}
colnames(data) <- c(’id’, ’t’, ’t0’, ’d’, ’z1’)
data <- data.frame(data)
model <- coxph(Surv(t0, t, d) ~ z1, data = data)
schoenfeld <- cox.zph(model, transform = ’identity’)
#RESULT
data
summary(model)
schoenfeld