# ===============================
# About: Chaos functions
# Dependences: none
# Author: Lucas Pagliosa
# Last revision: 04/07/17
# ==============================

source("~/Canon/R/utils.r")
loadPackages("forecast")
cls()

myACF <- function(ts, tau = 2)
{
  range = 1:(length(ts) - tau)
  s = ts[range]
  d = ts[range + tau]
  m = mean(ts)

  return(mean((s - m) * (d - m)) / var(ts))
}

findLE <- function(ts = createCombinedSine(nop = 2e3))
{
	output = lyap_k(series = ts, m = 2, d = 1, s = 200, t = 40,
	  ref = 1700, k = 2, eps = 4)
	return(lyap(output, 0.73, 2.47))
}

findMaxLyapunovExponents <- function()
{
  my.ts = createLorenz()
  ml=maxLyapunov(time.series=my.ts,
  min.embedding.dim=3,
  max.embedding.dim=10,
  time.lag=1,
  radius=50,theiler.window=100,
  min.neighs=2,min.ref.points=500,
  max.time.steps=40,do.plot=FALSE)
  plot(ml)
  ml.estimation = estimate(ml,regression.range = c(0,30),
  use.embeddings=4:5,
  do.plot = TRUE)
  # The max Lyapunov exponent of the Henon system is 0.41
  cat("expected: ",0.41," calculated: ",ml.estimation,"\n")
}