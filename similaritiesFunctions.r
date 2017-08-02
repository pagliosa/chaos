# ===============================
# About: Time series similarities methods
# Dependences: utils.r, rqa.r
# Author: Lucas Pagliosa
# Data of creation: 12/02/16
# Last revision 12/02/15
# ==============================

source("~/Canon/R/utils.r")
sourceFiles("~/Canon/R/RQA.r")

loadPackages("dtw")

dtwr <- function(ts1, ts2, ...)
{
  return(dtw(ts1, ts2)$distance)
}

dtwn <- function(ts1, ts2, ...)
{
  return(dtw(ts1, ts2)$normalizedDistance)
}

dtwd <- function(ts1, ts2, ...)
{
  return(dtwr(ts1, ts2) / euclidean(ts1, ts2))
}

distanceToDiagonal <- function(predict, test, size)
{
  linex = seq(1, size, length = length(predict))
  liney = seq(1, size, length = length(test))

  sum = sum(sqrt((predict - linex)^2 + (test - liney)^2))

  # printf("size/length: %d/%d\n", size, length(test))

  return (sum/length(test))
}

mddl <- function(predict, test, ..., plot = FALSE, shaded = FALSE, save = F,
  fileName = "mddl", xlab = "", ylab = "", retDTW = F)
{
  dtw = dtw(test, predict)

  mddlPlot <- function()
  {
    plot(dtw, xlab = xlab, ylab = ylab, ...)
    lines(x = seq(1, length(predict), length = length(predict)),
      y = seq(1,length(test),length = length(predict)), lty = 2)

    if(shaded)
    {
      x = seq(1, dtw$N, length = length(dtw$index1))
      y = seq(1, dtw$N, length = length(dtw$index2))

      for(i in 1:length(dtw$index1))
        lines(c(dtw$index1[i], x[i]), c(dtw$index2[i], y[i]))
    }
  }

  if (save)
    savePDF(mddlPlot(), fileName)
  if (plot)
    myPlot(mddlPlot())

  dist = distanceToDiagonal(dtw$index1, dtw$index2, length(test))

  if (!retDTW)
    return(dist)
  return(list(MDDL = dist, DTW = dtw$distance, DTWN = dtw$normalizedDistance))
}

distanceComparation <- function(nop = 200)
{
  sine = sin(2*-pi*seq(0, 2.5, len = nop))
  noise0 = sine + rnorm(nop, 0, 0.1)
  noise1 = sine + rnorm(nop, 0, 0.4)
  noise2 = sine + rnorm(nop, 0, 0.7)
  noise3 = sine + rnorm(nop, 0, 0.9)
  mean = rep(mean(sine), nop)

  plot(noise3, col = 5, type = "l")
  lines(noise2, col = 6)
  lines(noise1, col = 4)
  lines(noise0, col = 2)
  lines(sine, col = 3)
  lines(mean)

  e1 = euclidean(sine, noise0)
  e2 = euclidean(sine, noise1)
  e3 = euclidean(sine, noise2)
  e4 = euclidean(sine, noise3)
  e5 = euclidean(sine, mean)

  cat("ED(sine, noise0): ", e1, "\n")
  cat("ED(sine, noise1): ", e2, "\n")
  cat("ED(sine, noise2): ", e3, "\n")
  cat("ED(sine, noise3): ", e4, "\n")
  cat("ED(sine, mean): ", e5, "\n")

  cat("=====================\n")

  d1 = dtw(sine, noise0)$distance
  d2 = dtw(sine, noise1)$distance
  d3 = dtw(sine, noise2)$distance
  d4 = dtw(sine, noise3)$distance
  d5 = dtw(sine, mean)$distance

  cat("dtw(sine, noise0): ", d1,  "\n")
  cat("dtw(sine, noise1): ", d2, "\n")
  cat("dtw(sine, noise2): ", d3, "\n")
  cat("dtw(sine, noise3): ", d4, "\n")
  cat("dtw(sine, mean): ", d5, "\n")

  cat("=====================\n")

  cat("dtw-d(sine, noise0): ", d1/e1,  "\n")
  cat("dtw-d(sine, noise1): ", d2/e2, "\n")
  cat("dtw-d(sine, noise2): ", d3/e3, "\n")
  cat("dtw-d(sine, noise3): ", d4/e4, "\n")
  cat("dtw-d(sine, mean): ", d5/e5, "\n")

  cat("=====================\n")

  cat("dtwn(sine, noise0): ", dtwn(sine, noise0), "\n")
  cat("dtwn(sine, noise1): ", dtwn(sine, noise1), "\n")
  cat("dtwn(sine, noise2): ", dtwn(sine, noise2), "\n")
  cat("dtwn(sine, noise3): ", dtwn(sine, noise3), "\n")
  cat("dtwn(sine, mean): ", dtwn(sine, mean), "\n")

  cat("=====================\n")

  cat("mddl(sine, noise0): ", mddl(sine, noise0), "\n")
  cat("mddl(sine, noise1): ", mddl(sine, noise1), "\n")
  cat("mddl(sine, noise2): ", mddl(sine, noise2), "\n")
  cat("mddl(sine, noise3): ", mddl(sine, noise3), "\n")
  cat("mddl(sine, mean): ", mddl(sine, mean), "\n")

  cat("=====================\n")

  es = embedd(sine, 2, 1)
  e1 = embedd(noise0, 2, 1)
  e2 = embedd(noise1, 2, 1)
  e3 = embedd(noise2, 2, 1)
  e4 = embedd(noise3, 2, 1)
  e5 = embedd(mean, 2, 1)

  eis = eig(cov(es))
  ei1 = eig(cov(e1))
  ei2 = eig(cov(e2))
  ei3 = eig(cov(e3))
  ei4 = eig(cov(e4))
  ei5 = eig(cov(e5))

  cat("aed(sine, noise1): ", euclidean(eis, ei1), "\n")
  cat("aed(sine, noise2): ", euclidean(eis, ei2), "\n")
  cat("aed(sine, noise3): ", euclidean(eis, ei3), "\n")
  cat("aed(sine, noise3): ", euclidean(eis, ei4), "\n")
  cat("aed(sine, mean): ", euclidean(eis, ei5), "\n")

  cat("=====================\n")

  # Over embedding
  om = 6
  od = 1
  f = 0.2

  es = embedd(sine, om, od)
  e1 = embedd(noise0, om, od)
  e2 = embedd(noise1, om, od)
  e3 = embedd(noise2, om, od)
  e4 = embedd(noise3, om, od)
  e5 = embedd(mean, om, od)

  cat("rqa(sine, noise1): ", 1 / smwp(rqa(es, e1, f))$maxline, "\n")
  cat("rqa(sine, noise2): ", 1 / smwp(rqa(es, e2, f))$maxline, "\n")
  cat("rqa(sine, noise3): ", 1 / smwp(rqa(es, e3, f))$maxline, "\n")
  cat("rqa(sine, noise3): ", 1 / smwp(rqa(es, e4, f))$maxline, "\n")
  cat("rqa(sine, mean): ", 1 / smwp(rqa(es, e5, f))$maxline, "\n")
}

showMDDLIsGood <- function(nop = 500, plot = F, save = F, legend = F, numberOfTries = 30)
{
  ret = zeros(numberOfTries, 4)

  for (i in 1:numberOfTries)
  {
    sine = sin(2*-pi*seq(0, 2.5, len = nop))
    noise = sine + rnorm(nop, 0, 0.72)
    mean = rep(mean(sine), nop)

    mddlPlot <- function()
    {
      plot(noise, type = "l", xlab = "Time (t)", ylab = "Variable (x)", col = 2)
      lines(sine)
      lines(mean, col = 4)

      if (legend)
        legend(1, max(noise), c("Sine","noise 0.72", "mean"), lty = c(1, 4, 2));
    }

    if (save)
      savePDF(mddlPlot(), fileName = "MDDL-comparation")
    if (plot)
      myPlot(mddlPlot())

    edi = euclidean(sine, noise)
    ret[i, 1] = edi

    printf("ed(sine, noise): %f\n", edi)
    printf("ed(sine, mean): %f\n", euclidean(sine, mean))

    printf("==============\n")

    dtwni = dtwn(sine, noise)
    ret[i, 2] = dtwni

    printf("dtwn(sine, noise): %f\n", dtwni)
    printf("dtwn(sine, mean): %f\n", dtwn(sine, mean))

    printf("==============\n")

    dtwdi = dtwd(sine, noise)
    ret[i, 3] = dtwdi

    printf("dtwd(sine, noise): %f\n", dtwdi)
    printf("dtwd(sine, mean): %f\n", dtwd(sine, mean))

    printf("==============\n")

    mddl(sine, noise, plot = plot, xlab = "Sine", ylab = "Sine + N(0, 0.72)",
         save = save, fileName = "MDDL-match")
    mddl(sine, noise, plot = plot, xlab = "Sine", ylab = "Sine + N(0, 0.72)",
         shaded = T, save = save, fileName = "MDDL-diagonal")

    mddli = mddl(sine, noise)
    ret[i, 4] = mddli

    printf("mddl(sine, noise): %f\n", mddli)
    printf("mddl(sine, mean): %f\n", mddl(sine, mean))

    printf("==============\n")

  }

  perf = zeros(4, 3)

  for (i in 1:4)
  {
    vec = ret[,i]
    perf[i,] = c(min(vec), mean(vec), max(vec))
  }

  return(perf)
}
