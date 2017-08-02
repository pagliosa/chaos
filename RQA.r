# ===============================
# About: Recurrence Quantification Analysis
# Dependences: utils.r, timeSeries.r
# Author: Yule Vaz, Lucas Pagliosa
# Last revision: 08/09/2016
# ===============================

source("~/Canon/R/utils.r")
sourceFiles("~/Canon/R/timeSeries.r")

loadPackages("pracma", "seewave", "gridGraphics")

heaviside <- function(x, a = 0)
{
	ret = x

	ret[x > a] = 1
	ret[x == a] = 1/2
	ret[x < a] = 0
	return(ret)
}

gammaCRQA <- function(R_i_j)
{
  penalty = 5
	if (R_i_j == 0)
	  penalty = 0.5
	return(penalty)
}

mydist <- function(M1,M2)
{
	nc1 = ncol(M1)
	nc2 = ncol(M2)
	D = zeros(nc1, nc2)

	for (i in 1:nc1)
		for (j in 1:nc2)
			D[i,j] = sqrt(sum((M1[,i] - M2[,j])^2))
	return(D)
}

rqa <- function(M1, M2, fraction = 0.1, filename)
{
	print("Calculating RQA")

  M1 = t(M1)
  M2 = t(M2)

	nc1 = ncol(M1)
	nc2 = ncol(M2)

	print("Calculating distances")
	distances = mydist(M1, M2)

	pos_i = ceil(nc2 * fraction)
	pos_j = ceil(nc1 * fraction)

	if (pos_i == 0)
		pos_i = 1;
	if (pos_j == 0)
		pos_j = 1;

	R = zeros(nc1, nc2)
	Q = zeros(nc1, nc2)

	eps_i = NULL

	print("Calculating epsilons")
	for (i in 1:nc1)
  {
		epsilons_i = sort(distances[i,], decreasing = FALSE)
		eps_i = rbind(eps_i, epsilons_i[pos_i])
	}

	eps_j = NULL

	for (j in 1:nc2)
  {
		epsilons_j = sort(distances[,j], decreasing = FALSE)
		eps_j = rbind(eps_j, epsilons_j[pos_j])
	}

	# plot(M1[1,], M1[2,])
	# points(M2[1,], M2[2,], pch = 2)
	# plot(1, type = "n", xlim = range(M1[1,]), ylim = range(M1[2,]))

	print("Producing diagonals")
	for (i in 3:nc1)
  {
	  # points(M1[1,i], M1[2,i], col = 2)
		for (j in 3:nc2)
	  {
		  # printf("%g * %g\n", eps_i[i] - distances[i, j], eps_j[j] - distances[i, j])
		  theta_i = heaviside(eps_i[i] - distances[i, j])
		  theta_j = heaviside(eps_j[j] - distances[i, j])

		  # points(M2[1,j], M2[2,j], pch = 2, col = 3)

      R_i_j = theta_i * theta_j

			R[i,j] = R_i_j

			if (R_i_j == 1)
				Q[i,j] = max(Q[i - 1, j - 1], Q[i - 2, j - 1], Q[i - 1, j - 2]) + 1
			else
		  {
				a = Q[i - 1, j - 1] - gammaCRQA(R[i - 1, j - 1])
				b = Q[i - 2, j - 1] - gammaCRQA(R[i - 2, j - 1])
				c = Q[i - 1, j - 2] - gammaCRQA(R[i - 1, j - 2])
				Q[i, j] = max(0, a, b, c)
			}
		}
	}

	maxQ = max(Q) + 3
	list = list()
	list$R = R
	list$Q = Q
	list$maxline = round(maxQ)
	list$dist = 1 / maxQ

	#image(Q)

  # printMat(R, "R")
  # printMat(Q, "Q")
	# printf("mQ: %g\n", max(Q))

	# write(R, sprintf("%s-R.txt", filename), append = F, ncol = nc2, sep = "\t")
	# write(Q, sprintf("%s-Q.txt", filename), append = F, ncol = nc2, sep = "\t")

	return(list)
}

testRQA <- function(nop = 1000, fraction = 0.1)
{
  det = sin(2 * pi * seq(0, 9, length = nop))
  series = det + rnorm(mean = 0, sd =  0.1, n = nop)

  data1 = embedd(det, m = 2, d = 20)
  data2 = embedd(series, m = 2, d = 20)

  plot(data1)
  points(data2, col = 2)

  return(rqa(data1, data2, fraction = fraction))
}

testRQA2 <- function(nop = 1000, fraction = 0.1, centralize = T)
{
  det = createLogistic(nop)
  series = det + 4

  data1 = embedd(det, m = 2, d = 1)
  data2 = embedd(series, m = 2, d = 1)

  if (centralize)
    data1 = centralizeEmbddings(data1, data2)

  plotInterval(list(data1, data2))
  points(data1[,1], data1[,2])
  points(data2[,1], data2[,2], col = 2)

  rqa = rqa(data1, data2, fraction = fraction)
  image(rqa$Q)
  return(rqa)
}

testRQA3 <- function(nop = 200)
{
  s0 = createLogistic(plot = T, r = 3.8, nop = nop, noise = 0.00)
  s1 = createLogistic(plot = T, r = 3.8, nop = nop, noise = 0.02)
  s2 = createLogistic(plot = T, r = 3.9, nop = nop, noise = 0.08)
  s3 = createLogistic(plot = T, r = 3.3, nop = nop, noise = 0.05)
  s4 = createLogistic(plot = T, r = 3.3, nop = nop, noise = 0.00)
  s5 = noise = runif(nop, 0, 1);

  e0 = embedd(s0, m = 2, d = 1)
  e1 = embedd(s1, m = 2, d = 1)
  e2 = embedd(s2, m = 2, d = 1)
  e3 = embedd(s3, m = 2, d = 1)
  e4 = embedd(s4, m = 2, d = 1)
  e5 = embedd(s5, m = 2, d = 1)

  # centralize(list(e0, e1, e2, e3, e4, e5))

  myPlot({
    plotInterval(list(e0, e1, e2, e3, e4), 0)
    points(e0[,1], e0[,2])
    points(e1[,1], e1[,2], pch = 2, col = 2) # vermelho
    points(e2[,1], e2[,2], pch = 3, col = 3) # verde
    points(e3[,1], e3[,2], pch = 4, col = 4) # azul forte
    points(e4[,1], e4[,2], pch = 19, col = 5, cex = 3) # azul fraco
    points(e5[,1], e5[,2], pch = 6, col = 6) # magenta
  })

  for (i in seq(0.1, 0.5, 0.1))
  {
    printf("%f) 3.8 - 3.8 + N(0, 0.02): %d\n", i, smwp(rqa(e0, e1, i)$maxline))
    printf("%f) 3.8 - 3.9 + N(0, 0.08): %d\n", i, smwp(rqa(e0, e2, i)$maxline))
    printf("%f) 3.8 - 3.3 + N(0, 0.05): %d\n", i, smwp(rqa(e0, e3, i)$maxline))
    printf("%f) 3.8 - 3.3 + N(0, 0.00): %d\n", i, smwp(rqa(e0, e4, i)$maxline))
    printf("%f) 3.8 - U(0, 1): %d\n\n", i, smwp(rqa(e0, e5, i)$maxline))

    printf("%f) 3.9 + N(0, 0.08) - 3.8: %d\n", i, smwp(rqa(e2, e0, i)$maxline))
    printf("%f) 3.9 + N(0, 0.08) - 3.8 + N(0, 0.02): %d\n", i, smwp(rqa(e2, e1, i)$maxline))
    printf("%f) 3.9 + N(0, 0.08) - 3.3 + N(0, 0.05): %d\n", i, smwp(rqa(e2, e3, i)$maxline))
    printf("%f) 3.9 + N(0, 0.08) - 3.3 + N(0, 0.00): %d\n", i, smwp(rqa(e2, e4, i)$maxline))
    printf("%f) 3.9 + N(0, 0.08) - U(0, 1): %d\n\n", i, smwp(rqa(e2, e5, i)$maxline))
  }
}

testRQA4 <- function(nop = 1000, fraction = 0.1, centralize = T)
{
  ts1 = createSunspot(nop = nop, nos = nop)
  ts2 = createRossler(nop)

  data1 = embedd(ts1, m = 2, d = 8)
  data2 = embedd(ts2, m = 2, d = 8)

  if (centralize)
    data1 = centralizeEmbddings(data1, data2)

  plotInterval(list(data1, data2))
  lines(data1[,1], data1[,2])
  lines(data2[,1], data2[,2], col = 2)

  rqa = rqa(data1, data2, fraction = fraction)
  image(rqa$Q)
  return(rqa)
}
# Used in article
compareRQA <- function(nop = 500, m = 2, d = 1, fraction = 0.1, numberOfTries = 30)
{
  clearPlots()

  # ED, DTW, DTW-D, MDDL, RQA (x2)
  ret = zeros(numberOfTries, 10)

  for (i in 1:numberOfTries)
  {
    # l1 = createLogistic(plot = T, r = 3.8, nop = nop, noise = 0.01)
    # l2 = createLogistic(plot = T, r = 3.6, nop = nop, noise = 0.01)
    # noise = (rnorm(nop, mean = 0.6, sd = 0.01)); plot(noise, type = "l")
    # noise = runif(nop, 0, 1); plot(noise, type = "l")

    twiceNop = 2 * nop
    l1 = sin(2 * pi * seq(0, 9, len = nop))
    l2 = cos(2 * pi * seq(0, 9, len = twiceNop)) + runif(min = -0.1, max = 0.1,
      n = twiceNop)
    l2 = blackman.w(n = twiceNop) * l2
    l2 = l2[(nop + 1):twiceNop]
    noise  = runif(min = -0.5, max = 0.5, n = nop)

    col1 = 1
    col2 = 4
    coln = 2

    x = c(l1, l2, noise)

    l2_str = nop + 1
    l2_end = 2 * nop

    noise_str = l2_end + 1
    noise_end = l2_end + nop

    e1 = embedd(x[1:nop], m, d)
    e2 = embedd(x[l2_str:l2_end], m, d)
    e3 = embedd(x[noise_str:noise_end], m, d)

    savePDF(function()
    {
      plot(x, type = "l")
      lines(1:nop, x[1:nop], col = col1)
      lines(l2_str:l2_end, x[l2_str:l2_end], col = col2)
      lines(noise_str:noise_end, x[noise_str:noise_end], col = coln)
    }, "Series")

    # myPlot(plot(x[1:nop], type = "l", col = col1), save, "Logistic 1")
    # myPlot(plot(x[l2_str:l2_end], type = "l", col = col2), save, "Logistic 2")
    # myPlot(plot(x[noise_str:noise_end], type = "l", col = coln), save, "Noise")

    savePDF(function()
    {
      plotInterval(list(e1, e2, e3))
      points(e1[,1], e1[,2], col = col1, pch = 16)
      points(e2[,1], e2[,2], col = col2, pch = 2)
      points(e3[,1], e3[,2], col = coln, pch = 3)
    }, "Embeddings")

    ed1 = euclidean(l1, l2) / length(l1)
    ed2 = euclidean(l1, noise) / length(l1)
    ret[i, 1] = ed1
    ret[i, 2] = ed2

    printf("%d/%d\n", i, numberOfTries)

    printf("ED(logisti1, logistic2): %g\n", ed1)
    printf("ED(logistic1, noise): %g\n", ed2)
    printf("ED(logistic2, noise): %g\n", euclidean(l2, noise) / length(l1))

    printf("\n")

    dtw1 = dtwn(l1, l2)
    dtw2 = dtwn(l1, noise)
    ret[i, 3] = dtw1
    ret[i, 4] = dtw2

    printf("DTW(logisti1, logistic2): %g\n", dtw1)
    printf("DTW(logistic1, noise): %g\n", dtw2)
    printf("DTW(logistic2, noise): %g\n", dtwn(l2, noise))

    printf("\n")

    dtwd1 = dtwd(l1, l2)
    dtwd2 = dtwd(l1, noise)
    ret[i, 5] = dtwd1
    ret[i, 6] = dtwd2

    printf("DTW-D(logisti1, logistic2): %g\n", dtwd1)
    printf("DTW-D(logistic1, noise): %g\n", dtwd2)
    printf("DTW-D(logistic2, noise): %g\n", dtwd(l2, noise))

    printf("\n")

    mddl1 = mddl(l1, l2)
    mddl2 = mddl(l1, noise)
    ret[i, 7] = mddl1
    ret[i, 8] = mddl2

    printf("MDDL(logisti1, logistic2): %g\n", mddl1)
    printf("MDDL(logistic1, noise): %g\n", mddl2)
    printf("MDDL(logistic2, noise): %g\n", mddl(l2, noise))

    printf("\n")

    e1 = embedd(x[1:nop], m, d)
    e2 = embedd(x[l2_str:l2_end], m, d)
    e3 = embedd(x[noise_str:noise_end], m, d)

    rqa = rqa(e1, e2, fraction, "RQA1")
    rqa1 = rqa$dist
    rqa2 = (rqa(e1, e3, fraction, "RQA2"))$dist
    ret[i, 9] = rqa1
    ret[i, 10] = rqa2

    printf("RQA(logisti1, logistic2): %g\n", rqa1)
    printf("RQA(logistic1, noise): %g\n", rqa2)
    printf("RQA(logistic2, noise): %g\n", smwp(rqa(e2, e3, fraction,"RQA3"))$dist)
  }

  perf = zeros(10, 3)

  for (i in 1:10)
  {
    vec = ret[,i]
    perf[i,] = c(mean(vec), min(vec), max(vec))
  }

  write.table(perf, file = "seewave.txt", row.names = F, col.names = F)

  clearPlots()

  return(perf)

}
