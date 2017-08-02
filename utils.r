# ===============================
# About: Utilities functions
# Dependences: none
# Author: Lucas Pagliosa
# Last revision: 29/01/17
# ==============================

canonPath <- function()
{
  return ("~/Canon/R")
}

reset <- function(dir = "~/Canon/R", env = globalenv())
{
  clearPlots()
  rm(list = ls(all = TRUE, envir = env), envir = env)
  source(sprintf("%s/utils.r", dir))
  loadGraphics()
  loadPackages("FNN")
  sourceDir()
  sourceDir(dir)
  cls()
}

cls <- function()
{
  cat("\014")
}

allFuncs <- function()
{
  return(ls(all = TRUE, envir = globalenv()))
}

# declared function:
# f <- function(...)
# myPlot(f())
# ------------------
# implicit function
# myPlot(function(){...})

myPlot <- function(plotFunc, savePlot = F, filename, newPlot = F)
{
  if (newPlot)
    dev.new()
  par(mar = c(4.4, 5, 1, 0.4) + 0.1, cex.lab = 2.5, cex.axis = 1.8)
  plotFunc()
  if (savePlot)
  {
    savePDF(plotFunc, filename)
    dev.set(dev.prev())
  }
}

plotInterval <- function(listOfSpaces, margin = 0.2)
{
  nod = ncol(listOfSpaces[[1]])
  ranges = repmat(c(Inf, -Inf), nod, 1)
  margin = c(-margin, margin)

  for (space in listOfSpaces)
  {
    for (d in 1:nod)
    {
      r = range(space[,d])
      if (ranges[d, 1] > r[1])
        ranges[d, 1] = r[1]
      if (ranges[d, 2] < r[2])
        ranges[d, 2] = r[2]
    }
  }
  plot(1, type = "n", xlim = ranges[1,] + margin, ylim = ranges[2,] + margin)
}

savePDF <- function(plotFunc, filename = "R-image",
  pathDir = "C:/Users/pagliosa/Desktop", width = 7, height = 5)
{
  path = sprintf("%s/%s.pdf", pathDir, filename)
  pdf(file = path, width = width, height = height)
  par(mar = c(4.4, 5, 1, 0.4) + 0.1, cex.lab = 2.5, cex.axis = 1.8)
  plotFunc()
  dev.off()
}

getFiles <- function(path = "./", fileFormat = "[.][Rr]$", recursive = T)
{
  fileData = list()
  dirs = list.dirs(path = path, full.names = TRUE, recursive = recursive)
  count = 1

  for (dir in dirs)
  {
    printf("Searching dir: %s\n", dir)

    files = NULL
    allFiles = list.files(dir, pattern = fileFormat, all.files = FALSE)
    nof = length(allFiles)

    for (rFiles in allFiles)
    {
      files = c(files, file.path(dir, rFiles))
      printf("> File found: %s\n", rFiles)
    }

    if (nof > 0)
    {
      fileData[[count]] = list(files = files, name = getPathFile(dir),
        path = dir, nof = length(allFiles))
      count = count + 1
    }
  }

  return(fileData)
}

sourceFiles <- function(...)
{
  files = c(...);

  for (i in 1:length(files))
    source(files[i])
}

sourceDir <- function(path = "./")
{
  for (rFiles in list.files(path, pattern = "[.][Rr]$"))
  {
    file = file.path(path, rFiles)
    cat("Loading ", file, "\n")
    source(file)
  }
}

clearPlots <- function()
{
  graphics.off()
  par(mfrow = c(1, 1))
}

loadPackages <- function(...)
{
  packages = c(...)

  for (p in packages)
  {
    if (p %in% rownames(installed.packages()) == FALSE)
      suppressWarnings(suppressMessages(install.packages(p)))
    suppressWarnings(suppressMessages(require(p, character.only = TRUE)))
  }
}

loadGraphics <- function()
{
  loadPackages("plotrix", "scatterplot3d", "gplots", "lattice", "rgl")
}

error <- function(message)
{
  message(sprintf("ERROR: %s", message))
  return(0)
}

printfCounter = 0

openHeader <- function(...)
{
  printf(...)
  printfCounter <<- printfCounter + 1
}

closeHeader <- function(...)
{
  pfc = printfCounter - 1
  if (pfc >= 0)
    printfCounter <<- pfc
  printf(...)
}

printfInline <- function(...)
{
  invisible(cat(sprintf(...)))
}

printf <- function(...)
{
  for (i in seq(1, 1, length.out = printfCounter))
    cat(">")
  if (printfCounter > 0)
    cat(" ")
  printfInline(...)
}

highlightPrintf <- function(..., borderSpace = 0)
{
  for (i in seq(1, 1, length.out = borderSpace))
    printf("\n")
  printf("================\n")
  printf(...)
  printf("\n")
  printf("================\n")
  for (i in seq(1, 1, length.out = borderSpace))
    printf("\n")
}

printMat <- function(mat, message, onlySize = F)
{
  if (is.vector(mat))
    mat = matrix(mat, nrow = 1)

  if (missing(message))
    message = "mat"

  printf("%s\n(%d, %d)\n", message, nrow(mat), ncol(mat))

  if (onlySize)
    return()

  for (i in 1:nrow(mat))
  {
    for (j in 1:ncol(mat))
      printf("%g ", mat[i, j])
    printf("\n")
  }
}

addList <- function(list, newItem)
{
  list[[length(list) + 1]] = newItem
  return(list)
}

# Suppres messages, warnings and prints
smwp <- function(func, debug = F)
{
  sw <- function(func) suppressMessages(suppressWarnings(func))

  if (debug)
    res = sw(func)
  else
    sw(capture.output({
      res = func;
    }))

  return(res)
}

getOnlyResultsThatWorked <- function()
{
  set.seed(123)
  x <- stats::rnorm(50)
  doit <- function(x)
  {
    x <- sample(x, replace = TRUE)
    if(length(unique(x)) > 30) mean(x)
    else stop("too few unique points")
  }
  ## alternative 1
  res = lapply(1:100, function(i)
  {
    cat(i)
    try(doit(x), TRUE)
  })
  return(res)
}

numDigits <- function(data, number = 4)
{
  return(format(round(data, number), nsmall = number))
}

bnorm <- function(p, mu = 0, sd = 1)
{
  pn = pnorm(p, mu, sd)

  return(2 * pn - 1)
}

normalize <- function(min, max, VALUE, MIN, MAX)
{
  # value - min       VALUE - MIN
  # _____________ =  _____________
  #  max - min         MAX - MIN
  return((VALUE - MIN) * (max - min) / (MAX - MIN) + min)
}

normalizeVec <- function(min, max, vec)
{
  MIN = min(vec)
  MAX = max(vec)
  
  norm = c()
  for (i in 1:length(vec))
    norm = c(norm, normalize(min, max, vec[i], MIN, MAX))
  return(norm)
}

sqtDist <- function(a, b)
{
  return(sum((a - b)^2))
}

euclidean <- function(a, b, ...)
{
  return(sqrt(sqtDist(a, b)))
}

init <- function(n, m, value, min = 0, max = 1, byrow = T)
{
  if (missing(n))
    n = 1
  if (missing(m))
    m = 1

  # Fill matrix with increasing order
  if (missing(value))
    return(matrix(1:(n * m), n, m, byrow = byrow))

  # Fill matrix with random numbers
  if (is.logical(value) && !value)
    return(matrix(runif(n * m, min, max), n, m))

  # Fill matrix with repeated value
  return(repmat(value, n, m))
}

contains <- function(v1, v2)
{
  for (i in 1:length(v1))
    if (!any(v2 == v1[i]))
      return(FALSE)
  return(TRUE)
}

computeStatistics <- function(vec)
{
  n = length(vec)

  # Statistical moments
  mean = mean(vec)
  variance = var(vec)

  # Shannon
  sum = sum(vec)
  if (sum == 0)
    shannon = 0
  else
  {
    probs = vec / sum
    shannon = -sum(probs %*% log2(abs(probs) + 1))
  }

  # Energy
  energy = sum(abs(vec))/n

  return(list("mean" = mean, "var" = variance, "shannon" = shannon,
    "energy" = energy))
}

zeros <- function(n, m = 1)
{
  ret = matrix(0, n, m)
  if (missing(m))
    ret = as.vector(ret)
  return(ret)
}

lastOccurence <- function(string, char)
{
  indices = gregexpr(char, string)[[1]]
  lastOccurence = indices[length(indices)]
  return(lastOccurence)
}

getPathFile <- function(path)
{
  lastSlide = lastOccurence(path, "/") + 1
  file = substr(path, lastSlide, nchar(path))
  return(file)
}

globalClock = proc.time()

startClock <- function()
{
  # <<- sets global variables
  clock = proc.time()
  globalClock <<- clock
  return(clock)
}

endClock <- function(clock)
{
  if (missing(clock))
    clock = globalClock
  diff = proc.time() - clock
  return(diff)
}

getOS <- function()
{
  return(as.vector(Sys.info()["sysname"]))
}

fatorial <- function(arg)
{
  fat = 1
  for (i in 1:arg)
    fat = fat * i
  return(fat)
}

centralizeEmbddings <- function(emb1, emb2)
{
  for (i in 1:ncol(emb1))
    emb1[,i] = emb1[,i] + mean(emb2[,i]) - mean(emb1[,i])
  return(emb1)
}

computeSlope <- function(xvec, yvec, plot = F)
{
  lm = lm(yvec ~ xvec)
  y = as.numeric(lm$coefficients[1])
  slope = as.numeric(lm$coefficients[2])
  
  if (plot)
    abline(y, slope)
  
  return(list("y" = y, "slope" = slope))
}

computeKnnData <- function(space, k)
{
  knn = get.knn(space, k, algorithm = "kd_tree")
  return(list("id" = knn$nn.index, "dist" = knn$nn.dist))
}

differentiate <- function(vec, h = 1e-3)
{
  len = length(vec)
  derivative = (vec[3:len] - vec[1:(len - 2)])/(2 * h)
  return(derivative)
}
