# ===============================
# About: Dummy Performance Test for R Parallel Methods
# Dependences: utils.r
# Author: Lucas Pagliosa
# Last revision: 29/01/17
# ==============================

source("~/Canon/R/utils.r")

loadPackages("foreach", "doParallel", "parallel")

numberOfCores <- function()
{
  noc = detectCores()
  return(noc)
}

startParallel <- function(noc = -1, logFile)
{
  openHeader("Starting parallelism...\n")

  # Get number of cores
  realCores = numberOfCores()

  if (noc == -1)
    noc = realCores

  printf("Allocating %d/%d (%.2f%%) of cores\n", noc, realCores,
    noc / realCores * 100)

  # Create clusters
  clusters =  makeCluster(noc)

  # Export variables and functions
  clusterExport(clusters, allFuncs())

  # Export libraries
  clusterEvalQ(clusters, c())

  # Register clusters (needed for foreach)
  registerDoParallel(clusters)

  # Start parallel log
  startParLog()

  return(clusters)
}

endParallel <- function(clusters = cl)
{
  # Stop clusters
  stopCluster(clusters)
  closeHeader("Ending parallelism\n\n")
}

parLog = "Debug/parLog.txt"

saveInLog <- function(message, log = parLog, append = F)
{
  sink(log, append)
  message
  sink()
}

startParLog <- function(log = parLog)
{
  writeLines(c(""), log)
}

saveInLog <- function(message, log = parLog, append = F)
{
  sink(log, append)
  message
  sink()
}

printElapsedTime <- function(clock)
{
  time = endClock(clock)
  printf("Time elapsed: %.3f seconds\n", time["elapsed"])
}

cl = NULL

# cc stands for closeCluster
run <- function(func)
{
  clock = startClock()
  result = func
  printElapsedTime(clock)
  return(result)
}

runPar <- function(func, noc)
{
  cl <<- startParallel(noc)
  ret = run(func)
  endParallel(cl)
  return(ret)
}

runSequencial <- function(func = toy, args = c(1:5e3), method = "lapply")
{
  openHeader("Running sequencial '%s'\n", method)

  # ====================
  # Sequencial methods
  # ====================

  # Return a list
  if (method == "lapply")
    ret = run({lapply(args, func)})

  # Return a vector
  else if (method == "sapply")
    ret = run({sapply(args, func)})

  # For version of lapply
  else if (method == "foreach")
  {
    ffunc <- function(func, args)
    {
      printf("Executing foreach...")
      data = foreach(fat = args, .packages = c(), .combine = rbind) %do%
          func(fat)
      printfInline("Done\n")
      return(data)
    }
    ret = run(ffunc(func, args))
  }
  else
    error("Invalid method")

  closeHeader("Done\n\n")
  return(ret)
}

runParallel <- function(func = toy, args = c(1:10), method = "lapply", noc = -1)
{
  printf("Running parallel '%s'\n", method)

  # ====================
  # Parallel methods
  # ====================

  # Multicore lapply: ONLY WORKS ON LINUX OR MAC
  # The “mc” stands for “multicore,” and as you might gather,
  # this function distributes the lapply tasks across multiple
  # CPU cores to be executed in parallel. The downside is that
  # this shared memory approach to parallelism in R is limited
  # by how many cores your computer has.

  if (getOS() != "Windows")
    if (method == "mclapply")
      return(run({mclapply(args, func)}))

  # The beauty of mclapply is that the worker processes are all
  # created as clones of the master right at the point that mclapply
  # is called, so you don't have to worry about reproducing your
  # environment on each of the cluster workers. Unfortunately, that
  # isn't possible on Windows.

  # When using parLapply, you generally have to perform the following
  # additional steps:
  # > Create a PSOCK cluster: cl = makeCluster(noc)
  # > Register the cluster if desired: setDefaultCluster(cl)
  # > Load necessary packages on the cluster workers: clusterExport(cl, c(funcs))
  # > Export necessary data and functions to the global environment of
  #   the cluster workers: clusterEvalQ(cl, c(libs))

  # Parallel lapply
  if (method == "lapply")
    return(runPar({parLapply(cl, args, func)}, noc))

  # Parallel sapply
  if (method == "sapply")
    return(runPar({parSapply(cl, args, func)}, noc))

  # Foreach
  # We can do:
  # > cl = makeCluster(noc)
  # > registerDoParallel(cl)
  # > stopCluster
  # or:
  # > registerDoParallel(noc)
  # > stopImplicitCluster()
  #
  # The idea behind the foreach package is to create
  # ‘a hybrid of the standard for loop and lapply function’
  if (method == "foreach")
  {
    ffunc <- function(func, args)
    {
      printf("Executing foreach...")
      data = foreach(fat = args, .packages = c(), .combine = rbind) %dopar%
          func(fat)
      printfInline("Done\n")
      return(data)
    }
    return(runPar(ffunc(func, args), noc))
  }
  error("Invalid method")
}

benchmark <- function(noc = -1)
{
  # lapply
  runSequencial(method = "lapply")
  runParallel(method = "lapply", noc = noc)

  # foreach
  runSequencial(method = "foreach")
  runParallel(method = "foreach", noc = noc)
  return(0)
}
