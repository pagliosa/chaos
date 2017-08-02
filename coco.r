# library(deSolve)
# library(scatterplot3d)
#
#
# # Parameters for the solver
# param <- c(alpha = 10,
#           beta = 8/3,
#           c = 26.48)
#
# # Initial state
# yini <- c(x = 0.01, y = 0.0, z = 0.0)
#
# # Lorenz function
# lorenz <- function(Time, State, Param) {
#   with(as.list(c(State, Param)), {
#     xdot <- alpha * (y - x)
#     ydot <- x * (c - z) - y
#     zdot <- x*y - beta*z
#     return(list(c(xdot, ydot, zdot)))
#   })
# }
#
# # Run function
# runIt <- function(times) {
#   out <- as.data.frame(ode(func = lorenz, y = yini, parms = param, times = times))
#
#   plot3d(x=out[,2],
#                 y=out[,3],
#                 z=out[,4],
#                 type="l",
#                 box=T,
#                 grid=F,
#                 axis=F,
#                 xlab=NULL,
#                 ylab=NULL,
#                 zlab=NULL,
#                 main=NULL)
#
# }
#
# # Run All function combining functions
# runAll <- function() {
#   runIt(seq(0, 100, by=0.01))
# }
