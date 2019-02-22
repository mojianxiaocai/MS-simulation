# define the markovchain object
  statesNames <- c("0", "1")
  mcB <- new("markovchain", states = statesNames,
             transitionMatrix = matrix(c(p11,1-p11,1-p22,p22),
                                       nrow = 2, byrow = TRUE, dimnames = list(statesNames, statesNames)))
  outs <- rmarkovchain(n = 10, object = mcB, what = "list")
as.numeric(outs[2])  
