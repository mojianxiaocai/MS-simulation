source('MS-NAR_helper.R')
Nsize=100 #node number 
Time=100 #sample number 
beta=c(0.3,0.2,0.3,0.7) #parameter
p11=0.7   #Transition probability
p22=0.4   
Beta = beta[2:4]
sig=1
theta=c(beta,p11,p22,sig,sig)
W = getDyadW(N = Nsize, N1 = 10, delta = -0.8, normalize = TRUE)
G1 = beta[2] * W + beta[3] * diag(1, nrow(W))
G2 = beta[2] * W + beta[4] * diag(1, nrow(W))
G=list(G1,G2)
mu = getMu(G[[1]], rep(beta[1],Nsize))
Gamma0 = getGamma0Approx(G[[1]], sigma = sig)
Y0 = mvrnorm(n = 1, mu = mu, Sigma = Gamma0)


Ymat = matrix(0, nrow = Nsize, ncol = Time + 1)  ### use Ymat to store the simulated data
Ymat[, 1] = as.vector(Y0)  ### the first column of Ymat is assigned Y0
statesNames <- c("0", "1")
mcB <- new("markovchain", states = statesNames,
           transitionMatrix = matrix(c(p11,1-p11,1-p22,p22),
                                     nrow = 2, byrow = TRUE, dimnames = list(statesNames, statesNames)))
st <- rmarkovchain(n = Time, object = mcB, what = "list")
st=as.numeric(st)
for (i in 1:Time) {
  Ymat[, i + 1] = as.vector((rep(beta[1],Nsize))+ Beta[1] * W %*% Ymat[, i] + (Beta[2]+(Beta[3]-Beta[2])*st[i]) * 
                              Ymat[, i] + rnorm(Nsize, sd = sig))  ### follow the NAR model to simulate Y time series
}

