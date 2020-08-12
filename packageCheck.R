#This script is for internal use and was a sanity check. 

library(gtools)
library(rstan)
library(CNVRG)

set.seed(1245)
options(scipen=99)

#Simulate some data
notus <- 50
nsamples <- 5000
nreps <- 100
intensity <- 100 #high intensity means deviates vary less and will be easier to predict

comprop <- matrix(0, ncol = notus, nrow = 2)
indprop <- matrix(0, ncol = notus, nrow = nreps)

#Assemblage 1
comprop[1, ] <- rdirichlet(1, c(rep(15, 5), rep(1, notus - 5)))
#Assemblage 2
comprop[2, ] <- rdirichlet(1, c(rep(1, notus - 5), rep(15, 5)))

#Construct data matrix
com <- matrix(0, ncol = notus, nrow = nreps)
for (i in 1:(nreps / 2)) {
indprop[i, ] <- rdirichlet(1, comprop[1, ] * intensity)
com[i, ] <- rmultinom(1, nsamples, prob = indprop[i, ])
}
for (i in (1 + nreps / 2):nreps) {
indprop[i, ] <- rdirichlet(1, comprop[2, ] * intensity)
com[i, ] <- rmultinom(1, nsamples, prob = indprop[i, ])
}
com <- com + 1   
nsamples <- nsamples + 50

#make a dummy sample vector
samps <- rep("sampleX", nreps)

modelOut <- varHMC(countData = data.frame(samps,com),
                   starts = c(1,51),
                   ends = c(50,100),
                   chains = 2, 
                   burn = 100, 
                   samples = 250, 
                   thinning_rate = 2,
                   cores = 4,
                   params_to_save = c("pi", "p"))

mns <- extract_point_estimate(modelOut = modelOut, countData = data.frame(samps,com), treatments = 2)

par(mfrow = c(1,2))
plot(comprop[1,], mns$pointEstimates_pi[1,2:51])
abline(0,1)
plot(comprop[2,], mns$pointEstimates_pi[2,2:51])
abline(0,1)

plot(indprop[1,], mns$pointEstimates_p[1,2:51])
abline(0,1)
plot(indprop[2,], mns$pointEstimates_p[2,2:51])
abline(0,1)

modelOutVI <- varInf(countData = data.frame(samps,com),
                     starts = c(1,51),
                     ends = c(50,100),
                                 algorithm = "meanfield",
                                 output_samples = 500,
                                 params_to_save = c("p"))

point_est <- extract_point_estimate(modelOut = modelOutVI, countData = data.frame(samps,com), treatments = 2)

par(mfrow = c(1,2))
plot(comprop[1,], point_est$pointEstimates_pi[1,2:dim(data.frame(samps, com))[2]])
abline(0,1)
plot(comprop[2,], point_est$pointEstimates_pi[2,2:dim(data.frame(samps, com))[2]])
abline(0,1)

plot(indprop[1,], point_est$pointEstimates_p[1,2:51])
abline(0,1)
plot(indprop[2,], point_est$pointEstimates_p[2,2:51])
abline(0,1)

ts <- isd_transform(model_output = modelOut, isd_index = 3, countData = data.frame(samps, com))
str(ts)
test <- rstan::extract(modelOut, "pi")
test$pi[1,1,]
ts$pi[1,1,]
0.156/0.1097
#Seems to be working

#now try with VB
ts <- isd_transform(model_output = modelOutVI, isd_index = 3, countData = data.frame(samps, com))
str(ts)
test <- rstan::extract(modelOutVI, "pi")
test$pi[,1,]
ts[[1]]$pi.1.1[1:5]
0.14659900 /0.12746500
#seems to work

diff_abund(model_output = modelOut, countData = data.frame(samps,com))
diff_abund(model_output = modelOutVI, countData = data.frame(samps,com))

shan_HMC <- diversity_calc(model_output = modelOut, countData = data.frame(samps, com), params = c("p","pi"), entropy_measure = "shannon")

test <- rstan::extract(modelOut, "pi")
vegan::diversity(test$pi[1,1,], "shannon")
exp(vegan::diversity(test$pi[1,2,], "shannon")) #ssems to work

simp_HMC <- diversity_calc(model_output = modelOut, countData = data.frame(samps, com), params = c("p","pi"), entropy_measure = "simpson")
1/(vegan::diversity(test$pi[1,2,], "simpson")) #ssems to work

simp_VI <- diversity_calc(model_output = modelOutVI, countData = data.frame(samps, com), params = c("p","pi"), entropy_measure = "simpson", equivalents = F)
test <- rstan::extract(modelOutVI, "pi")
1/(vegan::diversity(test$pi[1,2,], "simpson")) #ssems to work

#now with p values
shan_HMC <- diversity_calc(model_output = modelOut, countData = data.frame(samps, com), params = "p", entropy_measure = "shannon")
test <- rstan::extract(modelOut, "p")
exp(vegan::diversity(test$p[3,7,], "shannon")) #seems good still


#now with p values
shan_Hvi <- diversity_calc(model_output = modelOutVI, countData = data.frame(samps, com), params = "p", entropy_measure = "shannon")
test <- rstan::extract(modelOutVI, "p")
exp(vegan::diversity(test$p[2,2,], "shannon")) #seems good still

#Make sure I get sensible answers from diff abund, using other more variable modeling approaches (intensity of 100 is to high)
#A: we generally do
t.seed(1245)
options(scipen=99)

#Simulate some data
notus <- 50
nsamples <- 5000
nreps <- 100
intensity <- 3 #high intensity means deviates vary less and will be easier to predict

comprop <- matrix(0, ncol = notus, nrow = 2)
indprop <- matrix(0, ncol = notus, nrow = nreps)

#Assemblage 1
comprop[1, ] <- rdirichlet(1, c(rep(15, 5), rep(1, notus - 5)))
#Assemblage 2
comprop[2, ] <- rdirichlet(1, c(rep(1, notus - 5), rep(15, 5)))

#Construct data matrix
com <- matrix(0, ncol = notus, nrow = nreps)
for (i in 1:(nreps / 2)) {
  indprop[i, ] <- rdirichlet(1, comprop[1, ] * intensity)
  com[i, ] <- rmultinom(1, nsamples, prob = indprop[i, ])
}
for (i in (1 + nreps / 2):nreps) {
  indprop[i, ] <- rdirichlet(1, comprop[2, ] * intensity)
  com[i, ] <- rmultinom(1, nsamples, prob = indprop[i, ])
}
com <- com + 1   
nsamples <- nsamples + 50

#make a dummy sample vector
samps <- rep("sampleX", nreps)

modelOutVI <- varInf(countData = data.frame(samps,com),
                     starts = c(1,51),
                     ends = c(50,100),
                     algorithm = "meanfield",
                     output_samples = 500,
                     params_to_save = c("pi", "p"))


diff_abund(model_output = modelOutVI, countData = data.frame(samps,com))

