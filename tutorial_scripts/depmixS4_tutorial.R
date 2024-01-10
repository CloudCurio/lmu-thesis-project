library(depmixS4)
data(speed)
set.seed(1)
mod <- depmix(response = rt ~ 1, data = speed, nstates = 2, trstart = runif(4))

fm <- fit(mod, emc = em.control(rand = FALSE))

#covariates on transition parameters
set.seed(1)
mod <- depmix(rt ~ 1, data = speed, nstates = 2, family = gaussian(), 
              transition = ~ scale(Pacc), instart = runif(2))
fm <- fit(mod, verbose = FALSE, emc = em.control(rand = FALSE))

summary(fm, which = "transition")
#Multivariate data
set.seed(1)
mod <- depmix(list(rt ~ 1, corr ~ 1), data = speed, nstates = 2,
              family = list(gaussian(), multinomial("identity")),
              transition = ~scale(Pacc), instart = runif(2))
fm <- fit(mod, verbose = FALSE, emc = em.control(rand = FALSE))

summary(fm, which = "response")

#fixing and constraining parameters
#first, fit a model without constraints
trst <- c(0.9, 0.1, 0, 0, 0.1, 0.9, 0, 0)
mod <- depmix(list(rt ~ 1, corr ~ 1), data = speed, transition = ~ Pacc,
              nstates = 2, family = list(gaussian(), multinomial("identity")),
              trstart = trst, instart = c(0.99, 0.01))
fm1 <- fit(mod, verbose = FALSE, emc = em.control(rand = FALSE))

#now, use the fitted values from previous model to contrain the regression
#coefficients on the transition matrix (parameters №6 and 10)
pars <- c(unlist(getpars(fm1)))
pars[6] <- pars[10] <- 11
pars[1] <- 0
pars[2] <- 1
pars[13] <- pars[14] <- 0.5
fm1 <- setpars(mod, pars)
conpat <- c(0, 0, rep(c(0, 1), 4), 1, 1, 0, 0, 1, 1, 1, 1)
conpat[6] <- conpat[10] <- 2
fm2 <- fit(fm1, equal = conpat)

summary(fm2)

#adding covariates on the prior probabilities
#mix model instead of a depmix model is used since data forms 
#independent observations
data("balance")
set.seed(1)
mod <- mix(list(d1 ~ 1, d2 ~ 1, d3 ~ 1, d4 ~ 1), data = balance,
           nstates = 3, family = list(multinomial("identity"),
                                      multinomial("identity"), 
                                      multinomial("identity"),
                                      multinomial("identity")), 
           respstart = runif(24), prior = ~ age,
           initdata = balance)
fm <- fit(mod, verbose = FALSE, emc = em.control(rand = FALSE))
fm
summary(fm, which = "prior")

#implementing viterbi algorithm for predictions
# 2-state model on rt and corr from speed data set
# with Pacc as covariate on the transition matrix
# ntimes is used to specify the lengths of 3 separate series
mod <- depmix(list(rt~1,corr~1),data=speed,transition=~Pacc,nstates=2,
              family=list(gaussian(),multinomial("identity")),ntimes=c(168,134,137))
fmod <- fit(mod)
# result of viterbi is stored in a depmix-fitted object in slot "posterior"
identical(viterbi(fmod),fmod@posterior)
vit1 <- viterbi(fmod)

new_data_fmod <- fmod
new_data_fmod@posterior