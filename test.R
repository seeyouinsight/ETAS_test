#{r setup, include=FALSE}
# load the functions needed for this analysis
#install.packages("inlabru")
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#install.packages("bayesianETAS")
#install.packages("ggstar")
#region install.packages("R.filesets")


knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(inlabru)
library(ggplot2)
library(dplyr)
library(foreach)
library(parallel)
library(INLA)

set.seed(1)
# extract sample from standard gaussian distribution
st.gaussian.sample <- rnorm(10000)
# transform in Gamma(1,2)
gamma.values <- gamma.t(st.gaussian.sample, 1, 2)
# transform in Uniform(0,1)
unif.values <- unif.t(st.gaussian.sample, 0, 5)
# transform in LogGaussian(0.5, 1)
loggaus.values <- loggaus.t(st.gaussian.sample, 0.5, 0.5)

ggplot() + 
  geom_density(aes(x = gamma.values, color = 'Gamma(1,2)')) +
  geom_density(aes(x = unif.values, color = 'Unif(0,5)')) +
  geom_density(aes(x = loggaus.values, color = 'Log-Gaussian(0.5, 0.5)')) +
  xlim(0,7) + 
  xlab(~eta(theta))


link.f.be <- list(mu = \(x) gamma.t(x, 0.1, 0.1), 
                  K = \(x) loggaus.t(x, -1, 2), 
                  alpha = \(x) unif.t(x, 0, 10), 
                  c_ = \(x) unif.t(x, 0, 10), 
                  p = \(x) unif.t(x, 1, 10))

M0 = 2.99
T1 = 0; T2 = 357

# import data
dd.ama <- read.csv2(file = 'data_M3.0.csv', header = TRUE, sep = ',') %>%
  mutate(time_date = as.POSIXct(paste0(year,'-',month,'-',day,' ',hr,':',min,':',sec)),
         time.diff = as.numeric(difftime(time_date, min(time_date) - 1, units = 'days')),
         Mw = as.numeric(Mw),
         Lon = as.numeric(Lon),
         Lat = as.numeric(Lat),
         Mw.class = cut(Mw, breaks = c(M0, 5, 7))) %>%
  arrange(time_date)

# data for Inlabru
data.bru <- data.frame(ts = dd.ama$time.diff, 
                       magnitudes = dd.ama$Mw) %>%
  mutate(idx.p = 1:nrow(dd.ama))


# Fitting the model
# set initial values for the parameters in the internal scale to get reasonable initial values (default is 0)
th.init <- list(th.mu = 0.5,
                th.K = 0.5,
                th.alpha = -2,
                th.c = -2,
                th.p = -2) 

# options for inlabru 
bru.opt.list <- list(bru_verbose = 4, # type of visual output 
                     bru_max_iter = 100, # maximum number of iterations
                     #bru_method = list(max_step = 0.5, rel_tol = 0.01),
                     inla.mode = 'experimental', # type of inla algorithm
                     bru_initial = th.init) # parameters initial values


fit_etas <- Hawkes.bru(sample.s = data.bru, # data 
                       M0 = M0, # magnitude of completeness
                       T1 = 0, T2 = T2, # time domain
                       link.functions = link.f.be, # link functions
                       coef.t. = 1, # binning parameter (delta)
                       delta.t. = 0.1, # binning parameter (Delta)
                       N.max. = 3, # binning parameter (n.max)
                       bru.opt = bru.opt.list)

ggplot(fit_etas$bru_iinla$track, aes(x = iteration, y = mode)) + 
  geom_line() + 
  facet_wrap(facets = vars(effect), scales = 'free')


inlabru:::make_track_plots(fit_etas)$default


# substitute 'th.mu' with 'th.K', 'th.alpha', 'th.c' or 'th.p' to explore the posterior of the others parameters
plot(fit_etas, 'th.mu')

# posterior of parameter mu in ETAS scale
post.mu <- data.frame(inla.tmarginal(link.f.be$mu, 
                                     fit_etas$marginals.fixed$th.mu),
                      param = 'mu')
# posterior of parameter mu in ETAS scale
post.K <- data.frame(inla.tmarginal(link.f.be$K, 
                                    fit_etas$marginals.fixed$th.K),
                     param = 'K')
# posterior of parameter mu in ETAS scale
post.alpha <- data.frame(inla.tmarginal(link.f.be$alpha, 
                                        fit_etas$marginals.fixed$th.alpha),
                         param = 'alpha')
# posterior of parameter mu in ETAS scale
post.c <- data.frame(inla.tmarginal(link.f.be$c_, 
                                    fit_etas$marginals.fixed$th.c),
                     param = 'c')
# posterior of parameter mu in ETAS scale
post.p <- data.frame(inla.tmarginal(link.f.be$p, 
                                    fit_etas$marginals.fixed$th.p),
                     param = 'p')

ggplot(rbind(post.mu, post.K, post.alpha, post.c, post.p), aes(x,y)) + 
  geom_line() + 
  facet_wrap(facets = vars(param), scales = 'free', labeller = label_parsed) + 
  xlab('param') + 
  ylab('pdf')


inla.rmarginal(10, inla.tmarginal(link.f.be$mu, 
                                  fit_etas$marginals.fixed$th.mu))

# predict function
lambda.N.post <- predict(fit_etas, # model fit 
                         data.frame(), # data (empty because the history of the process is passed to the function below directly)
                         ~ lambda.N(th.mu, th.K, th.alpha, th.c, th.p, 
                                    T1, T2, M0, 
                                    data.bru,
                                    link.f.be)) # target function

c(lambda.N.post[1:5], true = nrow(data.bru))

N.post <- predict(fit_etas, data.frame(), 
                  ~ data.frame(N = 800:1200, 
                               pdf = dpois(800:1200, 
                                           lambda.N(th.mu, th.K, 
                                                    th.alpha, th.c, th.p, 
                                                    T1, T2, M0, 
                                                    data.bru,
                                                    link.f.be)) ))

head(N.post)


ggplot(N.post, aes(x = N, y = mean)) + 
  geom_line(color = 'darkblue') + 
  geom_ribbon(aes(xmax = N, xmin = N, ymin = q0.025, ymax = q0.975), alpha = 0.2,
              fill = 'blue') + 
  geom_vline(xintercept = nrow(data.bru), linetype = 3) + 
  geom_line(data = data.frame(x = 800:1200,
                              y = dpois(800:1200, 967.6161)),
            aes(x,y), color = 'red', linetype = 2) + 
  ylab('pdf')

# Sample from the posterior distribution of functions of multiple parameters

N.samp <- generate(fit_etas, data.frame(), 
                   ~ data.frame(N = 800:1200, 
                                pdf = dpois(800:1200, 
                                            lambda.N(th.mu, th.K, 
                                                     th.alpha, th.c, th.p, 
                                                     T1, T2, M0, 
                                                     data.bru,
                                                     link.f.be)) ), 
                   n.samples = 1, seed = 1)

ggplot(N.post, aes(x = N, y = mean)) + 
  geom_line(color = 'darkblue') + 
  geom_line(data = data.frame(x = 800:1200,
                              y = dpois(800:1200, 967.6161)),
            aes(x,y), color = 'red', linetype = 2) +
  geom_ribbon(aes(xmax = N, xmin = N, ymin = q0.025, ymax = q0.975), alpha = 0.2,
              fill = 'blue') + 
  geom_line(data = N.samp[[1]], 
            aes(x = N, y = pdf), linetype = 3) + 
  ylab('pdf') 










