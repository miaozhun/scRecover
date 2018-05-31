# MLE of parameters of ZINB
mleZINB <- function(counts_1){
  counts_1 <- as.numeric(counts_1)
  if(sum(counts_1 == 0) > 0){
    if(sum(counts_1 == 0) == length(counts_1)){
      theta_1 <- 1
      mu_1 <- 0
      size_1 <- 1
      prob_1 <- size_1/(size_1 + mu_1)
    }else{
      options(show.error.messages = FALSE)
      zinb_try <- try(gamlssML(counts_1, family="ZINBI"), silent=TRUE)
      options(show.error.messages = TRUE)
      if('try-error' %in% class(zinb_try)){
        zinb_try_twice <- try(zeroinfl(formula = counts_1 ~ 1 | 1, dist = "negbin"), silent=TRUE)
        if('try-error' %in% class(zinb_try_twice)){
          # print("MLE of ZINB failed!");
          parameters <- c(NA, NA, NA, NA)
          names(parameters) <- c("theta", "mu", "size", "prob")
          return(parameters)
        }else{
          zinb_1 <- zinb_try_twice
          theta_1 <- plogis(zinb_1$coefficients$zero);names(theta_1) <- NULL
          mu_1 <- exp(zinb_1$coefficients$count);names(mu_1) <- NULL
          size_1 <- zinb_1$theta;names(size_1) <- NULL
          prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
        }
      }else{
        zinb_1 <- zinb_try
        theta_1 <- zinb_1$nu;names(theta_1) <- NULL
        mu_1 <- zinb_1$mu;names(mu_1) <- NULL
        size_1 <- 1/zinb_1$sigma;names(size_1) <- NULL
        prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
      }
    }
  }else{
    op <- options(warn=2)
    nb_try <- try(glm.nb(formula = counts_1 ~ 1), silent=TRUE)
    options(op)
    if('try-error' %in% class(nb_try)){
      nb_try_twice <- try(fitdistr(counts_1, "Negative Binomial"), silent=TRUE)
      if('try-error' %in% class(nb_try_twice)){
        nb_try_again <- try(mle2(counts_1~dnbinom(mu=exp(logmu),size=1/invk), data=data.frame(counts_1), start=list(logmu=0,invk=1), method="L-BFGS-B", lower=c(logmu=-Inf,invk=1e-8)), silent=TRUE)
        if('try-error' %in% class(nb_try_again)){
          nb_try_fourth <- try(glm.nb(formula = counts_1 ~ 1), silent=TRUE)
          if('try-error' %in% class(nb_try_fourth)){
            # print("MLE of NB failed!");
            parameters <- c(NA, NA, NA, NA)
            names(parameters) <- c("theta", "mu", "size", "prob")
            return(parameters)
          }else{
            nb_1 <- nb_try_fourth
            theta_1 <- 0
            mu_1 <- exp(nb_1$coefficients);names(mu_1) <- NULL
            size_1 <- nb_1$theta;names(size_1) <- NULL
            prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
          }
        }else{
          nb_1 <- nb_try_again
          theta_1 <- 0
          mu_1 <- exp(nb_1@coef["logmu"]);names(mu_1) <- NULL
          size_1 <- 1/nb_1@coef["invk"];names(size_1) <- NULL
          prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
        }
      }else{
        nb_1 <- nb_try_twice
        theta_1 <- 0
        mu_1 <- nb_1$estimate["mu"];names(mu_1) <- NULL
        size_1 <- nb_1$estimate["size"];names(size_1) <- NULL
        prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
      }
    }else{
      nb_1 <- nb_try
      theta_1 <- 0
      mu_1 <- exp(nb_1$coefficients);names(mu_1) <- NULL
      size_1 <- nb_1$theta;names(size_1) <- NULL
      prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
    }
  }
  parameters <- c(theta_1, mu_1, size_1, prob_1)
  names(parameters) <- c("theta", "mu", "size", "prob")
  return(parameters)
}




