## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----how_to_install, eval=FALSE------------------------------------------
#  #stable CRAN version
#  install.packages("tcensReg")
#  
#  # or ----
#  
#  #active devel. GitHub version
#  install.packages("devtools")
#  devtools::install_github("williazo/tcensReg")

## ----single_pop_data_gen-------------------------------------------------
#loading the package
library(tcensReg)  
mu <- 0.5
sigma <- 0.5
a <- 0
#generate random values from the truncated normal distribution using tcensReg function
y_star <- rtnorm(n=1000, mu=mu, sd=sigma, a=a)
#note that the lowerbound will always be non-negative
round(range(y_star), 3)

## ------------------------------------------------------------------------
nu <- 0.25
y <- ifelse(y_star <= nu, nu, y_star)
#calculating the number of censored observations
sum(y == nu)/length(y) 
#collecting the uncensored and censored data together
dt <- data.frame(y_star, y) 

## ----echo=FALSE, warning=FALSE, message=FALSE, fig.align='center',fig.height=4, fig.width=6,tidy=FALSE----
library(ggplot2)
library(viridis)
ggplot(data=dt, aes(x=y_star, y=..density..))+
  geom_histogram(binwidth=0.25, col=viridis(1),fill=viridis(1), alpha=0.6)+
  stat_density(bw=0.5, col=viridis(1), fill=viridis(1), alpha=0.3)+
  scale_x_continuous(breaks=seq(0, 2, 0.25))+
  ylab("Density")+
  xlab(expression(Y^"*"))+
    theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), panel.grid.major=element_blank(), 
        panel.background=element_blank(), panel.border=element_rect(fill=NA, color="black"),
        axis.text.x=element_text(size=10), legend.title=element_text(size=15),
        legend.text=element_text(size=15),  legend.key.width=unit(15, units="mm"), strip.text=element_text(size=20),
        axis.title=element_text(size=20))

## ----echo=FALSE, warning=FALSE, message=FALSE, fig.align='center',fig.height=4, fig.width=6,tidy=FALSE----
ggplot(data=dt, aes(x=y, y=..density..))+
  geom_histogram(binwidth=0.25, col=viridis(1, begin=0.5),fill=viridis(1, begin=0.5), alpha=0.6)+
  stat_density(bw=0.25, col=viridis(1, begin=0.5), fill=viridis(1, begin=0.5), alpha=0.3)+
  scale_x_continuous(breaks=seq(0, 2, 0.25))+
  ylab("Density")+
  xlab(expression(Y))+
    theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), panel.grid.major=element_blank(), 
        panel.background=element_blank(), panel.border=element_rect(fill=NA, color="black"),
        axis.text.x=element_text(size=10), legend.title=element_text(size=15),
        legend.text=element_text(size=15),  legend.key.width=unit(15, units="mm"), strip.text=element_text(size=20),
        axis.title=element_text(size=20))

## ----message=FALSE-------------------------------------------------------
tcensReg(y ~ 1, data=dt, a=0, v=0.25)

## ----message=FALSE-------------------------------------------------------
#tcensReg model
output <- tcensReg(y ~ 1, data=dt, a=a, v=nu)
#extracting the point estimates
tcensReg_est <- output$theta 
#exponentiating the estimate of log_sigma to obtain sigma
tcensReg_est[2] <- exp(tcensReg_est[2]) 

#OLS model
lm_output <- lm(y ~ 1, data=dt) 
lm_est <- c(coef(lm_output), summary(lm_output)$sigma)
#censored only model, i.e., Tobit model
cens_output <- tcensReg(y ~ 1, data=dt, v=nu) 
cens_est <- cens_output$theta
cens_est[2] <- exp(cens_est[2])


results_df <- data.frame(rbind(c(mu, sigma),
                               t(tcensReg_est),
                               lm_est,
                               t(cens_est)))
names(results_df) <- c("mu", "sigma")
row.names(results_df) <- c("Truth", "tcensReg", "Normal MLE", "Tobit")
results_df$mu_bias <- abs(results_df$mu - mu)
results_df$sigma_bias <- abs(results_df$sigma - sigma)

knitr::kable(results_df, format="markdown", digits=4)

## ----gen_twopopdata------------------------------------------------------
mu_1 <- 0.5
mu_2 <- 1
sigma_1 <- 0.25
sigma_2 <- 2
a <- 0

y_1_star <- rtnorm(1000, mu = mu_1, sd = sigma_1, a = a)
y_2_star <- rtnorm(1000, mu = mu_2, sd = sigma_2, a = a)
df <- data.frame(y_star = c(y_1_star, y_2_star), 
                 group = c(rep("Population 1", length(y_1_star)),
                           rep("Population 2", length(y_2_star))))

## ----two_pop_graph, echo=FALSE, warning=FALSE, message = FALSE, fig.align='center',fig.height=4, fig.width=6,tidy=FALSE----
ggplot(data = df, aes(x = y_star, y = ..density.., group = group, fill = group, col = group))+
  stat_density(bw = 0.5, alpha = 0.3)+
  ylab("Density")+
  xlab(expression(Y^"*"))+
  geom_vline(xintercept = 0.5, lty = 2, col = viridis(1, begin = 0.25))+
  geom_vline(xintercept = 1, lty = 2, col = viridis(1, begin = 0))+
  scale_color_manual(name = "", values = viridis(2, begin = 0.25, end = 0))+
  scale_fill_manual(name = "", values = viridis(2, begin = 0.25, end = 0))+
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), panel.grid.major = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(size = 10), legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),  legend.key.width = unit(15, units = "mm"), strip.text = element_text(size = 20),
        axis.title = element_text(size = 20))

## ----two_pop_censoring---------------------------------------------------
nu <- 0.25
df$y <- ifelse(df$y_star<=nu, nu, df$y_star)

## ----fit_sepvarmod-------------------------------------------------------
mod_result <- tcensReg_sepvar(y ~ group, a=a, v=nu, group_var="group", method="maxLik", data=df)
mod_result
sepvar_est <- mod_result$theta
sepvar_est[3:4] <- exp(sepvar_est[3:4])

results_df <- data.frame(rbind(c(mu_1, mu_2, sigma_1, sigma_2),
                               t(sepvar_est)))
names(results_df) <- c("mu_1", "mu_2", "sigma_1", "sigma_2")
row.names(results_df) <- c("Truth", "tcensReg")
results_df$mu1_bias <- abs(results_df$mu_1 - mu_1)
results_df$mu2_bias <- abs(results_df$mu_2 - mu_2)
results_df$sigma1_bias <- abs(results_df$sigma_1 - sigma_1)
results_df$sigma2_bias <- abs(results_df$sigma_2 - sigma_2)

knitr::kable(results_df, format="markdown", digits=4)

## ----speed comparison, message=FALSE, warning=FALSE----------------------
library(microbenchmark)
#testing the censored-only regression
library(censReg)
cens <- microbenchmark(tcensReg_method = tcensReg(y ~ 1, data=dt, v=nu, method="Newton"),
               censReg_method = censReg(y ~ 1, left=nu, data=dt))
knitr::kable(summary(cens), format="markdown", digits=4)

#point estimates are equivalent
tcensReg_est <- as.numeric(tcensReg(y ~ 1, data=dt, v=nu, method="Newton")$theta)
censReg_est <- as.numeric(coef(censReg(y ~ 1, left=nu, data=dt)))
all.equal(tcensReg_est, censReg_est)

#testing the truncated-only regression
library(truncreg)
trunc <- microbenchmark(
  tcensReg_method = tcensReg(y_star ~ 1, data=dt, a=a, method="Newton"),
  truncreg_method = truncreg(y_star ~ 1, point=a, data=dt))
knitr::kable(summary(trunc), format="markdown", digits=4)
tcensReg_est <- as.numeric(tcensReg(y_star ~ 1, data=dt, a=a, method="Newton")$theta)
#note truncreg returns sigma not log_sigma so we need to exponentiate our value
tcensReg_est[2] <- exp(tcensReg_est[2])
truncreg_est <- as.numeric(coef(truncreg(y_star ~ 1, point=a, data=dt)))
all.equal(tcensReg_est, truncreg_est)

