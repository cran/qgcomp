## ----metals data, echo=TRUE, results='markup', message=FALSE------------------
library("qgcomp")
library("knitr")
library("ggplot2")
data("metals", package="qgcomp")
head(metals)

## ----linear model and runtime, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
# we save the names of the mixture variables in the variable "Xnm"
Xnm <- c(
    'arsenic','barium','cadmium','calcium','chromium','copper',
    'iron','lead','magnesium','manganese','mercury','selenium','silver',
    'sodium','zinc'
)
covars = c('nitrate','nitrite','sulfate','ph', 'total_alkalinity','total_hardness')



# Example 1: linear model
# Run the model and save the results "qc.fit"
system.time(qc.fit <- qgcomp.noboot(y~.,dat=metals[,c(Xnm, 'y')], family=gaussian()))
#   user  system elapsed 
#  0.011   0.002   0.018 

# contrasting other methods with computational speed
# WQS regression
#system.time(wqs.fit <- gwqs(y~NULL,mix_name=Xnm, data=metals[,c(Xnm, 'y')], family="gaussian"))
#   user  system elapsed 
# 35.775   0.124  36.114 

# Bayesian kernel machine regression (note that the number of iterations here would 
#  need to be >5,000, at minimum, so this underestimates the run time by a factor
#  of 50+
#system.time(bkmr.fit <- kmbayes(y=metals$y, Z=metals[,Xnm], family="gaussian", iter=100))
#   user  system elapsed 
# 81.644   4.194  86.520 


#first note that qgcomp is very fast

# View results: scaled coefficients/weights and statistical inference about
# mixture effect
qc.fit




## ----logistic qgcomp, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

qc.fit2 <- qgcomp.noboot(disease_state~., expnms=Xnm, 
          data = metals[,c(Xnm, 'disease_state')], family=binomial(), 
          q=4)
qcboot.fit2 <- qgcomp.boot(disease_state~., expnms=Xnm, 
          data = metals[,c(Xnm, 'disease_state')], family=binomial(), 
          q=4, B=10,# B should be 200-500+ in practice
          seed=125, rr=FALSE)
qcboot.fit2b <- qgcomp.boot(disease_state~., expnms=Xnm, 
          data = metals[,c(Xnm, 'disease_state')], family=binomial(), 
          q=4, B=10,# B should be 200-500+ in practice
          seed=125, rr=TRUE)


# Compare a qgcomp.noboot fit:
qc.fit2
# and a qgcomp.boot fit:
qcboot.fit2
# and a qgcomp.boot fit, where the risk/prevalence ratio is estimated, 
#  rather than the odds ratio:
qcboot.fit2b


## ----adjusting for covariates, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

qc.fit3 <- qgcomp.noboot(y ~ mage35 + arsenic + barium + cadmium + calcium + chloride + 
                           chromium + copper + iron + lead + magnesium + manganese + 
                           mercury + selenium + silver + sodium + zinc,
                         expnms=Xnm,
                         metals, family=gaussian(), q=4)
qc.fit3
plot(qc.fit3)
qcboot.fit3 <- qgcomp.boot(y ~ mage35 + arsenic + barium + cadmium + calcium + chloride + 
                           chromium + copper + iron + lead + magnesium + manganese + 
                           mercury + selenium + silver + sodium + zinc,
                         expnms=Xnm,
                         metals, family=gaussian(), q=4, B=10,# B should be 200-500+ in practice
                         seed=125)
qc.fit3
plot(qcboot.fit3)
plot(qcboot.fit3, pointwiseref = 3)


## ----non-linear non-hom intro, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

qcboot.fit4 <- qgcomp(y~. + .^2,
                         expnms=Xnm,
                         metals[,c(Xnm, 'y')], family=gaussian(), q=4, B=10, seed=125)
plot(qcboot.fit4)

## ----overall non-linearity, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

qcboot.fit5 <- qgcomp(y~. + .^2,
                         expnms=Xnm,
                         metals[,c(Xnm, 'y')], family=gaussian(), q=4, degree=2, 
                      B=10, rr=FALSE, seed=125)
plot(qcboot.fit5)

## ----graphical non-linearity, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
library(splines)
# find all correlations > 0.6 (this is an arbitrary choice)
cormat = cor(metals[,Xnm])
idx = which(cormat>0.6 & cormat <1.0, arr.ind = TRUE)
newXnm = unique(rownames(idx)) # iron, lead, and cadmium


qc.fit6lin <- qgcomp.boot(y ~ iron + lead + cadmium + 
                         mage35 + arsenic + magnesium + manganese + mercury + 
                         selenium + silver + sodium + zinc,
                         expnms=newXnm,
                         metals, family=gaussian(), q=8, B=10)

qc.fit6nonlin <- qgcomp.boot(y ~ bs(iron) + bs(cadmium) + bs(lead) +
                         mage35 + arsenic + magnesium + manganese + mercury + 
                         selenium + silver + sodium + zinc,
                         expnms=newXnm,
                         metals, family=gaussian(), q=8, B=10, degree=2)

qc.fit6nonhom <- qgcomp.boot(y ~ bs(iron)*bs(lead) + bs(iron)*bs(cadmium) + bs(lead)*bs(cadmium) +
                         mage35 + arsenic + magnesium + manganese + mercury + 
                         selenium + silver + sodium + zinc,
                         expnms=newXnm,
                         metals, family=gaussian(), q=8, B=10, degree=3)
# it helps to place the plots on a common y-axis, which is easy due 
#  to dependence of the qgcomp plotting functions on ggplot
pl.fit6lin <- plot(qc.fit6lin, suppressprint = TRUE, pointwiseref = 4)
pl.fit6nonlin <- plot(qc.fit6nonlin, suppressprint = TRUE, pointwiseref = 4)
pl.fit6nonhom <- plot(qc.fit6nonhom, suppressprint = TRUE, pointwiseref = 4)

pl.fit6lin + coord_cartesian(ylim=c(-0.75, .75)) + 
  ggtitle("Linear fit: mixture of iron, lead, and cadmium")
pl.fit6nonlin + coord_cartesian(ylim=c(-0.75, .75)) + 
  ggtitle("Non-linear fit: mixture of iron, lead, and cadmium")
pl.fit6nonhom + coord_cartesian(ylim=c(-0.75, .75)) + 
  ggtitle("Non-linear, non-homogeneous fit: mixture of iron, lead, and cadmium")

## ----graphical caution, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
qc.overfit <- qgcomp.boot(y ~ bs(iron) + bs(cadmium) + bs(lead) +
                         mage35 + bs(arsenic) + bs(magnesium) + bs(manganese) + bs(mercury) + 
                         bs(selenium) + bs(silver) + bs(sodium) + bs(zinc),
                         expnms=Xnm,
                         metals, family=gaussian(), q=8, B=10, degree=1)
qc.overfit
plot(qc.overfit, pointwiseref = 5)
plot(qc.overfit, flexfit = FALSE, pointwiseref = 5)

## ----non-linear examples, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----

# using indicator terms for each quantile
qc.fit7a <- qgcomp.boot(y ~ factor(iron) + lead + cadmium + 
                         mage35 + arsenic + magnesium + manganese + mercury + 
                         selenium + silver + sodium + zinc,
                         expnms=newXnm,
                         metals, family=gaussian(), q=8, B=10, deg=2)
# underlying fit
summary(qc.fit7a$fit)$coefficients
plot(qc.fit7a)

# interactions between indicator terms
qc.fit7b <- qgcomp.boot(y ~ factor(iron)*factor(lead) + cadmium + 
                         mage35 + arsenic + magnesium + manganese + mercury + 
                         selenium + silver + sodium + zinc,
                         expnms=newXnm,
                         metals, family=gaussian(), q=8, B=10, deg=3)
# underlying fit
#summary(qc.fit7b$fit)$coefficients
plot(qc.fit7b)

# breaks at specific quantiles (these breaks act on the quantized basis)
qc.fit7c <- qgcomp.boot(y ~ I(iron>3)*I(lead>3) + cadmium + 
                         mage35 + arsenic + magnesium + manganese + mercury + 
                         selenium + silver + sodium + zinc,
                         expnms=newXnm,
                         metals, family=gaussian(), q=8, B=10, deg=2)
# underlying fit
summary(qc.fit7b$fit)$coefficients
plot(qc.fit7c)



## ----splines, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
AIC(qc.fit6lin$fit)
AIC(qc.fit6nonlin$fit)
AIC(qc.fit6nonhom$fit)

BIC(qc.fit6lin$fit)
BIC(qc.fit6nonlin$fit)
BIC(qc.fit6nonhom$fit)

## ----time-to-event, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
# non-bootstrapped version estimates a marginal structural model for the 
# confounder-conditional effect
survival::coxph(survival::Surv(disease_time, disease_state) ~ iron + lead + cadmium + 
                         arsenic + magnesium + manganese + mercury + 
                         selenium + silver + sodium + zinc +
                         mage35,
                         data=metals)
qc.survfit1 <- qgcomp.cox.noboot(survival::Surv(disease_time, disease_state) ~ .,expnms=Xnm,
                         data=metals[,c(Xnm, 'disease_time', 'disease_state')], q=4)
qc.survfit1
plot(qc.survfit1)

# bootstrapped version estimates a marginal structural model for the population average effect
#library(survival)
qc.survfit2 <- qgcomp.cox.boot(Surv(disease_time, disease_state) ~ .,expnms=Xnm,
                         data=metals[,c(Xnm, 'disease_time', 'disease_state')], q=4, 
                         B=5, MCsize=1000, parallel=TRUE)
qc.survfit2
p2 = plot(qc.survfit2, suppressprint = TRUE)  
p2 + labs(title="Linear log(hazard ratio), overall and exposure specific")

qc.survfit3 <- qgcomp.cox.boot(Surv(disease_time, disease_state) ~ . + .^2,expnms=Xnm,
                         data=metals[,c(Xnm, 'disease_time', 'disease_state')], q=4, 
                         B=5, MCsize=1000, parallel=TRUE)
qc.survfit3
p3 = plot(qc.survfit3, suppressprint = TRUE) 
p3 + labs(title="Non-linear log(hazard ratio) overall, linear exposure specific ln-HR")

qc.survfit4 <- qgcomp.cox.boot(Surv(disease_time, disease_state) ~ . + .^2,expnms=Xnm,
                         data=metals[,c(Xnm, 'disease_time', 'disease_state')], q=4, 
                         B=5, MCsize=1000, parallel=TRUE, degree=2)
qc.survfit4
# examining the overall hazard ratio as a function of overall exposure
hrs_q = exp(matrix(c(0,0,1,1,2,4,3,9), ncol=2, byrow=TRUE)%*%qc.survfit4$msmfit$coefficients)
colnames(hrs_q) = "Hazard ratio"
print("Hazard ratios by quartiles (min-25%,25-50%, 50-75%, 75%-max)")
hrs_q

p4 = plot(qc.survfit4, suppressprint = TRUE) 
p4 + labs(title="Non-linear log(hazard ratio), overall and exposure specific") 


