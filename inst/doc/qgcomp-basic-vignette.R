## ----invisibles, echo=FALSE, results='markup', message=FALSE------------------
library("knitr")
#library("gWQS")

## ----first step, echo=TRUE, results='markup', message=FALSE-------------------
library("qgcomp")
set.seed(543210)
qdat = simdata_quantized(n=5000, outcomtype="continuous", cor=c(.95, -0.3), b0=0, coef=c(0.25, -0.1, 0.05), q=4)
head(qdat)
cor(qdat[,c("x1", "x2", "x3")])
qgcomp(y~x1+x2+x3, expnms=c("x1", "x2", "x3"), data=qdat)

## ----metals data, echo=TRUE, results='markup', message=FALSE------------------
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
system.time(qc.fit <- qgcomp.glm.noboot(y~.,dat=metals[,c(Xnm, 'y')], family=gaussian()))
#   user  system elapsed 
#  0.011   0.002   0.018 

# contrasting other methods with computational speed
# WQS regression (v3.0.1 of gWQS package)
#system.time(wqs.fit <- gWQS::gwqs(y~wqs,mix_name=Xnm, data=metals[,c(Xnm, 'y')], family="gaussian"))
#   user  system elapsed 
# 35.775   0.124  36.114 

# Bayesian kernel machine regression (note that the number of iterations here would 
#  need to be >5,000, at minimum, so this underestimates the run time by a factor
#  of 50+
#system.time(bkmr.fit <- kmbayes(y=metals$y, Z=metals[,Xnm], family="gaussian", iter=100))
#   user  system elapsed 
# 81.644   4.194  86.520 

## ----linear model and runtimeb, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
# View results: scaled coefficients/weights and statistical inference about
# mixture effect
qc.fit

## ----linear model and runtime c, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
# quantized data
head(qc.fit$qx)

## ----linear model and runtime d, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
# regression with quantized data
qc.fit$qx$y = qc.fit$fit$data$y # first bring outcome back into the quantized data
newfit <- lm(y ~ arsenic_q + barium_q + cadmium_q + calcium_q + chromium_q + copper_q + 
    iron_q + lead_q + magnesium_q + manganese_q + mercury_q + selenium_q + 
    silver_q + sodium_q + zinc_q, data=qc.fit$qx)
newfit

## ----linear model and runtime e, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
sum(newfit$coefficients[-1]) # sum of all coefficients excluding intercept and confounders, if any
coef(qc.fit) # overall effect and intercept from qgcomp fit

## ----logistic qgcomp, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

qc.fit2 <- qgcomp.glm.noboot(disease_state~., expnms=Xnm, 
          data = metals[,c(Xnm, 'disease_state')], family=binomial(), 
          q=4)
qcboot.fit2 <- qgcomp.glm.boot(disease_state~., expnms=Xnm, 
          data = metals[,c(Xnm, 'disease_state')], family=binomial(), 
          q=4, B=10,# B should be 200-500+ in practice
          seed=125, rr=FALSE)
qcboot.fit2b <- qgcomp.glm.boot(disease_state~., expnms=Xnm, 
          data = metals[,c(Xnm, 'disease_state')], family=binomial(), 
          q=4, B=10,# B should be 200-500+ in practice
          seed=125, rr=TRUE)

## ----logistic qgcompb, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
qc.fit2

## ----logistic qgcompc, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
qcboot.fit2

## ----logistic qgcompd, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
qcboot.fit2b

## ----adj4cov a, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

qc.fit3 <- qgcomp.glm.noboot(y ~ mage35 + arsenic + barium + cadmium + calcium + chloride + 
                           chromium + copper + iron + lead + magnesium + manganese + 
                           mercury + selenium + silver + sodium + zinc,
                         expnms=Xnm,
                         metals, family=gaussian(), q=4)
qc.fit3
plot(qc.fit3)

## ----adj4cov b, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
qcboot.fit3 <- qgcomp.glm.boot(y ~ mage35 + arsenic + barium + cadmium + calcium + chloride + 
                           chromium + copper + iron + lead + magnesium + manganese + 
                           mercury + selenium + silver + sodium + zinc,
                         expnms=Xnm,
                         metals, family=gaussian(), q=4, B=10,# B should be 200-500+ in practice
                         seed=125)
qcboot.fit3
qcee.fit3 <- qgcomp.glm.ee(y ~ mage35 + arsenic + barium + cadmium + calcium + chloride + 
                           chromium + copper + iron + lead + magnesium + manganese + 
                           mercury + selenium + silver + sodium + zinc,
                         expnms=Xnm,
                         metals, family=gaussian(), q=4)
qcee.fit3

## ----adj4cov c, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
plot(qcee.fit3, pointwiseref = 3, flexfit = FALSE)
plot(qcboot.fit3, pointwiseref = 3, flexfit = FALSE)

## ----adj4cov d, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
pointwisebound.boot(qcboot.fit3, pointwiseref=3)
qgcomp:::modelbound.boot(qcboot.fit3)

## ----n-lin non-hom intro, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

qcboot.fit4 <- qgcomp(y~. + .^2,
                         expnms=Xnm,
                         metals[,c(Xnm, 'y')], family=gaussian(), q=4, B=10, seed=125)
plot(qcboot.fit4)

## ----ovrl n-lin, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

qcboot.fit5 <- qgcomp(y~. + .^2,
                         expnms=Xnm,
                         metals[,c(Xnm, 'y')], family=gaussian(), q=4, degree=2, 
                      B=10, rr=FALSE, seed=125)
plot(qcboot.fit5)
qcee.fit5b <- qgcomp.glm.ee(y~. + .^2,
                         expnms=Xnm,
                         metals[,c(Xnm, 'y')], family=gaussian(), q=4, degree=2, 
                         rr=FALSE)
plot(qcee.fit5b)

## ----ovrl n-linb, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
modelbound.boot(qcboot.fit5)
pointwisebound.boot(qcboot.fit5)
pointwisebound.noboot(qcee.fit5b)

## ----ovrl n-lin psi interp, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
qcboot.fit5

## ----graf n-lin 1, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
library(splines)
# find all correlations > 0.6 (this is an arbitrary choice)
cormat = cor(metals[,Xnm])
idx = which(cormat>0.6 & cormat <1.0, arr.ind = TRUE)
newXnm = unique(rownames(idx)) # iron, lead, and cadmium


qc.fit6lin <- qgcomp.glm.boot(y ~ iron + lead + cadmium + 
                         mage35 + arsenic + magnesium + manganese + mercury + 
                         selenium + silver + sodium + zinc,
                         expnms=newXnm,
                         metals, family=gaussian(), q=8, B=10)

qc.fit6nonlin <- qgcomp.glm.boot(y ~ bs(iron) + bs(cadmium) + bs(lead) +
                         mage35 + arsenic + magnesium + manganese + mercury + 
                         selenium + silver + sodium + zinc,
                         expnms=newXnm,
                         metals, family=gaussian(), q=8, B=10, degree=2)

qc.fit6nonhom <- qgcomp.glm.boot(y ~ bs(iron)*bs(lead) + bs(iron)*bs(cadmium) + bs(lead)*bs(cadmium) +
                         mage35 + arsenic + magnesium + manganese + mercury + 
                         selenium + silver + sodium + zinc,
                         expnms=newXnm,
                         metals, family=gaussian(), q=8, B=10, degree=3)

## ----graf n-lin 1b, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
pl.fit6lin <- plot(qc.fit6lin, suppressprint = TRUE, pointwiseref = 4)
pl.fit6lin + coord_cartesian(ylim=c(-0.75, .75)) + 
  ggtitle("Linear fit: mixture of iron, lead, and cadmium")

## ----graf n-lin 2, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
pl.fit6nonlin <- plot(qc.fit6nonlin, suppressprint = TRUE, pointwiseref = 4)
pl.fit6nonlin + coord_cartesian(ylim=c(-0.75, .75)) + 
  ggtitle("Non-linear fit: mixture of iron, lead, and cadmium")

## ----graf n-lin 3, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
pl.fit6nonhom <- plot(qc.fit6nonhom, suppressprint = TRUE, pointwiseref = 4)
pl.fit6nonhom + coord_cartesian(ylim=c(-0.75, .75)) + 
  ggtitle("Non-linear, non-homogeneous fit: mixture of iron, lead, and cadmium")

## ----grafwarn, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
qc.overfit <- qgcomp.glm.boot(y ~ bs(iron) + bs(cadmium) + bs(lead) +
                         mage35 + bs(arsenic) + bs(magnesium) + bs(manganese) + bs(mercury) + 
                         bs(selenium) + bs(silver) + bs(sodium) + bs(zinc),
                         expnms=Xnm,
                         metals, family=gaussian(), q=8, B=10, degree=1)
qc.overfit
plot(qc.overfit, pointwiseref = 5)

## ----grafwarn 2, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
plot(qc.overfit, flexfit = FALSE, pointwiseref = 5)

## ----n-lin exs, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
qc.fit7a <- qgcomp.glm.boot(y ~ factor(iron) + lead + cadmium + 
                         mage35 + arsenic + magnesium + manganese + mercury + 
                         selenium + silver + sodium + zinc,
                         expnms=newXnm,
                         metals, family=gaussian(), q=8, B=20, deg=2)
# underlying fit
summary(qc.fit7a$fit)$coefficients
plot(qc.fit7a)

## ----n-lin exs 2, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
qc.fit7b <- qgcomp.glm.boot(y ~ factor(iron)*factor(lead) + cadmium + 
                         mage35 + arsenic + magnesium + manganese + mercury + 
                         selenium + silver + sodium + zinc,
                         expnms=newXnm,
                         metals, family=gaussian(), q=8, B=10, deg=3)
# underlying fit
#summary(qc.fit7b$fit)$coefficients
plot(qc.fit7b)

## ----n-lin exs 3, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
qc.fit7c <- qgcomp.glm.boot(y ~ I(iron>4)*I(lead>4) + cadmium + 
                         mage35 + arsenic + magnesium + manganese + mercury + 
                         selenium + silver + sodium + zinc,
                         expnms=newXnm,
                         metals, family=gaussian(), q=8, B=10, deg=2)
# underlying fit
summary(qc.fit7c$fit)$coefficients
plot(qc.fit7c)

## ----splines, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
AIC(qc.fit6lin$fit)
AIC(qc.fit6nonlin$fit)
AIC(qc.fit6nonhom$fit)

BIC(qc.fit6lin$fit)
BIC(qc.fit6nonlin$fit)
BIC(qc.fit6nonhom$fit)

## ----parend, echo=TRUE--------------------------------------------------------
# return to standard processing
future::plan(future::sequential) # return to standard evaluation

