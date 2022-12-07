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
system.time(qc.fit <- qgcomp.noboot(y~.,dat=metals[,c(Xnm, 'y')], family=gaussian()))
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

## ----logistic qgcompb, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
qc.fit2

## ----logistic qgcompc, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
qcboot.fit2

## ----logistic qgcompd, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
qcboot.fit2b

## ----adjusting for covariates a, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

qc.fit3 <- qgcomp.noboot(y ~ mage35 + arsenic + barium + cadmium + calcium + chloride + 
                           chromium + copper + iron + lead + magnesium + manganese + 
                           mercury + selenium + silver + sodium + zinc,
                         expnms=Xnm,
                         metals, family=gaussian(), q=4)
qc.fit3
plot(qc.fit3)

## ----adjusting for covariates b, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
qcboot.fit3 <- qgcomp.boot(y ~ mage35 + arsenic + barium + cadmium + calcium + chloride + 
                           chromium + copper + iron + lead + magnesium + manganese + 
                           mercury + selenium + silver + sodium + zinc,
                         expnms=Xnm,
                         metals, family=gaussian(), q=4, B=10,# B should be 200-500+ in practice
                         seed=125)
qcboot.fit3
p = plot(qcboot.fit3)

## ----adjusting for covariates c, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
plot(qcboot.fit3, pointwiseref = 3)

## ----adjusting for covariates d, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
pointwisebound.boot(qcboot.fit3, pointwiseref=3)
qgcomp:::modelbound.boot(qcboot.fit3)

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

## ----overall non-linearityb, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
qgcomp::pointwisebound.boot(qcboot.fit5)
qgcomp:::modelbound.boot(qcboot.fit5)

## ----overall non-linearity psi interp, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
qcboot.fit5

## ----graphical non-linearity 1, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
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

## ----graphical non-linearity 1b, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
pl.fit6lin <- plot(qc.fit6lin, suppressprint = TRUE, pointwiseref = 4)
pl.fit6lin + coord_cartesian(ylim=c(-0.75, .75)) + 
  ggtitle("Linear fit: mixture of iron, lead, and cadmium")

## ----graphical non-linearity 2, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
pl.fit6nonlin <- plot(qc.fit6nonlin, suppressprint = TRUE, pointwiseref = 4)
pl.fit6nonlin + coord_cartesian(ylim=c(-0.75, .75)) + 
  ggtitle("Non-linear fit: mixture of iron, lead, and cadmium")

## ----graphical non-linearity 3, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
pl.fit6nonhom <- plot(qc.fit6nonhom, suppressprint = TRUE, pointwiseref = 4)
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

## ----graphical caution 2, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
plot(qc.overfit, flexfit = FALSE, pointwiseref = 5)

## ----non-linear examples, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
qc.fit7a <- qgcomp.boot(y ~ factor(iron) + lead + cadmium + 
                         mage35 + arsenic + magnesium + manganese + mercury + 
                         selenium + silver + sodium + zinc,
                         expnms=newXnm,
                         metals, family=gaussian(), q=8, B=20, deg=2)
# underlying fit
summary(qc.fit7a$fit)$coefficients
plot(qc.fit7a)

## ----non-linear examples 2, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
qc.fit7b <- qgcomp.boot(y ~ factor(iron)*factor(lead) + cadmium + 
                         mage35 + arsenic + magnesium + manganese + mercury + 
                         selenium + silver + sodium + zinc,
                         expnms=newXnm,
                         metals, family=gaussian(), q=8, B=10, deg=3)
# underlying fit
#summary(qc.fit7b$fit)$coefficients
plot(qc.fit7b)

## ----non-linear examples 3, results='markup', fig.show='hold', fig.height=3, fig.width=7.5, cache=FALSE----
qc.fit7c <- qgcomp.boot(y ~ I(iron>4)*I(lead>4) + cadmium + 
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

## ----time-to-event1, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
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


## ----time-to-event2, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
# bootstrapped version estimates a marginal structural model for the population average effect
#library(survival)
qc.survfit2 <- qgcomp.cox.boot(Surv(disease_time, disease_state) ~ .,expnms=Xnm,
                         data=metals[,c(Xnm, 'disease_time', 'disease_state')], q=4, 
                         B=5, MCsize=1000, parallel=TRUE, parplan=TRUE)
qc.survfit2

# testing proportional hazards (note that x=TRUE is not needed (and will cause an error if used))
survival::cox.zph(qc.survfit2$fit)
p2 = plot(qc.survfit2, suppressprint = TRUE)  
p2 + labs(title="Linear log(hazard ratio), overall and exposure specific")



## ----time-to-event3, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
qc.survfit3 <- qgcomp.cox.boot(Surv(disease_time, disease_state) ~ . + .^2,expnms=Xnm,
                         data=metals[,c(Xnm, 'disease_time', 'disease_state')], q=4, 
                         B=5, MCsize=1000, parallel=TRUE, parplan=TRUE)
qc.survfit3
p3 = plot(qc.survfit3, suppressprint = TRUE) 
p3 + labs(title="Non-linear log(hazard ratio) overall, linear exposure specific ln-HR")


## ----time-to-event3b, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
future::plan(future::multisession)# parallel evaluation
qc.survfit3 <- qgcomp.cox.boot(Surv(disease_time, disease_state) ~ . + .^2,expnms=Xnm,
                         data=metals[,c(Xnm, 'disease_time', 'disease_state')], q=4, 
                         B=5, MCsize=1000, parallel=TRUE, parplan=FALSE)
qc.survfit3
p3 = plot(qc.survfit3, suppressprint = TRUE) 
p3 + labs(title="Non-linear log(hazard ratio) overall, linear exposure specific ln-HR")


## ----time-to-event4, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
qc.survfit4 <- qgcomp.cox.boot(Surv(disease_time, disease_state) ~ . + .^2,expnms=Xnm,
                         data=metals[,c(Xnm, 'disease_time', 'disease_state')], q=4, 
                         B=5, MCsize=1000, parallel=TRUE, parplan=FALSE, degree=2)
qc.survfit4
# examining the overall hazard ratio as a function of overall exposure
hrs_q = exp(matrix(c(0,0,1,1,2,4,3,9), ncol=2, byrow=TRUE)%*%qc.survfit4$msmfit$coefficients)
colnames(hrs_q) = "Hazard ratio"
print("Hazard ratios by quartiles (min-25%,25-50%, 50-75%, 75%-max)")
hrs_q

p4 = plot(qc.survfit4, suppressprint = TRUE) 
p4 + labs(title="Non-linear log(hazard ratio), overall and exposure specific") 


## ----time-to-event5, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

# testing proportional hazards (must set x=TRUE in function call)
qc.survfit1ph <- qgcomp.cox.noboot(survival::Surv(disease_time, disease_state) ~ .,expnms=Xnm,
                         data=metals[,c(Xnm, 'disease_time', 'disease_state', "mage35")], q=4,
                         x=TRUE)
survival::cox.zph(qc.survfit1ph$fit)


## ----time-to-event6, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

# testing global proportional hazards for model (note that x=TRUE is not needed (and will cause an error if used))
phtest3 = survival::cox.zph(qc.survfit3$fit)
phtest3$table[dim(phtest3$table)[1],, drop=FALSE]


## ----clustering, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
set.seed(2123)
N = 250
t = 4
dat <- data.frame(row.names = 1:(N*t))
dat <- within(dat, {
  id = do.call("c", lapply(1:N, function(x) rep(x, t)))
  u =  do.call("c", lapply(1:N, function(x) rep(runif(1), t)))
  x1 = rnorm(N, u)
  y = rnorm(N) + u + x1
})

# pre-quantize
expnms = c("x1")
datl = quantize(dat, expnms = expnms)

qgcomp.noboot(y~ x1, data=datl$dat, family=gaussian(), q = NULL)

# neither of these ways yields appropriate clustering
#qgcomp.noboot(y~ x1, data=datl$dat, id="id", family=gaussian(), q = NULL)
#qgcomp.boot(y~ x1, data=datl$dat, family=gaussian(), q = NULL, MCsize=1000)

## ----clustering 2, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

# clustering by specifying id parameter on
qgcomp.boot(y~ x1, data=datl$dat, id="id", family=gaussian(), q = NULL, MCsize=1000, B = 5)
#qgcomp.boot(y~ x1, data=datl$dat, id="id", family=gaussian(), q = NULL, MCsize=1000, B = 500)
#   Mixture slope parameters (bootstrap CI):
#   
#               Estimate Std. Error Lower CI Upper CI t value
#   (Intercept)  -0.4632     0.0730   -0.606    -0.32 3.3e-10
#   psi1          0.9550     0.0398    0.877     1.03       0

# This can be verified using the `sandwich` package 
#fitglm = glm(y~x1, data=datl$dat)
#sw.cov = sandwich::vcovCL(fitglm, cluster=~id, type = "HC0")[2,2]
#sqrt(sw.cov)
# [1] 0.0409

## ----pe1, fig.height=5, fig.width=7.5-----------------------------------------
(qc.fit.adj <- qgcomp.noboot(y~.,dat=metals[,c(Xnm, covars, 'y')], expnms=Xnm, family=gaussian()))
plot(qc.fit.adj)

## ----pe2----------------------------------------------------------------------
# 40/60% training/validation split
set.seed(123211)
trainidx <- sample(1:nrow(metals), round(nrow(metals)*0.4))
valididx <- setdiff(1:nrow(metals),trainidx)
traindata <- metals[trainidx,]
validdata <- metals[valididx,]
dim(traindata) # 40% of total
dim(validdata) # 60% of total

## ----pe3a, fig.height=5, fig.width=7.5----------------------------------------
    

splitres <- qgcomp.partials(
  fun="qgcomp.noboot", f=y~., q=4, 
  traindata=traindata[,c(Xnm, covars, "y")],validdata=validdata[,c(Xnm, covars, "y")], expnms=Xnm,
  bayes=FALSE, 
  .fixbreaks = TRUE, .globalbreaks=FALSE
  )
splitres

## ----pe3b, fig.height=5, fig.width=7.5----------------------------------------

plot(splitres$pos.fit)


## ----pe3c, fig.height=5, fig.width=7.5----------------------------------------
    

splitres_alt <- qgcomp.partials(
  fun="qgcomp.noboot", f=y~., q=4, 
  traindata=traindata[,c(Xnm, covars, "y")],validdata=validdata[,c(Xnm, covars, "y")], expnms=Xnm,
  bayes=FALSE, 
  .fixbreaks = TRUE, .globalbreaks=TRUE
  )
splitres_alt

## ----pe4a---------------------------------------------------------------------


nonessentialXnm <- c(
    'arsenic','barium','cadmium','chromium','lead','mercury','silver'
)
essentialXnm <- c(
  'sodium','magnesium','calcium','manganese','iron','copper','zinc','selenium'
)
covars = c('nitrate','nitrite','sulfate','ph', 'total_alkalinity','total_hardness')


(qc.fit.essential <- qgcomp.noboot(y~.,dat=metals[,c(Xnm, covars, 'y')], expnms=essentialXnm, family=gaussian()))


## ----pe4b---------------------------------------------------------------------
(qc.fit.nonessential <- qgcomp.noboot(y~.,dat=metals[,c(Xnm, covars, 'y')], expnms=nonessentialXnm, family=gaussian()))

## ----md1a---------------------------------------------------------------------
Xnm <- c(
    'arsenic','barium','cadmium','calcium','chromium','copper',
    'iron','lead','magnesium','manganese','mercury','selenium','silver',
    'sodium','zinc'
)
covars = c('nitrate','nitrite','sulfate','ph', 'total_alkalinity','total_hardness')
asmiss = metals
set.seed(1232)
asmiss$arsenic = ifelse(runif(nrow(metals))>0.7, NA, asmiss$arsenic)
cc = asmiss[complete.cases(asmiss[,c(Xnm, covars, "y")]),] # complete.cases gives a logical index to subset rows
dim(metals) # [1] 452  26
dim(cc) # [1] 320  26

## ----md1b---------------------------------------------------------------------
qc.base <- qgcomp.noboot(y~.,expnms=Xnm, dat=metals[,c(Xnm, covars, 'y')], family=gaussian())
cat("Full data\n")
qc.base

## ----md1c---------------------------------------------------------------------

qc.cc  <- qgcomp.noboot(y~.,expnms=Xnm, dat=cc[,c(Xnm, covars, 'y')], family=gaussian())
cat("Complete case analyses\n")
cat("  #1 explicitly remove observations with missing values\n")
qc.cc


## ----md1d---------------------------------------------------------------------
qc.cc2 <- qgcomp.noboot(y~.,expnms=Xnm, dat=asmiss[,c(Xnm, covars, 'y')], family=gaussian())



cat("  #1 rely on R handling of NA values\n")
qc.cc2


## ----md1e---------------------------------------------------------------------
# calculation of arsenic quantiles is identical
all.equal(qc.cc$qx$arsenic_q, qc.cc2$qx$arsenic_q[complete.cases(qc.cc2$qx$arsenic_q)])
# all are equal

all.equal(qc.cc$qx$cadmium_q, qc.cc2$qx$cadmium_q[complete.cases(qc.cc2$qx$arsenic_q)])
# not equal

## ----parend, echo=TRUE--------------------------------------------------------
# return to standard processing
future::plan(future::sequential) # return to standard evaluation

