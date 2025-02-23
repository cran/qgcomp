## ----invisibles, echo=FALSE, results='markup', message=FALSE------------------
library("knitr")
#library("gWQS")

## ----first step, echo=TRUE, results='markup', message=FALSE-------------------
library("qgcomp")
library("ggplot2")
library("splines")
data("metals", package="qgcomp")
# using data from the intro vignette
Xnm <- c(
    'arsenic','barium','cadmium','calcium','chromium','copper',
    'iron','lead','magnesium','manganese','mercury','selenium','silver',
    'sodium','zinc'
)
covars = c('nitrate','nitrite','sulfate','ph', 'total_alkalinity','total_hardness')

cormat = cor(metals[,Xnm])

idx = which(cormat>0.6 & cormat <1.0, arr.ind = TRUE)
newXnm = unique(rownames(idx)) # iron, lead, and cadmium



## ----tm2evnt1, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
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


## ----tm2evnt2, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
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



## ----tm2evnt3, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
qc.survfit3 <- qgcomp.cox.boot(Surv(disease_time, disease_state) ~ . + .^2,expnms=Xnm,
                         data=metals[,c(Xnm, 'disease_time', 'disease_state')], q=4, 
                         B=5, MCsize=1000, parallel=TRUE, parplan=TRUE)
qc.survfit3
p3 = plot(qc.survfit3, suppressprint = TRUE) 
p3 + labs(title="Non-linear log(hazard ratio) overall, linear exposure specific ln-HR")


## ----tm2evnt3b, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
future::plan(future::multisession)# parallel evaluation
qc.survfit3 <- qgcomp.cox.boot(Surv(disease_time, disease_state) ~ . + .^2,expnms=Xnm,
                         data=metals[,c(Xnm, 'disease_time', 'disease_state')], q=4, 
                         B=5, MCsize=1000, parallel=TRUE, parplan=FALSE)
qc.survfit3
p3 = plot(qc.survfit3, suppressprint = TRUE) 
p3 + labs(title="Non-linear log(hazard ratio) overall, linear exposure specific ln-HR")


## ----tm2evnt4, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----
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


## ----tm2evnt5, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

# testing proportional hazards (must set x=TRUE in function call)
qc.survfit1ph <- qgcomp.cox.noboot(survival::Surv(disease_time, disease_state) ~ .,expnms=Xnm,
                         data=metals[,c(Xnm, 'disease_time', 'disease_state', "mage35")], q=4,
                         x=TRUE)
survival::cox.zph(qc.survfit1ph$fit)


## ----tm2evnt6, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

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

qgcomp.glm.noboot(y~ x1, data=datl$dat, family=gaussian(), q = NULL)

# neither of these ways yields appropriate clustering
#qgcomp.glm.noboot(y~ x1, data=datl$dat, id="id", family=gaussian(), q = NULL)
#qgcomp.glm.boot(y~ x1, data=datl$dat, family=gaussian(), q = NULL, MCsize=1000)

## ----clustering 2, results='markup', fig.show='hold', fig.height=5, fig.width=7.5, cache=FALSE----

# clustering by specifying id parameter on
#qgcomp.glm.boot(y~ x1, data=datl$dat, id="id", family=gaussian(), q = NULL, MCsize=1000, B = 5)
eefit <- qgcomp.glm.ee(y~ x1, data=datl$dat, id="id", family=gaussian(), q = NULL)
eefit
#qgcomp.glm.boot(y~ x1, data=datl$dat, id="id", family=gaussian(), q = NULL, MCsize=1000, B = 500)
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
(qc.fit.adj <- qgcomp.glm.noboot(y~.,dat=metals[,c(Xnm, covars, 'y')], expnms=Xnm, family=gaussian()))
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
  fun="qgcomp.glm.noboot", f=y~., q=4, 
  traindata=traindata[,c(Xnm, covars, "y")],validdata=validdata[,c(Xnm, covars, "y")], expnms=Xnm,
  bayes=FALSE, 
  .fixbreaks = TRUE, .globalbreaks=FALSE
  )
splitres

## ----pe3b, fig.height=5, fig.width=7.5----------------------------------------

plot(splitres$pos.fit)


## ----pe3c, fig.height=5, fig.width=7.5----------------------------------------
    

splitres_alt <- qgcomp.partials(
  fun="qgcomp.glm.noboot", f=y~., q=4, 
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


(qc.fit.essential <- qgcomp.glm.noboot(y~.,dat=metals[,c(Xnm, covars, 'y')], expnms=essentialXnm, family=gaussian()))


## ----pe4b---------------------------------------------------------------------
(qc.fit.nonessential <- qgcomp.glm.noboot(y~.,dat=metals[,c(Xnm, covars, 'y')], expnms=nonessentialXnm, family=gaussian()))

## ----multinomial, results='markup', fig.show='hold', fig.height=10, fig.width=7.5, cache=FALSE----
data("metals") # from qgcomp package
# create categorical outcome from the existing continuous outcome (usually, one will already exist)
metals$ycat = factor(quantize(metals, "y",q=4)$data$y, levels=c("0", "1", "2", "3"), labels=c("cct", "ccg", "aat", "aag")) 
# restrict to smaller dataset for simplicity
smallmetals = metals[,c("ycat", "arsenic", "lead", "cadmium", "mage35")]

## ----multinomial-2, results='markup', fig.show='hold', fig.height=10, fig.width=7.5, cache=FALSE----

### 1: Define mixture and underlying model ####
mixture = c("arsenic", "lead", "cadmium")
f2 = ycat ~ arsenic + lead + cadmium + mage35

rr = qgcomp.multinomial.noboot(
 f2, 
 expnms = mixture,
 q=4, 
 data = smallmetals, 
 )


rr2 = qgcomp.multinomial.boot(
 f2, 
 expnms = mixture,
 q=4, 
 data = smallmetals, 
 B =2,  # set to higher values >200 in general usage
 MCSize=10000 # set to higher values in small samples
 )

summary(rr)
summary(rr2) # differs from `rr` primarily due to low `MCSize` value

 plot(rr) 
#plot(rr2) # not yet functional

## ----weighting, results='markup', fig.show='hold', fig.height=10, fig.width=7.5, cache=FALSE----

### 1: Define mixture and underlying model ####
set.seed(12321)
metals$samplingweights = exp(rnorm(nrow(smallmetals), -0.5, 1))

uw = qgcomp.glm.noboot(y~., expnms = Xnm,q=4, data = metals[,c(Xnm, covars, "y")])
wtd = qgcomp.glm.noboot(y~., expnms = Xnm,q=4, data = metals[,c(Xnm, covars, "y")], weights=metals$samplingweights)
wtd2 = qgcomp.glm.ee(y~., expnms = Xnm,q=4, data = metals[,c(Xnm, covars, "y")], weights=metals$samplingweights)
#wtd3 = qgcomp.glm.boot(y~., expnms = Xnm,q=4, data = metals[,c(Xnm, covars, "y")], weights=metals$samplingweights)

# unweighted
uw
# weighted with invalid standard error, confidence interval
wtd
# weighted with valid (robust) standard error, confidence interval
wtd2
# weighted with valid (bootstrap) standard error, confidence interval
#wtd3


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
qc.base <- qgcomp.glm.noboot(y~.,expnms=Xnm, dat=metals[,c(Xnm, covars, 'y')], family=gaussian())
cat("Full data\n")
qc.base

## ----md1c---------------------------------------------------------------------

qc.cc  <- qgcomp.glm.noboot(y~.,expnms=Xnm, dat=cc[,c(Xnm, covars, 'y')], family=gaussian())
cat("Complete case analyses\n")
cat("  #1 explicitly remove observations with missing values\n")
qc.cc


## ----md1d---------------------------------------------------------------------
qc.cc2 <- qgcomp.glm.noboot(y~.,expnms=Xnm, dat=asmiss[,c(Xnm, covars, 'y')], family=gaussian())



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

