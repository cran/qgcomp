# functions that implement sample splitting for estimating partial effects and other purposes

qgcomp.partials <- function(
    fun = c("qgcomp.glm.noboot", "qgcomp.cox.noboot", "qgcomp.zi.noboot"),
    traindata=NULL,
    validdata=NULL,
    expnms=NULL,
    .fixbreaks=TRUE,
    .globalbreaks=FALSE,
    ...
){
  #' @title Partial effect sizes, confidence intervals, hypothesis tests
  #' @description Obtain effect estimates for "partial positive" and "partial
  #' negative" effects using quantile g-computation. This approach uses sample
  #' splitting to evaluate the overall impact of a set of variables with 
  #' effects in a single direction, where, using training data, all variables
  #' with effects in the same direction are grouped.  
  #' 
  #' @details In the basic (non bootstrapped) `qgcomp` functions, the positive and 
  #' negative "sums
  #' of coefficients" or "partial effect sizes" are given, which equal the sum
  #' of the negative and positive coefficients in the underlying model. Unfortunately,
  #' these partial effects don't constitute variables for which we can derive confidence
  #' intervals or hypothesis tests, so they are mainly for exploratory purposes. By employing
  #' sample splitting, however, we can obtain better estimates of these partial effects.
  #' 
  #' Sample splitting proceeds by partitioning the data into two samples (40/60 training/validtion
  #' split seems acceptable in many circumstances). The "overall mixture effect" is then
  #' estimated in the training data, and the mixture variables with positive and negative coefficients
  #' are split into separate groups. These two different groups are then used as 
  #' "the mixture of interest" in two additional qgcomp fits, where the mixture of interest
  #' is adjusted for the other exposure variables. For example, if the "positive partial effect"
  #' is of interest, then this effect is equal to the sum of the coefficients in the 
  #' qgcomp model fit to the validation data, with the mixture of interest selected by the
  #' original fit to the training data (note that some of these coefficients may be negative
  #' in the fit to the validation data - this is expected and necessary for valid hypothesis tests).
  #' 
  #' The positive/negative partial effects are necessarily exploratory, but sample splitting preserves
  #' the statistical properties at the expense of wider confidence intervals and larger variances. The
  #' two resulting mixture groups groups should be inspected for 
  #' @param fun character variable in the set "qgcomp.glm.noboot" (binary, count, continuous outcomes), 
  #' "qgcomp.cox.noboot" (survival outcomes), 
  #' "qgcomp.zi.noboot" (zero inflated outcomes). This describes which qgcomp package
  #' function is used to fit the model. (default = "qgcomp.glm.noboot")
  #' @return A 'qgcompmultifit' object, which inherits from \code{\link[base]{list}}, which contains
  #' \describe{
  #' \item{posmix}{character vector of variable names with positive coefficients in the qgcomp model 
  #'   fit to the training data} 
  #' \item{negmix}{character vector of variable names with negative coefficients in the qgcomp model 
  #'   fit to the training data} 
  #' \item{pos.fit}{a qgcompfit object fit to the validation data, in which the exposures of
  #'  interest are contained in 'posmix'} 
  #' \item{neg.fit}{a qgcompfit object fit to the validation data, in which the exposures of
  #'  interest are contained in 'negmix'}
  #'  }
  #' @param traindata Data frame with training data
  #' @param validdata Data frame with validation data
  #' @param expnms Exposure mixture of interest
  #' @param .fixbreaks (logical, overridden by .globalbreaks) Use the same quantile cutpoints in the training and validation data (selected in the training data). As of version 2.8.11, the default is TRUE, whereas it was implicitly FALSE in prior verions. Setting to TRUE increases variance but greatly decreases bias in smaller samples.
  #' @param .globalbreaks (logical, if TRUE, overrides .fixbreaks) Use the same quantile cutpoints in the training and validation data (selected in combined training and validation data). As of version 2.8.11, the default is TRUE, whereas it was implicitly FALSE in prior verions. Setting to TRUE increases variance but greatly decreases bias in smaller samples.
  #' @param ... Arguments to \code{\link[qgcomp]{qgcomp.glm.noboot}}, 
  #'    \code{\link[qgcomp]{qgcomp.cox.noboot}}, or 
  #'    \code{\link[qgcomp]{qgcomp.zi.noboot}}
  #' @family qgcomp_methods
  #' @export
  #' @examples 
  #' set.seed(123223)
  #' dat = qgcomp::simdata_quantized(n=1000, outcomtype="continuous", cor=c(.75, 0), 
  #'                                 b0=0, coef=c(0.25,-0.25,0,0), q=4)
  #' cor(dat)
  #' # overall fit (more or less null due to counteracting exposures)
  #' (overall <- qgcomp.glm.noboot(f=y~., q=NULL, expnms=c("x1", "x2", "x3", "x4"), data=dat))
  #' 
  #' # partial effects using 40% training/60% validation split
  #' trainidx <- sample(1:nrow(dat), round(nrow(dat)*0.4))
  #' valididx <- setdiff(1:nrow(dat),trainidx)
  #' traindata = dat[trainidx,]
  #' validdata = dat[valididx,]
  #' splitres <- qgcomp.partials(fun="qgcomp.glm.noboot", f=y~., q=NULL, 
  #'     traindata=traindata,validdata=validdata, expnms=c("x1", "x2", "x3", "x4"))
  #' splitres
  #' \dontrun{
  #' # under the null, both should give null results
  #' set.seed(123223)
  #' dat = simdata_quantized(n=1000, outcomtype="continuous", cor=c(.75, 0), 
  #'                                 b0=0, coef=c(0,0,0,0), q=4)
  #' # 40% training/60% validation
  #' trainidx2 <- sample(1:nrow(dat), round(nrow(dat)*0.4))
  #' valididx2 <- setdiff(1:nrow(dat),trainidx2)
  #' traindata2 <- dat[trainidx2,]
  #' validdata2 <- dat[valididx2,]
  #' splitres2 <- qgcomp.partials(fun="qgcomp.glm.noboot", f=y~., 
  #'    q=NULL, traindata=traindata2,validdata=validdata2, expnms=c("x1", "x2", "x3", "x4"))
  #' splitres2
  #' 
  #' # 60% training/40% validation
  #' trainidx3 <- sample(1:nrow(dat), round(nrow(dat)*0.6))
  #' valididx3 <- setdiff(1:nrow(dat),trainidx3)
  #' traindata3 <- dat[trainidx3,]
  #' validdata3 <- dat[valididx3,]
  #' splitres3 <- qgcomp.partials(fun="qgcomp.glm.noboot", f=y~., q=NULL, 
  #'     traindata=traindata3,validdata=validdata3, expnms=c("x1", "x2", "x3", "x4"))
  #' splitres3
  #' 
  #' # survival outcome
  #' set.seed(50)
  #' N=1000
  #' dat = simdata_quantized(n=1000, outcomtype="survival", cor=c(.75, 0, 0, 0, 1), 
  #'                                 b0=0, coef=c(1,0,0,0,0,1), q=4)
  #' names(dat)[which(names(dat=="x5"))] = "z"
  #' trainidx4 <- sample(1:nrow(dat), round(nrow(dat)*0.6))
  #' valididx4 <- setdiff(1:nrow(dat),trainidx4)
  #' traindata4 <- dat[trainidx4,]
  #' validdata4 <- dat[valididx4,]
  #' expnms=paste0("x", 1:5)
  #' f = survival::Surv(time, d)~x1 + x2 + x3 + x4 + x5 + z
  #' (fit1 <- survival::coxph(f, data = dat))
  #' (overall <- qgcomp.cox.noboot(f, expnms = expnms, data = dat))
  #' (splitres4 <- qgcomp.partials(fun="qgcomp.cox.noboot", f=f, q=4,
  #'                traindata=traindata4,validdata=validdata4,
  #'                 expnms=expnms))
  #'  
  #'  # zero inflated count outcome
  #'  set.seed(50)
  #'  n=1000
  #'  dat <- data.frame(y= (yany <- rbinom(n, 1, 0.5))*(ycnt <- rpois(n, 1.2)), x1=runif(n)+ycnt*0.2, 
  #'                        x2=runif(n)-ycnt*0.2, x3=runif(n),
  #'                       x4=runif(n) , z=runif(n))
  #'  # poisson count model, mixture in both portions, but note that the qgcomp.partials
  #'  # function defines the "positive" variables only  by the count portion of the model
  #'  (overall5 <- qgcomp.zi.noboot(f=y ~ z + x1 + x2 + x3 + x4 | x1 + x2 + x3 + x4 + z, 
  #'                  expnms = c("x1", "x2", "x3", "x4"), 
  #'                   data=dat, q=4, dist="poisson"))
  #'
  #' trainidx5 <- sample(1:nrow(dat), round(nrow(dat)*0.6))
  #' valididx5 <- setdiff(1:nrow(dat),trainidx5)
  #' traindata5 <- dat[trainidx5,]
  #' validdata5 <- dat[valididx5,]
  #' splitres5 <- qgcomp.partials(fun="qgcomp.zi.noboot", 
  #'     f=y ~ x1 + x2 + x3 + x4 + z | x1 + x2 + x3 + x4 + z, q=4, 
  #'     traindata=traindata5, validdata=validdata5, 
  #'     expnms=c("x1", "x2", "x3", "x4"))
  #' splitres5
  #'                 
  #' }
  if(is.null(traindata) | is.null(validdata))
    stop("traindata and validdata must both be specified")
  #
  traincall <- validcall <- match.call(expand.dots = TRUE)
  droppers <- match(c("traindata", "validdata", ".fixbreaks", ".globalbreaks", "fun"), names(traincall), 0L) #index (will need to add names here if more arguments are added)
  traincall[["data"]] <- eval(traincall[["traindata"]], parent.frame())
  validcall[["data"]] <- eval(validcall[["validdata"]], parent.frame())
  traincall <- traincall[-c(droppers)]
  validcall <- validcall[-c(droppers)]
  hasbreaks = ifelse("breaks" %in% names(traincall), TRUE, FALSE)
  hasq = ifelse("q" %in% names(traincall), TRUE, FALSE)
  qnull = ifelse(is.null(traincall$q), TRUE, FALSE)
  if(hasbreaks && .fixbreaks)
    .fixbreaks=FALSE
  # if q is set to null, and no breaks are provided ensure that no breaks are created
  if(hasq){
    if(!hasbreaks && qnull && .fixbreaks)
      .fixbreaks=FALSE
    if(!hasbreaks && qnull && .globalbreaks)
      .globalbreaks=FALSE
  }
  
  if(is.function(fun)){
    traincall[[1L]] <- validcall[[1L]] <- fun
  }else{
    traincall[[1L]] <- validcall[[1L]] <- as.name(fun[1])
  }
  if(.globalbreaks){
    #
    globalcall <- traincall
    globalcall[["data"]] <- eval(rbind(traincall[["data"]], validcall[["data"]]), parent.frame())
    global.fit = eval(globalcall, parent.frame())
    #print(dim(global.fit$fit$data))
    #
    if(!is.null(global.fit$breaks)){
      validcall[["breaks"]] = global.fit$breaks
      validcall[["q"]] = NULL
      traincall[["breaks"]] = global.fit$breaks
      traincall[["q"]] = NULL
    }
  }
  train.fit = eval(traincall, parent.frame())
  #####
  if(.fixbreaks && !.globalbreaks){
    validcall$breaks = train.fit$breaks
    validcall$q = NULL
  }
  ######
  posnms = expnms[is.element(expnms, names(train.fit$pos.weights))]
  negnms = expnms[is.element(expnms, names(train.fit$neg.weights))]
  if(length(posnms)==1 && all(posnms==c("count", "zero"))){
    posnms = names(train.fit$pos.weights$count)
    negnms = names(train.fit$neg.weights$count)
  }
  res = list(train.fit=train.fit)
  res$negmix <- res$posmix <- "none"
  if(length(posnms)>0){
    res$posmix = posnms
    poscall <- validcall
    if(!is.null(poscall$breaks)  && (.fixbreaks || .globalbreaks)){
      posidx <- which(expnms %in% posnms)
      poscall$breaks <- poscall$breaks[posidx]
    }
    vc = as.list(poscall)
    vc$expnms = c(posnms)
    res$pos.fit <- eval(as.call(vc), parent.frame())
  }
  if(length(negnms)>0){
    res$negmix = negnms
    negcall <- validcall
    if(!is.null(negcall$breaks)  && (.fixbreaks || .globalbreaks)){
      #if(.fixbreaks || .globalbreaks){
      negidx <- which(expnms %in% negnms)
      negcall$breaks <- negcall$breaks[negidx]
    }
    vc = as.list(negcall)
    vc$expnms = c(negnms)
    res$neg.fit <- eval(as.call(vc), parent.frame())
    
  }
  class(res) <- "qgcompmultifit"
  res
}





print.qgcompmultifit <- function(x, ...){
  #' @export
  cat(paste0("\nVariables with positive effect sizes in training data: ", paste(x$posmix, collapse = ", ")))
  cat(paste0("\nVariables with negative effect sizes in training data: ", paste(x$negmix, collapse = ", ")))
  cat("\nPartial effect sizes estimated in validation data\n")
  cat("Positive direction ")
  if(!is.null(x$pos.fit)){
    print(x$pos.fit, showweights=FALSE)
  } else{
    cat("\n No positive coefficients in model fit to training data ")
  }
  cat("\nNegative direction ")
  if(!is.null(x$neg.fit)){
    print(x$neg.fit, showweights=FALSE)
  } else{
    cat("\n No negative coefficients in model fit to training data ")
  } 
}



split_data <- function(data,  cluster=NULL, prop.train=0.4){
  #' @title Perform sample splitting
  #' @description This is a convenience function to split the input data into 
  #'  two independent sets, possibly accounting for 
  #'  single level clustering. These two sets can be used with 
  #'  \code{\link[qgcomp]{qgcomp.partials}} to get "partial" positive/negative effect estimates
  #'  from the original data, where sample splitting is necessary to get valid confidence intervals
  #'  and p-values. Sample splitting is also useful for any sort of exploratory model selection, where
  #'  the training data can be used to select the model and the validation model used to 
  #'  generate the final estimates (this process should not be iterative - e.g. no "checking" the 
  #'  results in the validation data and then re-fitting, as this invalidates inference in the
  #'  validation set.) E.g. you could use the training data to select non-linear terms for the 
  #'  model and then re-fit in validation data to get unbiased estimates.
  #' @param data A data.frame for use in qgcomp fitting
  #' @param cluster NULL (default) or character value naming a cluster identifier in the data.
  #' This is to prevent observations from a single cluster being in both the training and
  #' validation data, which reduces the effectiveness of sample splitting.
  #' @param prop.train proportion of the original dataset (or proportion of the clusters 
  #' identified via the 'cluster' parameter) that are used in the training data (default=0.4)
  #'
  #' @return A list of the following type:
  #'   list(
  #'    trainidx = trainidx,
  #'    valididx = valididx,
  #'    traindata = traindata,
  #'    validdata = validdata
  #'   )
  #'   
  #'   e.g. if you call `spl = split_data(dat)`, then spl$traindata will contain
  #'   a 40% sample from the original data, spl$validdata will contain the other 60%
  #'   and spl$trainidx, spl$valididx will contain integer indexes that track the 
  #'   row numbers (from the original data `dat`) that have the training and validation
  #'   samples.
  #'   

  #' @export
  #'
  #' @examples 
  #' data(metals)
  #' set.seed(1231124)
  #' spl = split_data(metals)
  #' Xnm <- c(
  #'   'arsenic','barium','cadmium','calcium','chromium','copper',
  #'  'iron','lead','magnesium','manganese','mercury','selenium','silver',
  #'  'sodium','zinc'
  #' )
  #' dim(spl$traindata) # 181 observations = 40% of total
  #' dim(spl$validdata) # 271 observations = 60% of total
  #' splitres <- qgcomp.partials(fun="qgcomp.glm.noboot", f=y~., q=4, 
  #'   traindata=spl$traindata,validdata=spl$validdata, expnms=Xnm)
  #' splitres
  #' 
  #' # also used to compare linear vs. non-linear fits (useful if you have enough data)
  #' set.seed(1231)
  #' spl = split_data(metals, prop.train=.5)
  #' lin = qgcomp.glm.boot(f=y~., q=4, expnms=Xnm, B=5, data=spl$traindata)
  #' nlin1 = qgcomp.glm.boot(f=y~. + I(manganese^2) + I(calcium^2), expnms=Xnm, deg=2, 
  #'   q=4, B=5, data=spl$traindata)
  #' nlin2 = qgcomp.glm.boot(f=y~. + I(arsenic^2) + I(cadmium^2), expnms=Xnm, deg=2, 
  #'   q=4, B=5, data=spl$traindata)
  #' AIC(lin);AIC(nlin1);AIC(nlin2)
  #' # linear has lowest training AIC, so base final fit off that (and bootstrap not needed)
  #' qgcomp.glm.noboot(f=y~., q=4, expnms=Xnm, data=spl$validdata)
  
  if(is.null(cluster[1])){
    ret = .split.iid.data(data, prop.train=prop.train)
  }else{
    ret = .split.cluster.data(data, cluster=cluster, prop.train=prop.train)
  }
  ret
}

# experimental functions that may make it into permanent base files

.split.iid.data <- function(data, prop.train=0.4){
  trainidx <- sort(sample(seq_len(nrow(data)), round(nrow(data)*prop.train)))
  valididx <- setdiff(seq_len(nrow(data)),trainidx)
  traindata <- data[trainidx,]
  validdata <- data[valididx,]
  list(
    trainidx = trainidx,
    valididx = valididx,
    traindata = traindata,
    validdata = validdata
  )
}

.split.cluster.data <- function(data, cluster="id", prop.train=0.4){
  clustids = sort(unique(data[[cluster]]))
  trainidx = vector()
  for(clust in clustids){
    inclust = which(data[[cluster]] == clust)
    trainidx = c(trainidx, sort(sample(inclust, round(length(inclust)*prop.train))))
  }
  valididx <- setdiff(seq_len(nrow(data)),trainidx)
  traindata <- data[trainidx,]
  validdata <- data[valididx,]
  list(
    trainidx = trainidx,
    valididx = valididx,
    traindata = traindata,
    validdata = validdata
  )
}
