# generic functions, http://4dpiecharts.com/2011/01/10/really-useful-bits-of-code-that-are-missing-from-r/
geomean <- function(x, na.rm = FALSE, trim = 0, ...)
{
  exp(mean(log(x, ...), na.rm = na.rm, trim = trim, ...))
}

geosd <- function(x, na.rm = FALSE, ...)
{
  exp(sd(log(x, ...), na.rm = na.rm, ...))
}

# Yield response functions
## Mitscherlich
mitschFunc <- function(dose, asymptote, taux, environnement, rateExp = 'exp')  {
  if (rateExp == 'exp') {
    asymptote * (1 - exp(-exp(taux) * (dose + environnement)))
  } else {
    asymptote * (1 - exp(-taux * (dose + environnement)))
  }
} 

## Colwell (p. 39)
colwFunc <- function(dose, b0, b1, b2) b0 + b1 * sqrt(dose) + b2 * dose

## Polynomial quadratic
quadFunc <- function(dose, b0, b1, b2) b0 + b1 * dose + b2 * dose**2


# Mixed models
## Mitscherlich
mitscherlich_nlme <- function(constrain.asym = NULL, constrain.env = NULL, constrain.rate = NULL,
                              startVal, startVector, rhs, data,
                              rateExp = 'exp', # should rate be exp-transformed or ordinary?
                              scaleResp = 'ordinary', # should response be log-transformed or ordinary?
                              scaleDose = 'ordinary', # should dose be log-transformed or ordinary?
                              maxIter = 50, minScale = 1e-8,
                              msVerbose = FALSE) {
  
  if (scaleResp == "log10") data$Response <- log10(data$Response)
  if (scaleDose == "log10") {
    data$Dose[data$Dose == 0] <- 0.001
    data$Dose <- log10(data$Dose)
  }
  
  # Build formula according to specified constrains
  if (is.null(constrain.asym)) {
    asymTag <- 'Asymptote'
  } else {
    asymTag <- paste(constrain.asym[1], '+', constrain.asym[2], '/(1+exp(-Asymptote))')
  }
  
  if (is.null(constrain.env)) {
    envTag = 'Environment'
  } else {
    envTag <- paste(constrain.env[1], '+', constrain.env[2], '/(1+exp(-Environment))')
  }
  
  if (is.null(constrain.rate)) {
    rateTag = 'Rate' # default is exponential
  } else {
    rateTag <- paste(constrain.rate[1], '+', constrain.rate[2], '/(1+exp(-Rate))')
  }
  
  if (rateExp == 'exp') {
    mainMod <- as.formula(paste('Response ~', asymTag, '* ((1 - exp( -exp(', rateTag, ') * (Dose +', envTag, '))))'))
  } else {
    mainMod <- as.formula(paste('Response ~', asymTag, '* ((1 - exp( -', rateTag, '* (Dose +', envTag, '))))'))
  }
  
  mm <- nlme(mainMod,
             data = data, 
             start = c(Asymptote = startVal[1], startVector,
                       Rate = startVal[2], startVector,
                       Environment = startVal[3], startVector 
             ), 
             fixed = list(as.formula(paste("Asymptote ~ ", rhs)),
                          as.formula(paste("Rate ~ ", rhs)),
                          as.formula(paste("Environment ~ ", rhs))
             ), 
             random = Asymptote ~ 1 | NoEssai/NoBloc,
             control = list(maxIter = maxIter, returnObject=TRUE, 
                            msVerbose=msVerbose, minScale = minScale),
             method = 'REML')
  
  return(mm)
}

## Colwell
colwell_nlme <- function(startVal, startVector, rhs, data,
                         scaleResp = 'ordinary', maxIter = 50, minScale = 1e-8,
                         msVerbose = FALSE) {
  
  if (scaleResp == "log10") data$Response <- log10(data$Response)
  
  mm <- nlme(Response ~  b0 + b1 * sqrt(Dose) + b2 * Dose,
             data = data, 
             start = c(b0 = startVal[1], startVector,
                       b1 = startVal[2], startVector,
                       b2 = startVal[3], startVector
             ), 
             fixed = list(as.formula(paste("b0 ~ ", rhs)),
                          as.formula(paste("b1 ~ ", rhs)),
                          as.formula(paste("b2 ~ ", rhs))
             ), 
             random = b1 ~ 1 | NoEssai/NoBloc,
             control = list(maxIter = maxIter, returnObject=TRUE, 
                            msVerbose=msVerbose, minScale = minScale))
  
  return(mm)
}

## Polynomial quadratic
quad_nlme <- function(startVal, startVector, rhs, data,
                         scaleResp = 'ordinary', maxIter = 50, minScale = 1e-8,
                         msVerbose = FALSE) {
  
  if (scaleResp == "log10") data$Response <- log10(data$Response)
  
  mm <- nlme(Response ~  b0 + b1 * Dose + b2 * Dose^2,
             data = data, 
             start = c(b0 = startVal[1], startVector,
                       b1 = startVal[2], startVector,
                       b2 = startVal[3], startVector
             ), 
             fixed = list(as.formula(paste("b0 ~ ", rhs)),
                          as.formula(paste("b1 ~ ", rhs)),
                          as.formula(paste("b2 ~ ", rhs))
             ), 
             random = b1 ~ 1 | NoEssai/NoBloc,
             control = list(maxIter = maxIter, returnObject=TRUE, 
                            msVerbose=msVerbose, minScale = minScale))
  return(mm)
}



# Yield response prediction from new data, by hand
## Mitscherlich prediction
mitschPred <- function(mm, newdata, rhs, col_dose, 
                       constrain.asym = NULL, constrain.env = NULL, constrain.rate = NULL, 
                       logYield = FALSE, ranEf = 0, rateExp = 'exp') {
  library(stringi)
  # parameters
  ## collect parameters from the model
  parameter <- summary(mm)$tTable
  
  ## collect specific mean and standard error parameters for Asymptote, Taux and Environnement
  asymParam <- parameter[stri_detect_fixed(rownames(parameter), "Asymptote"), 1:2]
  tauxParam <- parameter[stri_detect_fixed(rownames(parameter), "Rate"), 1:2] 
  enviParam <- parameter[stri_detect_fixed(rownames(parameter), "Enviro"), 1:2]
  
  ## the linear combination of site specific condition (model.matrix) used to compute Mitscherlich parameters
  modmat <- model.matrix(as.formula(paste('~', paste(rhs, collapse = '+'))), data = newdata)
  asymMitsch <- modmat %*% asymParam[, 1]
  tauxMitsch <- modmat %*% tauxParam[, 1]
  enviMitsch <- modmat %*% enviParam[, 1]
  
  # predict
  doseMitsch <- newdata[col_dose]
  
  # Build formula according to specified constrains
  if (is.null(constrain.asym)) {
    asymTag <- asymMitsch[1] + ranEf
  } else {
    asymTag <- constrain.asym[1] + constrain.asym[2]/(1+exp(-(asymMitsch[1] + ranEf)))
  }
  
  if (is.null(constrain.env)) {
    envTag <- enviMitsch[1]
  } else {
    envTag <- constrain.env[1] + constrain.env[2]/(1+exp(-(enviMitsch[1])))
  }
  
  if (is.null(constrain.rate)) {
    rateTag <- tauxMitsch[1]
  } else {
    rateTag <- constrain.rate[1] + constrain.rate[2]/(1+exp(-(tauxMitsch[1])))
  }
  
  pred <- mitschFunc(doseMitsch, asymTag, rateTag, envTag, rateExp = rateExp)
  
  if(logYield) pred <- 10^perd
  
  # output
  return(list(pred = pred[, 1], 
              metaParam = list(A = asymMitsch[1], R = tauxMitsch[1], E = enviMitsch[1]),
              drag = list(A = t(t(modmat[1, ] * asymParam[, 1])),
                          R = t(t(modmat[1, ] * tauxParam[, 1])),
                          E = t(t(modmat[1, ] * enviParam[, 1])))))
}

## Colwell prediction
colwPred <- function(mm, newdata, rhs, col_dose) {
  library(stringi)
  # parameters
  ## collect parameters from the model
  parameter <- summary(mm)$tTable
  
  ## collect specific mean and standard error parameters for Asymptote, Taux and Environnement
  b0Param <- parameter[stri_detect_fixed(rownames(parameter), "b0"), 1:2]
  b1Param <- parameter[stri_detect_fixed(rownames(parameter), "b1"), 1:2] 
  b2Param <- parameter[stri_detect_fixed(rownames(parameter), "b2"), 1:2]
  
  ## the linear combination of site specific condition (model.matrix) used to compute Mitscherlich parameters
  modmat <- model.matrix(as.formula(paste('~', paste(rhs, collapse = '+'))), data = newdata)
  b0Colw <- modmat %*% b0Param[, 1]
  b1Colw <- modmat %*% b1Param[, 1]
  b2Colw <- modmat %*% b2Param[, 1]
  
  # predict
  doseColw <- newData[col_dose]
  pred <- colwFunc(doseColw, b0Colw, b1Colw, b2Colw)
  
  # output
  return(pred[, 1])
}


# Regression analysis
## R-squared
rsq <- function(y, y_hat) {
  sum((y_hat - mean(y))^2) / sum((y - mean(y))^2)
}

## Residuals
residu <- function(model, level = 0, computeResidues = TRUE, residues,
                   main = "") {
  if (computeResidues) {
    r = residuals(model, type="normalized", level)
  } else {
    r = residues
  }
  hist(r, freq=F, main = main)
  xfit<-seq(min(r),max(r),length=40)
  yfit<-dnorm(xfit)
  lines(xfit, yfit, col="red", lwd=2)
  print("Shapiro-Wilk:")
  print(shapiro.test(r))
  print(paste("Kurtosis:", kurtosis(r)))
  print(paste("Skewness:", skewness(r)))
}

## Normality
normality <- function(r) {
  hist(r,freq=F)
  xfit<-seq(min(r),max(r),length=40)
  yfit<-dnorm(xfit)
  lines(xfit, yfit, col="red", lwd=2)
  print("Shapiro-Wilk")
  print(shapiro.test(r))
  print("Kurtosis")
  print(kurtosis(r))
}

# Optimal dose
## Mitscherlich
doseOptMitsherlich <- function(asymptote, taux, environnement, 
                               prix_dose, prix_vendable, interets = 0, periode = 0, 
                               rateExp = 'exp') {
  if (rateExp == 'exp') {
    -exp(-taux) * (log(prix_dose * (1 + interets)**periode / (asymptote * prix_vendable)) + exp(taux) * environnement - taux)
  } else {
    (log((asymptote * taux * prix_vendable) / (prix_dose * (interets + 1)**periode)) - environnement * taux) / taux 
  }
}
## Colwell
doseOptColwell <- function(b0, b1, b2, prix_dose, prix_vendable, interets = 0, periode = 0) {
  En <- prix_dose * (1 + interets)**periode / prix_vendable
  (0.5 * b1 / (En - b2))**2
}

# Compute the sqrt of a matrix
## needed for varImp in nlmePlotConfTable
## http://realizationsinbiostatistics.blogspot.ca/2008/08/matrix-square-roots-in-r_18.html
denman.beavers <- function(mat, maxit = 50) {
  stopifnot(nrow(mat) == ncol(mat))
  niter <- 0
  y <- mat
  z <- diag(rep(1, nrow(mat)))
  for (niter in 1:maxit) {
    y.temp <- 0.5 * (y + solve(z))
    z <- 0.5 * (z + solve(y))
    y <- y.temp
  }
  return(list(sqrt = y, sqrt.inv = z))
}

# Plot fixed effects and their confidence intervals and 
# importance of variables using t-values
nlmePlotConfTable <- function(mm, varCat, varCatOrder, paramCat,
                              conf_level = 0.95, 
                              pval.breaks = c(0, 0.05, 0.10, 1.00),
                              numApart = FALSE,
                              varNamesFrom = NULL, varNamesTo = NULL,
                              varCatNamesTo = NULL, limits_df = NULL,
                              removeIntercept = FALSE,
                              varImpType = "t", # "t" or "cat", t is t-value, cat is correlation-adjusted t-scores
                              greyEnd = 0.7,
                              modify.param.names = NULL) {
  
  # Base table
  confTable <- as.data.frame(intervals(mm, level = conf_level)$fixed)
  
  # p-values
  confTable$p <- summary(mm)$tTable[, 5]
  confTable$pCat <- cut(confTable$p, breaks = pval.breaks)
  confTable$pCat <- factor(sub(',', ' - ', confTable$pCat)) # prettier format
  confTable$test <- as.factor(ifelse(confTable[, 1] <= 0 & confTable[, 3] >= 0, "not significant", "significant"))
  
  # t-values without (t) or with (cat) adjustment for correlation
  if (varImpType == "t") {
    confTable$t <- summary(mm)$tTable[, 4]
    varImpXlab <- "t-value"
  } else if (varImpType == "cat") {
    # http://arxiv.org/pdf/0903.2003v4.pdf, eq 3.1, p.5
    # denman.beavers, custom function (defined) in this file to compute the sqrt of a matrix
    confTable$t <- denman.beavers(cov2cor(vcov(mm)))$sqrt.inv %*% summary(mm)$tTable[, 4]
    varImpXlab <- "Correlation-adjusted t-scores"
  }
  
  # column for Variable category
  confTable$VarLong <- rownames(confTable) # to keep track of original names
  confTable$VarCat <- 'Numeric'
  for (i in 1:length(varCat)) {
    confTable$VarCat[grepl(varCat[i], confTable$VarLong)] <- varCatNamesTo[i]
  }
  rownames(confTable) <- 1:nrow(confTable)
  confTable$VarCat <- as.factor(confTable$VarCat)
  confTable$VarCat <- factor(confTable$VarCat, levels = levels(confTable$VarCat)[varCatOrder])
  
  # Variable names
  confTable$Variable <- confTable$VarLong
  
  ## remove category names from variable names which are categories
  for (i in 1:length(varCat)) {
    confTable$Variable <- sub(varCat[i], '', confTable$Variable)
  }
  
  ## translate names
  if(!is.null(varNamesFrom)) {
    for (i in 1:length(varNamesFrom)) {
      confTable$Variable <- sub(varNamesFrom[i], varNamesTo[i], confTable$Variable)
    }
  }
  
  ## remove parameter names from variable names
  for (i in 1:length(paramCat)) {
    confTable$Variable <- sub(paste0(paramCat[i], '.'), '', confTable$Variable)
  }
  confTable$Variable <- factor(confTable$Variable)
  
  # column for Parameter category
  confTable$ParamCat <- NA
  for (i in 1:length(paramCat)) {
    filterMPC <- grepl(paramCat[i], confTable$VarLong)
    confTable$ParamCat[filterMPC] <- paramCat[i]
  }
  confTable$ParamCat <- as.factor(confTable$ParamCat)
  
  ## change parameter names if specified
  if (!is.null(modify.param.names)) {
    levels(confTable$ParamCat) <- modify.param.names
  }
  
  if(!is.null(limits_df)) {
    if(!is.na(limits_df[1, 2])) confTable$lower[confTable$ParamCat == "Asymptote" & confTable$lower < limits_df[1, 2]] <- NA
    if(!is.na(limits_df[1, 3])) confTable$upper[confTable$ParamCat == "Asymptote" & confTable$upper > limits_df[1, 3]] <- NA
    if(!is.na(limits_df[2, 2])) confTable$lower[confTable$ParamCat == "Environment" & confTable$lower < limits_df[2, 2]] <- NA
    if(!is.na(limits_df[2, 3])) confTable$upper[confTable$ParamCat == "Environment" & confTable$upper > limits_df[2, 3]] <- NA
    if(!is.na(limits_df[3, 2])) confTable$lower[confTable$ParamCat == "Rate" & confTable$lower < limits_df[3, 2]] <- NA
    if(!is.na(limits_df[3, 3])) confTable$upper[confTable$ParamCat == "Rate" & confTable$upper > limits_df[3, 3]] <- NA
  }
  
  if(removeIntercept) {
    confTable <- confTable[confTable$VarCat != 'Intercept', ]
    confTable <- confTable[confTable$Variable != '(Intercept)', ] # temporary hack
  }
  
  if(numApart) {
    fixefPar <- list()
    numPar <- droplevels(confTable[confTable$VarCat %in% c("Numeric"), ]) # "Intercept", 
    facPar <- droplevels(confTable[!(confTable$VarCat %in% c("Numeric")), ])
    levels(facPar$pCat) <- levels(confTable$pCat) # to assure that all p-values levels are in the legend
    levels(numPar$pCat) <- levels(confTable$pCat) # to assure that all p-values levels are in the legend
    fixefPar[[1]] <- ggplot(numPar, aes(y = Variable, x = est., xmin = lower, xmax = upper)) +
      facet_grid(VarCat ~ ParamCat, scale = 'free', space = 'free_y') +
      geom_vline(xintercept=0, linetype = 3) +
      geom_errorbarh(aes(colour=pCat), height = 0, size = 0.6) +
      geom_point(aes(colour=pCat), size = 1, shape = 16) + # , shape = test
      scale_colour_grey(start = 0, end = greyEnd, name = 'p-value', drop = FALSE) +
      xlab('Coefficient') +
      ylab('') +
      theme_bw() +
      theme(strip.text.y = element_text(angle = 0))
    fixefPar[[2]] <- ggplot(facPar, aes(y = Variable, x = est., xmin = lower, xmax = upper)) +
      facet_grid(VarCat ~ ParamCat, scale = 'free', space = 'free_y') +
      geom_vline(xintercept=0, linetype = 3) +
      geom_errorbarh(aes(colour=pCat), height = 0, size = 1) +
      geom_point(aes(colour=pCat), size = 2, shape = 16) + # , shape = test
      scale_colour_grey(start = 0, end = greyEnd, name = 'p-value', drop = FALSE) +
      xlab('Coefficient') +
      ylab('') +
      theme_bw() +
      theme(strip.text.y = element_text(angle = 0))
  } else {
    fixefPar <- ggplot(confTable, aes(y = Variable, x = est., xmin = lower, xmax = upper)) +
      geom_vline(xintercept=0, linetype = 3) +
      geom_errorbarh(aes(colour=pCat), height = 0, size = 1) +
      geom_point(aes(colour=pCat), size = 2, shape = 16) + # , shape = test
      scale_colour_grey(start = 0, end = greyEnd, name = 'p-value') +
      xlab('Coefficient') +
      ylab('') +
      theme_bw() +
      theme(strip.text.y = element_text(angle = 0))
    
    if(length(unique(confTable$ParamCat)) == 1) {
      fixefPar <- fixefPar + facet_grid(. ~ ParamCat, scale = 'free', space = 'free_y')
    } else {
      fixefPar <- fixefPar + facet_grid(VarCat ~ ParamCat, scale = 'free', space = 'free_y')
    }
    
  }
  
  varImp <- ggplot(confTable, aes(x = abs(t), y = Variable)) +
    geom_segment(aes(x = 0, xend=abs(t), yend = Variable,
                     colour=pCat), lwd = 2) +
    scale_colour_grey(start = 0, end = greyEnd, name = 'p-value') +
    xlab(varImpXlab) +
    ylab('') +
    theme_bw() +
    theme(strip.text.y = element_text(angle = 0))
  
  if(length(unique(confTable$VarCat)) == 1) {
    varImp <- varImp + facet_grid(. ~ ParamCat, scale = 'free', space = 'free_y')
  } else {
    varImp <- varImp + facet_grid(VarCat ~ ParamCat, scale = 'free', space = 'free_y')
  }
  
  return(list(fixefPar, varImp))
}


optDoseMitschMC <- function(mm, rhs, col_dose = 'Dose',
                            constrain.asym = NULL, 
                            constrain.env = NULL, 
                            constrain.rate = NULL,
                            logYield = FALSE,
                            refTest = NULL, data = NULL, 
                            simCase = NULL, nbNd = 50, doseMax = 250,
                            nbMC = 50,
                            C = 1.20, V = 275,
                            alphaMC = 0.05,
                            ranEfLevel = 1, rateExp = 'exp') {
  
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(mvtnorm)
  library(stringi)
  
  if(!is.null(refTest)) {
    data <- tbl_df(data) # data must be a tbl_df
    essaiBase <- data %>% filter(NoEssai == refTest)
    col_class <- lapply(data, class)
    casBase <- essaiBase[1, ]
    filterBase <- col_class == "numeric" & !(names(casBase) %in% c("Response", col_dose))
    casBase[, filterBase] <- apply(essaiBase[, filterBase], 2, mean)
    newData <- casBase
    for (i in 2:nbNd) {
      newData <- rbind(newData, casBase)
    }
    newData$Dose <- seq(0, doseMax, length = nbNd)
    
    if (ranEfLevel == 0) {
      ranEf <- 0
    } else if (ranEfLevel == 1) {
      ranEf <- ranef(mm)[[1]][rownames(ranef(mm)[[1]]) == refTest, ]
    } else if (ranEfLevel > 1) {
      stop('ranEfLevel can not have level 2 or more effect')
    }
  }
  
  if(!is.null(simCase)) {
    casBase <- simCase[1, ]
    newData <- casBase
    for (i in 2:nbNd) {
      newData <- rbind(newData, casBase)
    }
    newData$Dose <- seq(0, doseMax, length = nbNd)
    ranEf <- 0
  }
  
  predTB <- mitschPred(mm, newData, rhs, col_dose = 'Dose', 
                       constrain.asym = constrain.asym, 
                       constrain.env = constrain.env, 
                       constrain.rate = constrain.rate, 
                       logYield = FALSE, ranEf = ranEf,
                       rateExp = rateExp)
  
  coefficients <- summary(mm)$tTable
  
  ### collect specific mean and standard error coefficients for Asymptote, Taux and Environnement
  asymParam <- coefficients[stri_detect_fixed(rownames(coefficients), "Asymptote"), 1:2]
  tauxParam <- coefficients[stri_detect_fixed(rownames(coefficients), "Rate"), 1:2] 
  enviParam <- coefficients[stri_detect_fixed(rownames(coefficients), "Enviro"), 1:2]
  
  ## Monte-Carlo simulation
  ### randomly generate coefficients from their estimated mean and standard error
  ### then predict yield and optimal economical dose
  modmat <- model.matrix(as.formula(paste('~', rhs)), data = newData)
  predMC_list <- list()
  
  asymRows <- stri_detect_fixed(rownames(coefficients), "Asymptote")
  tauxRows <- stri_detect_fixed(rownames(coefficients), "Rate")
  enviRows <- stri_detect_fixed(rownames(coefficients), "Enviro")
  
  coefficientsMC <- matrix(nrow = nbMC, ncol = nrow(coefficients))
  colnames(coefficientsMC) <- rownames(coefficients)
  mitschMC <-  matrix(nrow = nbMC, ncol = 5)
  colnames(mitschMC) <- c("asymMitsch", "tauxMitsch", "enviMitsch", "optDose", "optYield")
  mitschMC <- as.data.frame(mitschMC)
  
  for (i in 1:nbMC) {
    # coefficientss
    ## generate random
    coefficientsMC[i, ] <- rmvnorm(1, coefficients[, 1], vcov(mm))
    
    ## coefficientss of specific model matrix
    asymMitsch <- modmat %*% coefficientsMC[i, asymRows]
    tauxMitsch <- modmat %*% coefficientsMC[i, tauxRows]
    enviMitsch <- modmat %*% coefficientsMC[i, enviRows]
    doseMitsch <- newData$Dose
    
    
    # Create tags according to specified constrains
    if (is.null(constrain.asym)) {
      asymTag <- asymMitsch[1] + ranEf
    } else {
      asymTag <- constrain.asym[1] + constrain.asym[2]/(1+exp(-(asymMitsch[1] + ranEf)))
    }
    
    if (is.null(constrain.env)) {
      envTag <- enviMitsch[1]
    } else {
      envTag <- constrain.env[1] + constrain.env[2]/(1+exp(-(enviMitsch[1])))
    }
    
    if (is.null(constrain.rate)) {
      rateTag <- tauxMitsch[1]
    } else {
      rateTag <- constrain.rate[1] + constrain.rate[2]/(1+exp(-(tauxMitsch[1])))
    }
    
    # prediction rendement
    predMC_list[[i]] <- data.frame(Dose = doseMitsch,
                                   Response = mitschFunc(doseMitsch, 
                                                         asymTag,
                                                         rateTag, 
                                                         envTag,
                                                         rateExp = rateExp))
    
    # optimal dose
    doseOpt <- doseOptMitsherlich(asymptote = asymTag,
                                  taux = rateTag,
                                  environnement = envTag,
                                  prix_dose = C,
                                  prix_vendable = V,
                                  interets = 0,
                                  periode = 0,
                                  rateExp = rateExp)
    rdtOpt <- mitschFunc(doseOpt, 
                         asymTag,
                         rateTag, 
                         envTag,
                         rateExp = rateExp)
    
    mitschMC[i, ] <- c((asymMitsch[1] + ranEf), tauxMitsch[1], enviMitsch[1],
                       doseOpt, rdtOpt)
  }
  
  # Optimal
  asymMitsch <- modmat %*% coefficients[asymRows, 1]
  tauxMitsch <- modmat %*% coefficients[tauxRows, 1]
  enviMitsch <- modmat %*% coefficients[enviRows, 1]
  doseMitsch <- newData$Dose
  
  doseOpt <- doseOptMitsherlich(asymptote = asymMitsch[1] + ranEf,
                                taux = tauxMitsch[1],
                                environnement = enviMitsch[1],
                                prix_dose = C,
                                prix_vendable = V,
                                interets = 0,
                                periode = 0,
                                rateExp = rateExp)
  
  rdtOpt <- mitschFunc(doseOpt, 
                       asymMitsch[1] + ranEf,
                       tauxMitsch[1], 
                       enviMitsch[1],
                       rateExp = rateExp)
  
  names(predMC_list) <- 1:length(predMC_list)
  predMC <- ldply(predMC_list)
  predMC <- predMC[is.finite(predMC$Response), ]
  predMC <- predMC[predMC$Response > 1e-3, ]
  predMC$.id <- factor(predMC$.id)
  modelMC <- data.frame(newData$Dose, predTB$pred)
  
  if(!is.null(refTest)) {
    dataMC <- data.frame(essaiBase$Dose, essaiBase$Response, NoBloc = essaiBase$NoBloc)
  } 
  
  gg <- ggplot(predMC, aes(x = Dose, y = Response)) +
    geom_line(aes(colour = .id)) +
    geom_line(data = modelMC, aes(x = newData.Dose, y = predTB.pred), lwd = 1) +
    geom_point(data = mitschMC, aes(x = optDose, y = optYield), shape = 3, alpha = 0.1) +
    scale_color_manual(values = rep(rgb(0, 0, 0, alphaMC), nbMC)) +
    theme_bw() +
    theme(legend.position="none") 
  
  if(!is.null(refTest)) {
    gg <- gg + 
      geom_point(data = dataMC, 
                 aes(x = essaiBase.Dose, y = essaiBase.Response, shape = factor(NoBloc)),
                 size = 3)
  } 
  
  return(list(gg = gg, optDose = doseOpt, optYield = rdtOpt,
              optDoseMC = mitschMC$optDose, optYieldMC = mitschMC$optYield,
              pred = data.frame(Dose = newData$Dose, Response = predTB$pred),
              metaParam = predTB$metaParam, coefficientsMC = coefficientsMC,
              ranEf = ranEf, model.matrix = modmat))
}


paperPlotMC <- function(odm_o,
                        data,
                        ess,
                        image = "none",
                        xlab = expression("Dose P (kg⋅ha"^"-1"~")"), 
                        ylab = expression("Yield (Mg⋅ga"^"-1"~")"),
                        denslab = "Density",
                        xlim, ylim,
                        alphaAB = 0.5,
                        alphaMClines = 0.05,
                        alphaOptPts = 0.3,
                        alphaImage = 0.5,
                        text.size = 1,
                        text.axis = 1,
                        pch = 3, pch.cex = 1,
                        plotLines = TRUE,
                        plotOptPoints = TRUE,
                        kernel = "normal",
                        side.density = TRUE,
                        optLocation = "lines",
                        optShow = "both",
                        extra_points = NULL) {
  library(MASS)
  plot.new()
  if(side.density) {
    graphics::layout(matrix(c(2,0,1,3), ncol=2, byrow=TRUE),
                     widths = c(4/5, 1/5),
                     heights = c(1/5, 4/5))
  }
  par(mar = c(5,5,1,1), xaxt = "s", yaxt = "s")
  if(image == "image" | image == "contour") {
    if(kernel == "normal") {
      kd <- kde2d(odm_o$optDoseMC, odm_o$optYieldMC, n = 500)
      optx <- which.max(apply(kd$z, 1, max))
      opty <- which.max(apply(kd$z, 2, max))
    } else if(kernel == "lognormal") {
      cond <- odm_o$optDoseMC > 0 & odm_o$optYieldMC > 0
      if(sum(cond) <= 1) {
        kd <- NULL
      } else {
        kd <- kde2d(log10(odm_o$optDoseMC[cond]),
                    log10(odm_o$optYieldMC[cond]), 
                    n = 500)
        kd$x <- 10^kd$x
        kd$y <- 10^kd$y
        optx <- which.max(apply(kd$z, 1, max))
        opty <- which.max(apply(kd$z, 2, max))
      }
      
    } else {
      stop("Invalid kernel. Should be either 'normal' of 'lognormal'.")
    }
    
    if (image == "image" & !is.null(kd)) {
      image(kd,
            xlim = xlim, ylim = ylim,
            xlab = xlab,
            ylab = ylab,
            las = 1,
            col = rev(gray.colors(n = 10, alpha = alphaImage)),
            cex.lab = text.size, cex.axis = text.axis, cex.main = text.size, cex.sub = text.size)
      if(optShow == "both" | optShow == "kernel") {
        if (optLocation == "lines") {
          abline(v = kd$x[optx], col = rgb(0, 0, 0, alphaAB), lty = 3)
          abline(h = kd$y[opty], col = rgb(0, 0, 0, alphaAB), lty = 3)
        } else if (optLocation == "point") {
          points(x = kd$x[optx], y = kd$y[opty], bg = rgb(0.5, 0.5, 0.5, alphaAB), cex = 2, pch = 21)
        }
      }
    } else if (image == "contour" & !is.null(kd)) {
      contour(kd,
              xlim = xlim, ylim = ylim,
              xlab = xlab,
              ylab = ylab,
              las = 1,
              drawlabels = FALSE,
              col = rgb(0, 0, 0, alphaImage),
              cex.lab = text.size, cex.axis = text.axis, cex.main = text.size, cex.sub = text.size)
      if(optShow == "both" | optShow == "kernel") {
        if (optLocation == "lines") {
          abline(v = kd$x[optx], col = rgb(0, 0, 0, alphaAB), lty = 3)
          abline(h = kd$y[opty], col = rgb(0, 0, 0, alphaAB), lty = 3)
        } else if (optLocation == "point") {
          points(x = kd$x[optx], y = kd$y[opty], bg = rgb(0.5, 0.5, 0.5, alphaAB), cex = 2, pch = 21)
        }
      }
    }
    
  } else {
    kd <- NULL
    plot(x = odm_o$gg$data$Dose, y = odm_o$gg$data$Response,
         xlim = xlim, ylim = ylim,
         type = "n",
         xlab = xlab,
         ylab = ylab,
         las = 1,
         cex.lab = text.size, cex.axis = text.axis, cex.main = text.size, cex.sub = text.size)
    grid(lty = 1, col = "grey90")
  }
  
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = NULL)
  if (plotLines) {
    for(i in seq_along(levels(odm_o$gg$data$.id))) {
      dataTmp_pl <- odm_o$gg$data[odm_o$gg$data$.id == levels(odm_o$gg$data$.id)[i], ]
      lines(x = dataTmp_pl$Dose, y = dataTmp_pl$Response, col = rgb(0, 0, 0, alphaMClines))
    }
  }
  if (plotOptPoints) {
    points(x = odm_o$optDoseMC, y = odm_o$optYieldMC, pch = pch, cex = pch.cex,
           col = rgb(0, 0, 0, alphaOptPts))
  }
  if (!is.null(ess)) {
    df_refTest <- data %>% filter(NoEssai == ess) %>% droplevels()
    for (i in 1:length(levels(df_refTest$NoBloc))) {
      filtCond <- df_refTest$NoBloc == levels(df_refTest$NoBloc)[i]
      points(x = df_refTest$Dose[filtCond], 
             y = df_refTest$Response[filtCond],
             pch = 14+i, col = rgb(0.1, 0.1, 0.1, alphaAB))
    }
  }
  if (!is.null(extra_points)) {
    points(x = extra_points[, 1], y = extra_points[, 2],
           pch = 1)
  }
  lines(x = odm_o$pred$Dose, y = odm_o$pred$Response, lwd = 2)
  
  if(optShow == "both" | optShow == "center") {
    if (optLocation == "lines") {
      abline(h = odm_o$optYield, col = rgb(0, 0, 0, alphaAB), lty = 2)
      abline(v = odm_o$optDose, col = rgb(0, 0, 0, alphaAB), lty = 2)
    } else if (optLocation == "point") {
      points(x = odm_o$optDose, y = odm_o$optYield, col = 'white', cex = 2.0, pch = 16) # pch = 24, bg = rgb(0.5, 0.5, 0.5, alphaAB)
      points(x = odm_o$optDose, y = odm_o$optYield, col = 'black', cex = 1.8, pch = 16)
    }
  }
  
  if(side.density) {
    par(mar = c(0,5,1,1), xaxt = "n", yaxt = "s")
    plot(x = density(odm_o$optDoseMC)$x,
         y = density(odm_o$optDoseMC)$y,
         type = "l",
         xlim = xlim,
         bty = "n",
         ylab = denslab,
         las = 1,
         cex.lab = text.size, cex.axis = text.axis, cex.main = text.size, cex.sub = text.size)
    abline(v = odm_o$optDose, col = rgb(0, 0, 0, alphaAB), lty = 2)
    if(!is.null(kd)) abline(v = kd$x[optx], col = rgb(0, 0, 0, alphaAB), lty = 3)
    lines(x = density(odm_o$optDoseMC)$x,
          y = density(odm_o$optDoseMC)$y)
    par(mar = c(5,0,1,1), yaxt = "n", xaxt = "s")
    plot(x = density(odm_o$optYieldMC)$y,
         y = density(odm_o$optYieldMC)$x,
         type = "l",
         ylim = ylim,
         bty = "n",
         xlab = denslab,
         las = 1,
         cex.lab = text.size, cex.axis = text.axis, cex.main = text.size, cex.sub = text.size)
    abline(h = odm_o$optYield, col = rgb(0, 0, 0, alphaAB), lty = 2)
    if(!is.null(kd)) abline(h = kd$y[opty], col = rgb(0, 0, 0, alphaAB), lty = 3)
    lines(x = density(odm_o$optYieldMC)$y,
          y = density(odm_o$optYieldMC)$x)
  }
}


compareEffects <- function(data, refTests, mm, rhs, col_dose = 'Dose', logYield = FALSE,
                           rateExp = 'exp') {
  
  gg <- list()
  for (i in 1:length(refTests)) {
    essaiBase <- data %>% filter(NoEssai == refTests[i])
    newdata <- essaiBase[1, ]
    mm_mpred <- mitschPred(mm, newdata = newdata, rhs = rhs, col_dose = col_dose, 
                           logYield = logYield,
                           rateExp = rateExp)
    mm_mpred_gg <- data.frame(Parameter = rownames(mm_mpred$drag$A), Asymptote = mm_mpred$drag$A, 
                              Rate = mm_mpred$drag$R, Environment = mm_mpred$drag$E)
    rownames(mm_mpred_gg) <- 1:nrow(mm_mpred_gg)
    gg[[i]] <- melt(mm_mpred_gg)
  }
  names(gg) <- refTests
  gg <- ldply(gg)
  names(gg)[3] <- 'MetaParamCat'
  
  ggplot(filter(gg, Parameter != '(Intercept)'), aes(x = value, y = Parameter)) +
    theme_bw() +
    geom_point(aes(shape = .id), size = 3) +
    facet_grid(. ~ MetaParamCat, scale = 'free')
}


# Sensitivity analysis
cases_f <- function(data, nbCase, ref, quant.start = 2, quant.end = 4, fac.default = c("Superior", "LRC"),
                    custom.eval = NULL) {
  if (class(data[ref][[1]]) == 'factor' & nlevels(data[ref][[1]]) != nbCase) {
    nbCase <- nlevels(data[ref][[1]])
    warning("nbCase has been forced to equal nlevels of ref")
  }
  
  cases <- data.frame(mat.or.vec(nr = nbCase, nc = ncol(data)))
  colnames(cases) <- colnames(data)
  
  # fill factor columns
  for (i in 1:length(fac.default)) {
    cases[sapply(data, class) == 'factor'][, i] <- factor(rep(fac.default[i], times = nbCase),
                                                          levels = levels(data[sapply(data, class) == 'factor'][, i][[1]]))
  }
  
  # fill ref column
  if (class(data[ref][[1]]) == 'factor') {
    cases[ref] <- factor(levels(data[ref][[1]]), levels = levels(data[ref][[1]]))
  } else {
    if(is.null(custom.eval)) {
      cases[ref] <- seq(quantile(data[ref][[1]])[quant.start], quantile(data[ref][[1]])[quant.end], length = nbCase)
      #if (nbCase %% 2 != 0) { # if odd replace center by 0
      #  cases[nbCase/2 + 0.5, ref] <- 0
      #}
    } else {
      cases[ref] <- custom.eval
    }
  }
  return(cases)
}

sensitivity.mitsch.parameters <- function(cases, dose, mm, rhs,
                                          prix_dose = 1.20, prix_vendable = 275,
                                          rateExp = 'exp') {
  predL <- list()
  optL <- list()
  caseA <- c()
  caseR <- c()
  caseE <- c()
  for (i in 1:nrow(cases)) {
    newDat <- cases[rep(i, length(dose)), ]
    newDat$dose <- dose
    model_pred <- mitschPred(mm = mm, newdata = newDat, rhs = rhs, col_dose = "dose",
                             #constrain.asym = c(0, 60), 
                             logYield = FALSE, ranEf = 0, rateExp = rateExp)
    optDose <- doseOptMitsherlich(asymptote = model_pred$metaParam$A, #60/(1+exp(-model_pred$metaParam$A)), 
                                  taux = model_pred$metaParam$R,
                                  environnement = model_pred$metaParam$E, 
                                  prix_dose = prix_dose, prix_vendable = prix_vendable,
                                  rateExp = rateExp)
    optYield <- mitschFunc(dose = optDose,
                           asymptote = model_pred$metaParam$A,#  60/(1+exp(-model_pred$metaParam$A)), 
                           taux = model_pred$metaParam$R,
                           environnement = model_pred$metaParam$E,
                           rateExp = rateExp)
    optL[[i]] <- data.frame(Dose = optDose, Pred = optYield)
    predL[[i]] <- data.frame(Dose = dose, Pred = model_pred$pred)
    caseA <- c(caseA, model_pred$metaParam$A)
    caseR <- c(caseR, model_pred$metaParam$R)
    caseE <- c(caseE, model_pred$metaParam$E)
  }
  return(list(caseA = caseA, caseR = caseR, caseE = caseE, predL = predL, optL = optL))
}


sensitivity.mitsch.plot <- function(pred.cases, cases, class.case, tag,
                                    xlim = c(0, 320), ylim = c(25, 45)) {
  predL <- pred.cases$predL
  optL <- pred.cases$optL
  names(predL) <- rownames(cases)
  names(optL) <- rownames(cases)
  predL <- ldply(predL, data.frame)
  textPos_predL <- predL %>%
    group_by(.id) %>%
    summarise(maxY = max(Pred), maxD = max(Dose))
  optL <- ldply(optL, data.frame)
  
  ggplot(predL, aes(x = Dose, y = Pred)) +
    theme_bw() +
    theme(legend.position="none") +
    geom_line(aes(colour = .id)) +
    geom_label(data = textPos_predL, aes(x = maxD+3, y = maxY, label = .id),
               colour = "white", fill = 'grey15', hjust = 0) +
    #annotate('text', label = tag, 
    #         x = Inf, y = Inf, hjust = 1.05, vjust = 1.5) +
    geom_point(data = optL, aes(x = Dose, y = Pred), size = 2.5, colour = 'white') +
    geom_point(data = optL, aes(x = Dose, y = Pred), size = 2.0, colour = 'black') +
    xlab(expression("Dose N (kg ha"^"-1"~")")) +
    ylab(expression("Yield (Mg ha"^"-1"~")")) +
    xlim(xlim) +
    ylim(ylim) +
    scale_colour_manual(values = rep('black', times = length(unique(predL$.id))))
}

