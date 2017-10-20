##################
## Mixed model  ##
##################
# To be executed previously:
## source('0_init.R')
## source('1_climate.R')
## source('2.0_mmPreProcess.R')

# Variables to include in the model
var_mm <- c(#"growing.season",
  "Maturity.Order",
  #"ionomicGroup",
  #"Cultivar",
  "PrecCrop.Five.Classes",
  "soilTextIlr1", "soilTextIlr2", "soilTextIlr3",
  "soilTypeIlr1_3", "soilTypeIlr2_3",
  "soil_pH",
  #"pc1Temp_x45", "pc2Prec_x45", "pc1Temp_46x", "pc2Prec_46x",
  "pcpCum_0_45", "pcpCum_46_end",
  "DJCb_0_45", "DJCb_46_end",
  "SDI_0_45", "SDI_46_end"
)

# define table for mm analysis
df_mm_filt <- data_mm %>%
  filter(RendVendableMaxParEssai >= 28) %>%
  filter(!(NoEssai %in% c(32, 40:46, 123:129))) %>% # N tests mixed with P and K tests
  #filter(!(NoEssai == "39" & NoBloc == "1")) %>% # mixed P source
  #filter(!(NoEssai == "49" & NoBloc == "2")) %>% # mixed P source
  filter(!(NoEssai == "50" & NoBloc == "2")) %>% # probable field errors
  filter(N_numDoseLevels >= 4) %>%
  filter(N_maxDoseTrial >= 100) %>%
  filter(N_minDoseTrial <= 50) %>%
  # "RendVendable" "PoidsSpec" "RendMoy"
  select(one_of(c("RendVendable", "NtotDose", "NoEssai", "NoBloc", var_mm))) %>%
  na.omit() %>% # many NAs in pH and OM, and some in yield
  droplevels() # Make sure not having extra factors to have a start vector of correct length

# log transform coefficients: no log transformation
##log_transf_cols <- c("pcpCum_0_45", "pcpCum_46_end",
##                     "DJCb_0_45", "DJCb_46_end",
##                     "SDI_0_45", "SDI_46_end")#, "growing.season")
df_mm_filt.log <- df_mm_filt
##df_mm_filt.log[, log_transf_cols] <- log10(df_mm_filt[, log_transf_cols])

df_mm_filt.log$Maturity.Order -> a
# dummy code Maturity order
if("Maturity.Order" %in% var_mm) {
  df_mm_filt.log$Maturity.Order <- with(df_mm_filt.log, model.matrix(~ordered(factor(Maturity.Order)))[, 2])
}

# scale variables
df_mm_filt.log.sc <- df_mm_filt.log
varnumcols <- sapply(df_mm_filt.log.sc, class) == "numeric" &
  colnames(df_mm_filt.log.sc) %in% var_mm
df_mm_filt.log.sc[, varnumcols] <- apply(df_mm_filt.log.sc[, varnumcols], 2, scale)

# renaming columns
names(df_mm_filt.log.sc)[1] <- 'Response'
names(df_mm_filt.log.sc)[2] <- 'Dose'

# data table for analysis
df_mm <- df_mm_filt.log.sc

# number of sites, blocks and observations used
length(unique(df_mm$NoEssai))
length(unique(paste(df_mm$NoEssai, df_mm$NoBloc)))
nrow(df_mm)


# Preparing the multilevel model
## starting values
start_list <- list()
for (i in 1:length(var_mm)) {
  if (is.factor(df_mm[var_mm[i]][[1]])) {
    start_list[[i]] <- rep(0, length(levels(df_mm[var_mm[i]][[1]]))-1)
  } else {
    start_list[[i]] <- 0
  }
}
start_vector <- unlist(start_list)

## Right-hand side of the formula to construct the formulas in the nlme model
rhs <- paste(var_mm, collapse= "+")


# Multilevel models
## Mitscherlich ordinary model
### this approach was slected.
### because the model is sensible to starting values
### we need an adequate and independent way to select them
### we loop through a grid and find the lowest AIC
# grid <- as.matrix(expand.grid(A = seq(35, 42, by = 2),
#                               R = seq(0.01, 0.09, by = 0.02),
#                               E = seq(0, 200, by = 50)))
# aic_n <- c()
# for (i in 1:nrow(grid)) {
#   mitsch_nlme_No <- mitscherlich_nlme(
#     startVal = grid[i, ], startVector = start_vector, # c(40, 0.05, 50)
#     rhs = rhs, data = df_mm, scaleResp = 'ordinary', rateExp = 'ordinary',
#     maxIter = 50, minScale = 1e-8, msVerbose = FALSE)
#   aic_n[i] <- AIC(mitsch_nlme_No)
#   print(paste(i, '/', nrow(grid)))
# }
# grid[which.min(aic_n), ]
# plot(aic_n)

mitsch_nlme_No <- mitscherlich_nlme(
  startVal = c(39, 0.09, 100),
  startVector = start_vector,
  rhs = rhs, data = df_mm, scaleResp = 'ordinary', rateExp = 'ordinary',
  maxIter = 50, minScale = 1e-8, msVerbose = FALSE)
AIC(mitsch_nlme_No)
residu(mitsch_nlme_No, level = 0, main = "Réponse: rendement")

mean(ranef(mitsch_nlme_No)[[1]][, 1])
sd(ranef(mitsch_nlme_No)[[1]][, 1])
shapiro.test(ranef(mitsch_nlme_No)[[1]][, 1])

mean(ranef(mitsch_nlme_No)[[2]][, 1])
sd(ranef(mitsch_nlme_No)[[2]][, 1])
shapiro.test(ranef(mitsch_nlme_No)[[2]][, 1])

range(summary(mitsch_nlme_No)$tTable[2:18, 1])

## Mitscherlich log10 model
mitsch_nlme_Nl <- mitscherlich_nlme(
  startVal = c(log10(36), 0.03, 25), startVector = start_vector,
  rhs = rhs, data = df_mm, scaleResp = 'log10', rateExp = 'ordinary',
  maxIter = 50, minScale = 1e-8, msVerbose = FALSE)
AIC(mitsch_nlme_Nl)
residu(mitsch_nlme_Nl)

## Square-root ordinary model
colw_nlme_No <- colwell_nlme(
  startVal = c(30, 0.5, -0.01), startVector = start_vector,
  rhs = rhs, data = df_mm, scaleResp = 'ordinary',
  maxIter = 50, minScale = 1e-8, msVerbose = FALSE)
AIC(colw_nlme_No)
residu(colw_nlme_No)
rsq(y = df_mm$Response, y_hat = predict(colw_nlme_No, level = 0))
rsq(y = df_mm$Response, y_hat = predict(colw_nlme_No, level = 1))
rsq(y = df_mm$Response, y_hat = predict(colw_nlme_No, level = 2))

## Square-root log10 model
colw_nlme_Nl <- colwell_nlme(
  startVal = c(30, 0.5, -0.01), startVector = start_vector,
  rhs = rhs, data = df_mm, scaleResp = 'log10',
  maxIter = 50, minScale = 1e-8, msVerbose = FALSE)
AIC(colw_nlme_Nl)
residu(colw_nlme_Nl)

## Quadratic ordinary model
quad_nlme_No <- quad_nlme(
  startVal = c(30, 0.5, -0.01), startVector = start_vector,
  rhs = rhs, data = df_mm, scaleResp = 'ordinary',
  maxIter = 50, minScale = 1e-8, msVerbose = FALSE)
AIC(quad_nlme_No)
residu(quad_nlme_No)

## Quadratic log10 model
quad_nlme_Nl <- quad_nlme(
  startVal = c(30, 0.5, -0.01), startVector = start_vector,
  rhs = rhs, data = df_mm, scaleResp = 'log10',
  maxIter = 50, minScale = 1e-8, msVerbose = FALSE)
AIC(quad_nlme_Nl)
residu(quad_nlme_Nl)

## Linear to max model (piecewise function)
### try a grid of starting values to find convergence
# grid <- as.matrix(expand.grid(intercept = seq(20, 26, by = 1),
#                               slope = seq(0.1, 0.2, by = 0.05),
#                               max_yield = seq(35, 45, by = 2.5)))
# aic_n <- c()
# for (i in 1:nrow(grid)) {
#   try(lp_nlme_No <- lp_nlme(
#     startVal = grid[i, ], startVector = start_vector,
#     rhs = rhs, data = df_mm, scaleResp = 'ordinary',
#     maxIter = 5, minScale = 1e-8, msVerbose = TRUE))
#   aic_n[i] <- AIC(lp_nlme_No)
#   print(paste(i, '/', nrow(grid)))
# }
# grid[which.min(aic_n), ]
# plot(aic_n)

lp_nlme_No <- lp_nlme(
  startVal = c(20, 0.2, 35), startVector = start_vector,
  rhs = rhs, data = df_mm, scaleResp = 'ordinary',
  maxIter = 10, minScale = 1e-8, msVerbose = TRUE)
### no convergence

### Comments
### The Mitscherlich model is more sound compared to Colwell.
### Compared to its log10 counterpart, The ordinary Mitscherlich model
### returns better residues. The ordinary Mitscherlich model is selected.

## mmN: The selected N mixed model
mmN <- mitsch_nlme_No
summary(mmN)$tTable
rsq(y = df_mm$Response, y_hat = predict(mmN, level = 0))
rsq(y = df_mm$Response, y_hat = predict(mmN, level = 1))
rsq(y = df_mm$Response, y_hat = predict(mmN, level = 2))
residu(residues = ranef(mitsch_nlme_No)[[1]][, 1], computeResidues = FALSE)
mean(ranef(mitsch_nlme_No)[[1]][, 1])
sd(ranef(mitsch_nlme_No)[[1]][, 1])
sd(ranef(mitsch_nlme_No)[[1]][, 1]) / summary(mmN)$tTable[1, 1]
residu(residues = ranef(mitsch_nlme_No)[[2]][, 1], computeResidues = FALSE)
mean(ranef(mitsch_nlme_No)[[2]][, 1])
sd(ranef(mitsch_nlme_No)[[2]][, 1])

isAsym <- grepl('Asymptote', rownames(summary(mmN)$tTable))
isRate <- grepl('Rate', rownames(summary(mmN)$tTable))
isEnv <- grepl('Environment', rownames(summary(mmN)$tTable))

are_coef <- data.frame(Asymptote = summary(mmN)$tTable[isAsym, 1][-1],
                       Rate = rate_coef <- summary(mmN)$tTable[isRate, 1][-1],
                       Environement = summary(mmN)$tTable[isEnv, 1][-1])
pairs(are_coef)


## Export plot for publication
p_mmN <- nlmePlotConfTable(mm = mmN,
                           conf_level = 0.95,
                           pval.breaks = c(0, 0.05, 0.1, 1.00), # p-value breaks
                           varCat = c('PrecCrop.Five.Classes', 'Intercept'), # variables which are categorical
                           varCatNamesTo = c('Prec. crop', 'Intercept'),
                           varCatOrder = c(1, 2, 3),
                           paramCat = c('Asymptote', 'Rate', 'Environment'), # name of parameters in the model
                           varNamesFrom = c(translate_col$to_name_mm),
                           varNamesTo = c(translate_col$to_name_mm_long),
                           greyEnd = 0.65,
                           removeIntercept = TRUE,
                           numApart = F,
                           varImpType = "cat",
                           limits_df = NULL,
                           modify.param.names = NULL)
p_mmN[[1]]#[[1]]
#p_mmN[[1]][[2]]
p_mmN[[2]] # http://www.researchgate.net/publication/263037296_A_review_of_statistical_methods_for_determination_of_relative_importance_of_correlated_predictors_and_identification_of_drivers_of_consumer_liking

png("images/N_fixefEf.png",width=3000, height=1400, res = 300) #
p_mmN[[1]]
dev.off()


## Export variance covariance matrix of the mixed model
mm.vcov <- vcov(mmN)

### Use understandable column and row names
for (i in 1:nrow(translate_col)) {
  new <- gsub(pattern = translate_col$to_name_mm[i],
              replacement = translate_col$to_name_mm_long[i],
              x = colnames(mm.vcov))
  colnames(mm.vcov) <- new
  rownames(mm.vcov) <- new
}

### Export
#write.table(mm.vcov, 'mmN.vcov.csv', sep = ';', dec = '.')

sup_mat2 <- data.frame(intervals(mmN)[[1]], p.val = summary(mmN)$tTable[, 5])
write.table(sup_mat2, 'mmN.intervals_supp_mat_2.csv', sep = ';', dec = '.')

# Optimal dose distribution
# -------------------------
## How: Generate random parameters from a multinormal distribution
## with means of parameters and their covariance, then predict the
## output from these randomly generated parameters. Each output in form of
## (yield vs dose) is then subjected to an economical analysis to obtain an
## economically optimal dose.

## Create a dataset with conditions where optimal dose is foreseek
## Easy to do by taking a particular test. Then interesting to compare
## to model with actual measurements.

# R² at level 1 for all tests
rsqMM <- c()
for (i in 1:length(levels(df_mm$NoEssai))) {
  ess <- levels(df_mm$NoEssai)[i]
  rsqMM[i] <- rsq(y = df_mm$Response[df_mm$NoEssai == ess],
                  y_hat = predict(mmN, level = 1)[df_mm$NoEssai == ess])
}
data.frame(levels(df_mm$NoEssai), rsqMM)

metaParam_L <- list()
#trials <- levels(df_mm$NoEssai)
trials <- c('4', '71', '89', '342')
metaparam_trials <- matrix(ncol = 4, nrow = length(trials))
colnames(metaparam_trials) <- c('Asymptote', 'Environment', 'Rate', 'Optimal dose')
rownames(metaparam_trials) <- trials
for (i in 1:length(trials)) {
  set.seed(555)
  ess <- trials[i]
  mc <- optDoseMitschMC(mm = mmN, rhs = rhs, col_dose = 'Dose',
                        logYield = FALSE,
                        refTest = ess, data = df_mm, nbNd = 50, doseMax = 300,
                        nbMC = 1, C = 1.20, V = 275, alphaMC = 0.02, ranEfLevel = 1,
                        rateExp = 'ordinary')

  metaparam_trials[i, 1] <- mc$metaParam$A
  metaparam_trials[i, 2] <- mc$metaParam$E
  metaparam_trials[i, 3] <- mc$metaParam$R
  metaparam_trials[i, 4] <- mc$optDose

  mc$optDose

  png(filename = paste0("images/N_", ess, "MC.png"),
      width = 1000, height = 1000, res = 150)
  paperPlotMC(odm_o = mc, data = df_mm,
              ess = ess,
              image = "none",
              alphaImage = 0.6,
              alphaMClines = 0.02, # 0.05 for ~ 500 MC, 0.02 for ~ 1000 MC
              alphaOptPts = 0.2, # 0.3 for ~ 500 MC, 0.2 for ~ 1000 MC
              xlab = expression("Dose N (kg ha"^"-1"~")"),
              ylab = expression("Yield (Mg ga"^"-1"~")"),
              denslab = 'Density',
              xlim = c(-5, 300), ylim = c(10, 50),
              pch = 3, pch.cex = 0.5,
              text.size = 1.5,
              text.axis = 1.2,
              plotLines = FALSE,
              plotOptPoints = FALSE,
              kernel = "normal",
              side.density = FALSE,
              optLocation = "point",
              optShow = "center")
  dev.off()
  print(ess)
}

data.frame(levels(df_mm$NoEssai), rsqMM, metaparam_trials)

summary(metaparam_trials)
metaparam_trials[rownames(metaparam_trials) %in% trials, ]

gg_metaparam_trials <- melt(metaparam_trials)
names(gg_metaparam_trials) <- c('Trial', 'Parameter', 'Value')

gg_metaparam_trials %>% filter(Trial == '400')
compareEffects(data = df_mm, refTests = '400', mm = mmN,
               rhs, col_dose = 'Dose', logYield = FALSE,
               rateExp = 'ordinary')
gg_metaparam_trials %>%
  group_by(Parameter) %>%
  summarise(q0 = quantile(Value, probs = 0),
            q25 = quantile(Value, probs = 0.25),
            q50 = quantile(Value, probs = 0.5),
            q75 = quantile(Value, probs = 0.75),
            q100 = quantile(Value, probs = 1))

quantile(gg_metaparam_trials$Value)
png("images/param_hist.png",width=4000, height=1200, res = 300) #
ggplot(gg_metaparam_trials, aes(Value)) +
  theme_bw() +
  facet_grid(. ~ Parameter, scales = 'free') +
  geom_histogram() +
  ylab('Count')
dev.off()

library(GGally)
png("images/corparameters.png",width=4000, height=4000, res = 600) #
ggpairs(data.frame(metaparam_trials),
        diag = list(continuous = "barDiag")) +
  theme_bw()
dev.off()

# Boxplots
var_mm_num <- var_mm[-2]
ggBP <- unique(melt(df_mm[var_mm_num])) # df_mm # data_mm

for (i in 1:nrow(translate_col)) {
  ggBP$variable <- sub(translate_col$to_name_mm[i], translate_col$to_name_mm_long[i], ggBP$variable)
}

selectedTests <- c('4', '71', '89', '342')

data_df %>%
  filter(NoEssai %in% selectedTests) %>%
  select(one_of(c('NoEssai', 'NoBloc', 'Cultivar',
                  'PrecCrop.Five.Classes', 'Maturity.Order', 'Texture'))) %>%
  unique() %>%
  print(n = Inf)

ggBP_test <- df_mm %>% # data_mm # df_mm
  filter(NoEssai %in% selectedTests) %>%
  group_by(NoEssai) %>%
  select(one_of(var_mm_num)) %>%
  summarise_each(funs(mean))

ggBP_test <- melt(ggBP_test, id = 'NoEssai')
for (i in 1:nrow(translate_col)) {
  ggBP_test$variable <- sub(translate_col$to_name_mm[i],
                            translate_col$to_name_mm_long[i],
                            ggBP_test$variable)
}

N_BP <- ggplot(ggBP,
               aes(x = variable, y = value)) +
  theme_bw() +
  geom_boxplot(outlier.colour = 'grey80', colour = 'grey50') +
  #geom_jitter() +
  geom_point(data = ggBP_test,
             aes(shape = NoEssai), size = 4, fill = rgb(0, 0, 0, 0.2)) +
  scale_shape_manual(name="Test ID", values = c(21, 22, 23, 24, 25)) +
  coord_flip() +
  xlab('') +
  ylab("Scaled variables")
N_BP


png("images/N_BP.png",width=3000, height=1200, res = 300) #
N_BP
dev.off()
