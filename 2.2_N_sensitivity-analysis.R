########################
# Sensitivity analysis #
########################
# To be executed previously:
## source('0_init.R')
## source('1_climate.R')
## source('2.0_mmPreProcess.R')
## source('2.1_N_multilevel-model.R')

## All variables are constant except one, which is subjected
## to contrasts

## Optimal dose

### This is the same loop as ran previously
### added to avoid messing plot outputs
metaParamA <- c()
metaParamR <- c()
metaParamE <- c()

for (i in levels(df_mm$NoEssai)) { # find A, E and R for each field trial
  ess <- i
  mcN <- optDoseMitschMC(mm = mmN, rhs = rhs, col_dose = 'Dose',
                         logYield = FALSE,
                         refTest = ess, data = df_mm, nbNd = 50, doseMax = 300,
                         nbMC = 1, C = 1.20, V = 275, alphaMC = 0.02, ranEfLevel = 1,
                         rateExp = 'ordinary')
  metaParamA <- c(metaParamA, mcN$metaParam$A)
  metaParamR <- c(metaParamR, mcN$metaParam$R)
  metaParamE <- c(metaParamE, mcN$metaParam$E)
}
#### Find median of A, E and R
cteA <- median(metaParamA)
cteR <- median(metaParamR)
cteE <- median(metaParamE)

#### Sequence of n numbers from the min to the max of A, E and R
n <- 50
varA <- seq(min(metaParamA), max(metaParamA), length = n)
varR <- seq(min(metaParamR), max(metaParamR), length = n)
varE <- seq(min(metaParamE), max(metaParamE), length = n)

#### Compute the optimal dose for each sequence, where other parameters are constant
optA <- doseOptMitsherlich(asymptote = varA, #60/(1+exp(-(varA))),
                           taux = cteR,
                           environnement = cteE,
                           prix_dose = 1.2, prix_vendable = 275,
                           rateExp = 'ordinary')
optR <- doseOptMitsherlich(asymptote = cteA, #60/(1+exp(-(cteA))),
                           taux = varR,
                           environnement = cteE,
                           prix_dose = 1.2, prix_vendable = 275,
                           rateExp = 'ordinary')
optE <- doseOptMitsherlich(asymptote = cteA, #60/(1+exp(-(cteA))),
                           taux = cteR,
                           environnement = varE,
                           prix_dose = 1.2, prix_vendable = 275,
                           rateExp = 'ordinary')

ggdf <- data.frame(Type = rep(c('Asymptote', 'Environment', 'Rate'), each = n),
                   Parameter = c(varA, #60/(1+exp(-(varA))),
                                 varE,
                                 exp(varR)),
                   Dose = c(optA, optE, optR))
ggdf$Typelabels <- factor(ggdf$Type, labels = c(expression("Asymptote (Mg ha"^"-1"~")"),
                                                expression("Environment (kg ha"^"-1"~")"),
                                                'Rate'))

ggplot(ggdf, aes(x = Parameter, y = Dose)) +
  geom_line() +
  facet_grid(. ~ Typelabels, scales = 'free', labeller = label_parsed) + #) +
  ylab(expression("Optimal economic dose N (kg ha"^"-1"~")")) +
  ylim(c(0,300)) +
  theme_bw()
ggsave(filename = 'images/N_optDose.png', scale = 1.2, dpi = 300,
       height = 3, width = 6)
#ggsave(filename = 'images/N_optDose.tiff', scale = 1.2, dpi = 300, height = 3, width = 7)

## A grid of values to find the spread of the slope
gridParA <- expand.grid(E = varE, R = varR)
gridParE <- expand.grid(A = varA, R = varR)
intercA <- c()
slopeA <- c()
intercE <- c()
slopeE <- c()
for (i in 1:nrow(gridParA)) {
  # Asmptote
  optA <- doseOptMitsherlich(asymptote = varA,
                             taux = gridParA[i, 2],
                             environnement = gridParA[i, 1],
                             prix_dose = 1.5, prix_vendable = 200,
                             rateExp = 'ordinary')
  lmA <- lm(optA ~ varA)
  intercA[i] <- coef(lmA)[1]
  slopeA[i] <- coef(lmA)[2]

  # ENvironment
  optE <- doseOptMitsherlich(asymptote = gridParE[i, 1],
                             taux = gridParE[i, 2],
                             environnement = varE,
                             prix_dose = 1.5, prix_vendable = 200,
                             rateExp = 'ordinary')
  lmE <- lm(optE ~ varE)
  intercE[i] <- coef(lmE)[1]
  slopeE[i] <- coef(lmE)[2]
}
plot(density(slopeA))
plot(density(slopeE))
plot(density(intercA))
plot(density(intercE))

mean(slopeA)

slopedf <- data.frame(Asymptote = slopeA, Environment = slopeE)
intercdf <- data.frame(Asymptote = intercA, Environment = intercE)
gglmpar <- melt(bind_rows('Slope' = slopedf, 'Intercept' = intercdf, .id = 'id'))
boxplot(data = filter(gglmpar, id == 'Slope'), value ~ variable, main = 'Slope')
boxplot(data = filter(gglmpar, id == 'Intercept'), value ~ variable, main = 'Intercept')

names(gglmpar)
class(gglmpar$value)
ggplot(gglmpar, aes(x = value)) +
  facet_grid(. ~ id, scales = 'free') +
  geom_density(aes(fill = variable))



# Sensibility of variables in the model
## Climate: 9 cases
### build scaler matrices for mean and sd to put predictors to scale
cli_scaler.mean <- matrix(ncol = 3, nrow = 2)
colnames(cli_scaler.mean) <- c('ppt', 'temp', 'sdi')
rownames(cli_scaler.mean) <- c('Beginning of the season to 45 DAP', '46 DAP to harvest')
cli_scaler.sd <- cli_scaler.mean
cli_scaler.mean[1, 1] <- mean(df_mm_filt.log['pcpCum_0_45'][[1]])
cli_scaler.mean[2, 1] <- mean(df_mm_filt.log['pcpCum_46_end'][[1]])
cli_scaler.mean[1, 2] <- mean(df_mm_filt.log['DJCb_0_45'][[1]])
cli_scaler.mean[2, 2] <- mean(df_mm_filt.log['DJCb_46_end'][[1]])
cli_scaler.mean[1, 3] <- mean(df_mm_filt.log['SDI_0_45'][[1]])
cli_scaler.mean[2, 3] <- mean(df_mm_filt.log['SDI_46_end'][[1]])
cli_scaler.sd[1, 1] <- sd(df_mm_filt.log['pcpCum_0_45'][[1]])
cli_scaler.sd[2, 1] <- sd(df_mm_filt.log['pcpCum_46_end'][[1]])
cli_scaler.sd[1, 2] <- sd(df_mm_filt.log['DJCb_0_45'][[1]])
cli_scaler.sd[2, 2] <- sd(df_mm_filt.log['DJCb_46_end'][[1]])
cli_scaler.sd[1, 3] <- sd(df_mm_filt.log['SDI_0_45'][[1]])
cli_scaler.sd[2, 3] <- sd(df_mm_filt.log['SDI_46_end'][[1]])

### create predictors from cliClassMean (created in 1.0_climate.R)
cliClassMean_log <- cliClassMean
#cliClassMean_log[, 3:5] <- log10(cliClassMean[, 3:5])

cliClassMean_log.sc <- cliClassMean_log
for(i in 1:6) {
  for(j in 3:5) {
    mean_here <- cli_scaler.mean[rownames(cli_scaler.mean) == cliClassMean_log.sc$period[i], colnames(cliClassMean_log.sc)[j]]
    sd_here <- cli_scaler.sd[rownames(cli_scaler.sd) == cliClassMean_log.sc$period[i], colnames(cliClassMean_log.sc)[j]]
    cliClassMean_log.sc[i, j] <- (cliClassMean_log[i, j] - mean_here)/sd_here
  }
}

### dose and cases
doseN <- seq(0, 250, 5)
nbCase <- 9

### build cases (finally!)
cases <- data.frame(#growing.season = rep(0, nbCase),
                    #Cultivar = rep("Superior", nbCase),
                    Maturity.Order = rep(0, nbCase),
                    PrecCrop.Five.Classes = rep("LRC", nbCase),
                    soilTextIlr1 = rep(0, nbCase),
                    soilTextIlr2 = rep(0, nbCase),
                    soilTextIlr3 = rep(0, nbCase),
                    soilTypeIlr1_3 = rep(0, nbCase),
                    soilTypeIlr2_3 = rep(0, nbCase),
                    soil_pH = rep(0, nbCase),
                    ##### les cliClassMean_log.sc sont mal calculés ##############
                    # 0 to 45, dry-dry-dry-opt-opt-opt-wet-wet-wet
                    pcpCum_0_45 = unlist(c(cliClassMean_log.sc[1, 3], cliClassMean_log.sc[1, 3], cliClassMean_log.sc[1, 3],
                                           cliClassMean_log.sc[2, 3], cliClassMean_log.sc[2, 3], cliClassMean_log.sc[2, 3],
                                           cliClassMean_log.sc[3, 3], cliClassMean_log.sc[3, 3], cliClassMean_log.sc[3, 3])),
                    DJCb_0_45 = unlist(c(cliClassMean_log.sc[1, 4], cliClassMean_log.sc[1, 4], cliClassMean_log.sc[1, 4],
                                           cliClassMean_log.sc[2, 4], cliClassMean_log.sc[2, 4], cliClassMean_log.sc[2, 4],
                                           cliClassMean_log.sc[3, 4], cliClassMean_log.sc[3, 4], cliClassMean_log.sc[3, 4])),
                    SDI_0_45 = unlist(c(cliClassMean_log.sc[1, 5], cliClassMean_log.sc[1, 5], cliClassMean_log.sc[1, 5],
                                         cliClassMean_log.sc[2, 5], cliClassMean_log.sc[2, 5], cliClassMean_log.sc[2, 5],
                                         cliClassMean_log.sc[3, 5], cliClassMean_log.sc[3, 5], cliClassMean_log.sc[3, 5])),
                    # 46-end, dry-opt-wet-dry-opt-wet-dry-opt-wet
                    pcpCum_46_end = unlist(c(cliClassMean_log.sc[4, 3], cliClassMean_log.sc[5, 3], cliClassMean_log.sc[6, 3],
                                           cliClassMean_log.sc[4, 3], cliClassMean_log.sc[5, 3], cliClassMean_log.sc[6, 3],
                                           cliClassMean_log.sc[4, 3], cliClassMean_log.sc[5, 3], cliClassMean_log.sc[6, 3])),
                    DJCb_46_end = unlist(c(cliClassMean_log.sc[4, 4], cliClassMean_log.sc[5, 4], cliClassMean_log.sc[6, 4],
                                           cliClassMean_log.sc[4, 4], cliClassMean_log.sc[5, 4], cliClassMean_log.sc[6, 4],
                                           cliClassMean_log.sc[4, 4], cliClassMean_log.sc[5, 4], cliClassMean_log.sc[6, 4])),
                    SDI_46_end = unlist(c(cliClassMean_log.sc[4, 5], cliClassMean_log.sc[5, 5], cliClassMean_log.sc[6, 5],
                                           cliClassMean_log.sc[4, 5], cliClassMean_log.sc[5, 5], cliClassMean_log.sc[6, 5],
                                           cliClassMean_log.sc[4, 5], cliClassMean_log.sc[5, 5], cliClassMean_log.sc[6, 5]))
)

cases$PrecCrop.Five.Classes <- factor(cases$PrecCrop.Five.Classes,
                                      levels = levels(df_mm$PrecCrop.Five.Classes))
#cases$Cultivar <- factor(cases$Cultivar,
#                                      levels = levels(df_mm$Cultivar))

rownames(cases) <- c("Dry-Dry", "Dry-Moderate", "Dry-Wet",
                     "Moderate-Dry", "Moderate-Moderate", "Moderate-Wet",
                     "Wet-Dry", "Wet-Moderate", "Wet-Wet")

predL <- list()
optL <- list()
caseA <- c()
caseR <- c()
caseE <- c()
for (i in 1:nrow(cases)) {
  newDat <- data.frame(cases[rep(i, length(doseN)), ], dose = doseN)
  model_pred <- mitschPred(mm = mmN, newdata = newDat, rhs = rhs, col_dose = "dose",
                           #constrain.asym = c(0, 60),
                           logYield = FALSE, ranEf = 0, rateExp = 'ordinary')
  optDose <- doseOptMitsherlich(asymptote = model_pred$metaParam$A, #60/(1+exp(-model_pred$metaParam$A)),
                                taux = model_pred$metaParam$R,
                                environnement = model_pred$metaParam$E,
                                prix_dose = 1.20, prix_vendable = 275, rateExp = 'ordinary')
  optYield <- mitschFunc(dose = optDose,
                         asymptote = model_pred$metaParam$A,#  60/(1+exp(-model_pred$metaParam$A)),
                         taux = model_pred$metaParam$R,
                         environnement = model_pred$metaParam$E, rateExp = 'ordinary')
  optL[[i]] <- data.frame(Dose = optDose, Pred = optYield)
  predL[[i]] <- data.frame(Dose = doseN, Pred = model_pred$pred)
  caseA <- c(caseA, model_pred$metaParam$A)
  caseR <- c(caseR, model_pred$metaParam$R)
  caseE <- c(caseE, model_pred$metaParam$E)
}

summary(caseA)
summary(exp(caseR))
summary(caseE)

names(predL) <- rownames(cases)
names(optL) <- rownames(cases)
predL <- ldply(predL, data.frame)
textPos_predL <- predL %>%
  group_by(.id) %>%
  summarise(maxY = max(Pred), maxD = max(Dose))

optL <- ldply(optL, data.frame)

ggplot9climates <- function(predL, highlighted) {
  gg_palette <- rep("gray70", 9) ; gg_palette[highlighted] <- 'black'
  gg_linewidth <- rep(0.5, 9) ; gg_linewidth[highlighted] <- 1.5
  gg_ptshape <- rep(1, 9) ; gg_ptshape[highlighted] <- 16

  ggplot(predL, aes(x = Dose, y = Pred)) +
    theme_bw() +
    geom_line(aes(colour = .id,size = .id)) +
    #geom_text(data = textPos_predL, aes(x = maxD, y = maxY, label = .id),
    #          hjust=0, vjust=0.25, size = 4) +
    geom_point(data = optL, aes(x = Dose, y = Pred), size = 4.5, colour = 'white') +
    geom_point(data = optL, aes(x = Dose, y = Pred, shape = .id, colour = .id), size = 4) +
    geom_text(data = textPos_predL[highlighted, ], aes(x = maxD, y = maxY, label = .id),
              hjust=0, vjust=0.25, size = 4) +
    xlab(expression("Dose N (kg ha"^"-1"~")")) +
    ylab(expression("Yield (Mg ha"^"-1"~")")) +
    xlim(c(0, 320)) +
    ylim(c(25, 45)) +
    scale_colour_manual(values = gg_palette) +
    scale_size_manual(values = gg_linewidth) +
    scale_shape_manual(values = gg_ptshape) +
    theme(legend.position="none")
}

#svg("images/N_9climatsD.svg", width = 5, height = 5, pointsize = 12)
png("images/N_9climatsD.png",width=1500, height=1500, res = 300) #
ggplot9climates(predL, highlighted = 1:3)
dev.off()

#svg("images/N_9climatsM.svg", width = 5, height = 5, pointsize = 12)
png("images/N_9climatsM.png",width=1500, height=1500, res = 300) #
ggplot9climates(predL, highlighted = 4:6)
dev.off()

#svg("images/N_9climatsW.svg", width = 5, height = 5, pointsize = 12)
png("images/N_9climatsW.png",width=1500, height=1500, res = 300) #
ggplot9climates(predL, highlighted = 7:9)
dev.off()

## 3 cases
sensibility.analysis.tags <- read.csv('csv/sensibility-analysis-tags.csv', sep = ';',
                                      header = TRUE, row.names = 1, stringsAsFactors = FALSE)
### change tags for low, median, high for numeric variables
sensibility.analysis.tags[with(sensibility.analysis.tags, grep(x = Low, pattern = 'Low')), 1] <- 'Low'
sensibility.analysis.tags[with(sensibility.analysis.tags, grep(x = Medium, pattern = 'Median')), 2] <- 'Median'
sensibility.analysis.tags[with(sensibility.analysis.tags, grep(x = High, pattern = 'High')), 3] <- 'High'


### all except maturity order
gg_sensitivity <- list()
rangecases_l <- list()
rangecases <- matrix(ncol = 4, nrow = 2)
colnames(rangecases) <- c('A', 'E', 'R', 'DoseOpt')
for (i in var_mm[-1]) {
  cases <- cases_f(data = select(df_mm, one_of(var_mm)), nbCase = 3, ref = i, quant.start = 2, quant.end = 4,
                   fac.default = c("Small grain")) # 'Superior',
  if (class(df_mm[i][[1]]) == 'numeric') {
    rownames(cases) <- as.matrix(sensibility.analysis.tags)[rownames(sensibility.analysis.tags) == i, ]
  } else {
    rownames(cases) <- levels(df_mm[i][[1]])
  }
  pred.cases <- sensitivity.mitsch.parameters(cases = cases, dose = seq(0, 250, 1), mm = mmN,
                                              rhs = rhs, prix_dose = 1.20, prix_vendable = 275,
                                              rateExp = 'ordinary')
  rangecases[, 1] <- range(pred.cases$caseA)
  rangecases[, 2] <- range(pred.cases$caseE)
  rangecases[, 3] <- range(pred.cases$caseR)
  rangecases[, 4] <- range(unlist(pred.cases$optL)[c(1, 3, 5)])
  rangecases_l[[i]] <- rangecases
  gg_sensitivity[[i]] <- sensitivity.mitsch.plot(pred.cases = pred.cases, cases = cases,
                                                 xlim = c(0, 320), ylim = c(20, 45),
                                                 class.case = class(df_mm[i][[1]]),
                                                 tag = translate_col[translate_col$to_name_mm == i, 3])
  #ggsave(paste0('images/sensitivity/', i, '.png')) # .svg
}

### maturity order
cases <- cases_f(data = select(df_mm, one_of(var_mm)), nbCase = 5, ref = var_mm[1],
                 custom.eval = sort(unique(df_mm$Maturity.Order)),
                 fac.default = c("Small grain")) # 'Superior',
rownames(cases) <- 1:5#paste('Maturity order', 1:5)
pred.cases <- sensitivity.mitsch.parameters(cases = cases, dose = seq(0, 250, 1), mm = mmN,
                                            rhs = rhs, prix_dose = 1.20, prix_vendable = 275,
                                            rateExp = 'ordinary')
range(pred.cases$caseA)
range(pred.cases$caseE)
range(pred.cases$caseR)
pred.cases$optL
gg_sensitivity$Maturity.Order <- sensitivity.mitsch.plot(pred.cases = pred.cases, cases = cases,
                                                         xlim = c(0, 320), ylim = c(20, 45),
                                                         tag = translate_col[translate_col$to_name_mm == 'Maturity.Order', 3])
#ggsave(paste0('images/sensitivity/', var_mm[1], '.png'))

### grid of plots
remx <- var_mm[-c(6,  11)]
remy <- var_mm[-c(2:6)]
for (i in remx) {
  gg_sensitivity[[i]] <- gg_sensitivity[[i]] +
    xlab('') +
    theme(axis.text.x = element_blank())
}
for (i in remy) {
  gg_sensitivity[[i]] <- gg_sensitivity[[i]] +
    ylab('') +
    theme(axis.text.y = element_blank())
}


svg(filename = 'images/all2.svg', width = 10, height = 14, pointsize = 12)
multiplot(plotlist = gg_sensitivity, cols = 3) # export RStudio 900 X 1200
dev.off()

i = 10
gg_sensitivity[[i]]
names(gg_sensitivity)[i]

## Sold price
### does not work for now...
cases <- cases_f(data = select(df_mm, one_of(var_mm)), nbCase = 1, ref = var_mm[1],
                 fac.default = c("Small grain")) # 'Superior',
cases[1,1] <- sort(unique(df_mm$Maturity.Order))[3]
cases$PrecCrop.Five.Classes <- factor(cases$PrecCrop.Five.Classes,
                                      levels = levels(df_mm$PrecCrop.Five.Classes))
doseN <- seq(0, 250, 1)
rownames(cases) <- 1

priceCases <- 100
sold_price <- seq(250, 300, length = priceCases)
dose_price <- seq(1, 1.5, length = priceCases)
prices <- expand.grid(sold_price, dose_price)

predL <- list()
optL <- list()
newDat <- data.frame(cases[rep(1, length(doseN)), ], dose = doseN)
model_pred <- mitschPred(mm = mmN, newdata = newDat, rhs = rhs, col_dose = "dose",
                         logYield = FALSE, ranEf = 0, rateExp = 'ordinary')
optDose <- c()
for (i in 1:nrow(prices)) {
  optDose[i] <- doseOptMitsherlich(asymptote = model_pred$metaParam$A,
                                   taux = model_pred$metaParam$R,
                                   environnement = model_pred$metaParam$E,
                                   prix_dose = prices[i, 2], prix_vendable = prices[i, 1],
                                   rateExp = 'ordinary')
}

ggdf <- data.frame(prices, optDose)
names(ggdf) <- c("Potato.price", "Dose.N.price", "Optimal.Dose")

library(directlabels)
ggPrices <- ggplot(ggdf, aes (x = Potato.price, y = Dose.N.price, z = Optimal.Dose)) +
  theme_bw() +
  geom_tile(aes(fill = Optimal.Dose)) +
  stat_contour(aes(colour=..level..)) +
  scale_colour_gradient(high = 'black', low = 'black') +
  scale_fill_gradient(high = 'black', low = 'white') +
  theme(legend.position="none")
direct.label(ggPrices, method = "first.points") # "maxvar.points"
