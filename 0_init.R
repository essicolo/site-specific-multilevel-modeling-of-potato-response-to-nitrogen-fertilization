###################
# Load librairies #
###################
library(MASS) # version 7.3-47 ; kde2d: for kernel density estimation
library(plyr) # version 1.8.4 ; generic package
library(dplyr)# version 0.7.4 ; generic package
library(reshape2) # version 1.4.2 ; melt
library(ggplot2) # version 2.2.1 ; generic package
library(grid) # version 3.4.2 ; arrows in ggplots
library(compositions) # version 1.40-1 ; acomp, ilr, ilrInv: for CoDa
library(robCompositions) # version 2.0.6 ; impKNNa: for data imputation
library(mvoutlier) # version 2.0.8 ; sign1: for outlier detection
library(psych) # version 1.7.8 ; principal: for rotated PCA
library(vegan) # version 2.4-4 ; rda
library(nlme) # version 3.1-131 ; nlme: for multilevel mixed-effect modeling
library(e1071) # version 1.6-8 ; kurtosis: for Kurtosis test
library(caret) # version 6.0-77 ; machine learning
library(maptools) # version 0.9-2 ; map
library(mapproj) # version 1.2-5 ; map
gpclibPermit()

# Correct MASS and dplyr possible interference
filter <- dplyr::filter
select <- dplyr::select


###############
# Import data #
###############
load("potato_nitrogen_data_set.RData")


#######################
# Map of observations #
#######################
canmap@data$IDUU <- rownames(canmap@data)
canmap.points <- fortify(canmap, region="IDUU")
names(canmap@data)[3] <- 'id'
canmap.df <- join(canmap.points, canmap@data, by="id")

quemap.df <- canmap.df %>%
  filter(NAME == "Quebec") %>%
  droplevels()

insetrect <- data.frame(xmin = -75, xmax = -60,
                        ymin = 44.5, ymax = 50)

textSize <- 14
canada.plot <- ggplot(canmap.df) +
  geom_polygon(aes(x = long, y = lat, group=group), fill = 'grey70') +
  geom_polygon(aes(x = long, y = lat, group=group), data = quemap.df, fill = 'grey25') +
  geom_path(aes(x = long, y = lat, group=group), color = "white") +
  geom_rect(data = insetrect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            alpha = 0, colour = "black", size = 1, linetype = 1) +
  theme(panel.background = element_rect(fill = 'white')) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  #theme(axis.ticks.length = unit(0.001, "mm")) + labs(x=NULL, y=NULL) +
  theme(plot.margin = unit(c(-1, 0, -1, -1), "cm")) +
  labs(x = '', y = '') +
  coord_equal()

quebec.plot <- ggplot(canmap.df, aes(long, lat, group = group)) +
  geom_polygon(fill = 'grey75') +
  geom_polygon(data = quemap.df, fill = 'grey25') +
  geom_path(color = "white") +
  geom_point(data = data_N, aes(x = LonDD, y = LatDD, group = 1),
             size = 4, pch = 21, fill = 'white') +
  theme(panel.background = element_rect(fill = 'white')) +
  theme(axis.text.x = element_text(colour = 'black', size = textSize),
        axis.text.y = element_text(colour = 'black', size = textSize),
        axis.title.x = element_text(colour = 'black', size = textSize),
        axis.title.y = element_text(colour = 'black', size = textSize),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x = 'Longitude', y = 'Latitude') +
  coord_fixed() +
  coord_map(xlim = c(insetrect$xmin, insetrect$xmax),
            ylim = c(insetrect$ymin, insetrect$ymax))

png("images/map.png",width=3000, height=1600, res = 300) #
grid.newpage()
vp_quebec <- viewport(width = 1, height = 1, x = 0.5, y = 0.5)  # the larger map
vp_canada <- viewport(width = 0.4, height = 0.4, x = 0.743, y = 0.753)  # the inset in upper right
print(quebec.plot, vp = vp_quebec)
print(canada.plot, vp = vp_canada)
dev.off()


##################
# Arrange tables #
##################
# Preceding crop classes (5 groups)
data_N <- left_join(x = data_N, y = amalgamation_precCrop, by = 'PrecedentCultural') # join
data_N$PrecCrop.Five.Classes <- relevel(data_N$PrecCrop.Five.Classes, ref = "Small grain") # LRC as first level

# Cultivar classes, by maturity order
data_N$Cultivar = relevel(factor(data_N$Cultivar), ref = "Superior") # Correct cultivar factors. Factor "SUperior" should come first.
amalgamation_Cultivar <- left_join(amalgamation_Cultivar, amalgamation_Cultivar_Order5, "Maturity.5")
data_N <- left_join(x = data_N,
                     y = amalgamation_Cultivar[c('Cultivar', 'Maturity.Order')],
                     by = 'Cultivar') # join

# Texture categories
## amalgamate categories
amalgamation_Texture <- amalgamation_Texture %>%
  filter(Texture %in% levels(data_N$Texture)) %>%
  droplevels() # filter to take only the levels present in data_N$Texture

## push to data_N
data_N <- left_join(x = data_N, y = amalgamation_Texture, by = 'Texture')

## Numeric values instead of classes: use centroid of the class
### Table including centroids
qcTextCentroid <- qcTextCentroid %>%
  filter(Texture %in% levels(data_N$Texture)) %>%
  droplevels() # filter to take only the levels present in data_N$Texture

### Left join based on the Texture column
data_N <- left_join(data_N, qcTextCentroid, by = "Texture")

### if SableMesure is NA, replace by the centroid of the textural category. same for ArgileMesure
data_N <- data_N %>%
  mutate(SableMesure = ifelse(is.na(SableMesure), SableCentroid, SableMesure)) %>%
  mutate(ArgileMesure = ifelse(is.na(ArgileMesure), ArgileCentroid, ArgileMesure))
data_N <- data_N %>% mutate(LimonDiff = 100 - SableMesure - ArgileMesure) # compute silt from clay and sand
### Note: compositional transformation of textures are placed after soil elements imputation

# Climate categories (Thornthwait class)
ClasseClim_func <- function(ETP, PPT, irrigation) {
  classeCli <- rep(NA, length(ETP))
  budget <- (PPT - ETP) / ETP
  classeCli[budget < -0.1] <- "dry"
  classeCli[budget > 0.167] <- "wet"
  classeCli[budget >= -0.1 & budget <= 0.167] <- "moderate"
  classeCli[irrigation] <- "moderate"
  return(factor(classeCli))
}

## CLimate category 0-45 days
data_N <- data_N %>%
  mutate(ClasseClimatique_.45 = ClasseClim_func(ETP = etpKcCum_0_45, PPT = pcpCum_0_45, irrigation = Irrigation))

## Climate category 46-end of season
data_N <- data_N %>%
  mutate(ClasseClimatique_46. = ClasseClim_func(ETP = etpKcCum_46_end, PPT = pcpCum_46_end, irrigation = Irrigation))


#########################
# N-Fertilization tests #
#########################
tmpNsemisDose <- data_N$NsemisDose
tmpNfracDose <- data_N$NfracDose
tmpNsemisDose[is.na(data_N$NsemisDose)] <- 0
tmpNfracDose[is.na(data_N$NfracDose)] <- 0
data_N$NtotDose <- tmpNsemisDose + tmpNfracDose


###########
# Soil Al #
###########
tmpAlSolTabi <- data_N$AlSolTabi
tmpAlSol <- data_N$AlSol
tmpAlSolTabi[is.na(data_N$AlSolTabi)] <- 0
tmpAlSol[is.na(data_N$AlSol)] <- 0
data_N$AlSolTot <- tmpAlSolTabi + tmpAlSol
data_N$AlSolTot[data_N$AlSolTot == 0] <- NA


#################
# Soil elements #
#################
# Data imputation of soil analysis results
soil_col <- c("CaSol", "MgSol", "Ksol", "FeSol", "AlSolTot", "Psol") ### Exclude "Bsol",  "CuSol", "ZnSol", and "MnSol" because too few data available

## Portrait of missing data # matrixplot(data_N[, soil_col]) # data_N[2514:2785, soil_col]
## no soil analysis data from Origine == "Inventaire provincial" and
## SourceProject == "Giroux, M. 1982. Can. J Soil Sc., 62: 503-517.".
## The imputation will be made by the mean of full data after imputation

data_N$soilAnalysis_allNA <- apply(data_N[soil_col], 1, function(X) all(is.na(X))) # keep track of rows where soil analysis is empty
data_N$soilAnalysis_anyNA <- apply(data_N[soil_col], 1, function(X) any(is.na(X))) # keep track of rows where soil analysis contains any NA
data_N$soilAnalysis_countNA <- apply(data_N[soil_col], 1, function(X) sum(is.na(X))) # number of NAs in soil data
data_N$soilAnalysis_hasTooMuchNA <- data_N$soilAnalysis_countNA > 3 # there is 6 components. not imputing more than half

# Impute with kNN (robCompositions) and adding the organic matter (long process)
set.seed(23165)
soil_imputeElemRob <- impCoda(as.matrix(data_N[!data_N$soilAnalysis_hasTooMuchNA, soil_col]),
                               maxit = 10, eps = 0.5, method = "ltsReg",
                               closed = TRUE, init = "KNN", k = 10, noise = 0)
colnames(soil_imputeElemRob$xImp) <- paste0(soil_col, '_imp')

# push imputed soil elements into the main table
data_N <- left_join(x = data_N,
                     y = data.frame(NoEssai = data_N$NoEssai[!data_N$soilAnalysis_hasTooMuchNA],
                                    NoBloc = data_N$NoBloc[!data_N$soilAnalysis_hasTooMuchNA],
                                    NoTraitement = data_N$NoTraitement[!data_N$soilAnalysis_hasTooMuchNA],
                                    soil_imputeElemRob$xImp),
                     by = keys_col) # join


###########
# Soil pH #
###########
# correct pH with method CaCl2 to H2O
data_N$pH1dans1_H2O <- data_N$pH1dans1
data_N$pH1dans1_H2O[data_N$pHactifMethode == 'CaCl2'] <- 0.27 + 1.03 * data_N$pH1dans1[data_N$pHactifMethode == 'CaCl2']

# Impute pH from soil analysis, but apart because the Aitchison metric is inappropriate
data_N$pHOM_allNA <- apply(data_N[c("pH1dans1_H2O", "Csol")], 1, function(X) all(is.na(X))) # keep track of rows where soil analysis is empty
data_N$pHOM_anyNA <- apply(data_N[c("pH1dans1_H2O", "Csol")], 1, function(X) any(is.na(X))) # keep track of rows where soil analysis contains any NA
data_N$pHOM_countNA <- apply(data_N[c("pH1dans1_H2O", "Csol")], 1, function(X) sum(is.na(X))) # number of NAs in soil data

data_impOMpH <- data_N %>%
  select(one_of(c(keys_col, soil_col, "pH1dans1_H2O", "Csol"))) %>%
  filter(!data_N$soilAnalysis_hasTooMuchNA, !data_N$pHOM_allNA)

data_impOMpH[c(soil_col, 'Csol')] <- log10(data_impOMpH[c(soil_col, 'Csol')]) # impute on log

soil_imputeOMpHRob <- impKNNa(as.matrix(select(data_impOMpH, -one_of(keys_col))), metric = "Euclidean",
                              k = 6, primitive = TRUE,
                              normknn = TRUE, adj = 'median') # impute pH and C
soil_imputeOMpHRob$xImp[, 8] <- 10^soil_imputeOMpHRob$xImp[, 8] # C back to original scale

soilImputedOMpH <- tbl_df(data.frame(data_impOMpH[keys_col],
                                     soil_imputeOMpHRob$xImp[, c(7, 8)]))
colnames(soilImputedOMpH)[c(4, 5)] <- c("pH1dans1_H2O_imp", "Csol_imp")

# push to data_N
data_N <- left_join(x = data_N,
                     y = soilImputedOMpH,
                     by = keys_col) # join

# Soil textural balances: transform texture proportions to balances
# [C | Mineral], [Sand,Silt | Clay] and [Sand | Silt]
soilTextComp <- acomp(data_N[c('Csol_imp', 'SableMesure', 'LimonDiff', 'ArgileMesure')])
soilTextSBP <- matrix(c(-1, 1, 1, 1,
                        0, 1, 1,-1,
                        0, 1,-1, 0),
                      byrow = TRUE,
                      ncol = 4)
soilTextBal <- ilr(soilTextComp, V = gsi.buildilrBase(t(soilTextSBP)))
soilTextBal <- ilrNA(soilTextComp, soilTextSBP, soilTextBal)

### store balances in the main data frame
data_N$soilTextIlr1 <- soilTextBal[, 1]
data_N$soilTextIlr2 <- soilTextBal[, 2]
data_N$soilTextIlr3 <- soilTextBal[, 3]


#########
# Yield #
#########
rdt_max <- rep(NA, nrow(data_N)) # max yield per trial (empty vector).
for (i in 1:nlevels(data_N$NoEssai)) {
  filtre <- data_N$NoEssai == levels(data_N$NoEssai)[i] # filter for test i
  rdt_max[filtre] <- max(data_N$RendVendable[filtre], na.rm=TRUE)
}
rdt_max[!is.finite(rdt_max)] <- NA
data_N$RendVendableMaxParEssai <- rdt_max


###################
# Soil morphology #
###################
# Pedological soil type
soilTypeSBP.3 <- matrix(c(-1, 1,-1,
                          -1, 0, 1),
                        byrow = TRUE,
                        ncol = 3)

colnames(soilTypeSBP.3) <- c('loam.gley', 'sand.podzol', 'sand.gley')
soilTypeComp.3 <- acomp(data_N[c('m.3a', 'm.3b', 'm.3c')])
soilTypeBal.3 <- ilr(soilTypeComp.3, V = gsi.buildilrBase(t(soilTypeSBP.3)))
soilTypeBal.3 <- ilrNA(soilTypeComp.3, soilTypeSBP.3, soilTypeBal.3)

data_N$soilTypeIlr1.3 <- soilTypeBal.3[, 1] # ilr1 [gley | podzol]
data_N$soilTypeIlr2.3 <- soilTypeBal.3[, 2] # ilr2 gley.[fine | coarse]

soilTypeBal.3.NA <- apply(soilTypeBal.3, 1, anyNA)
soilTypeBal.3.pca <- princomp(soilTypeBal.3[!soilTypeBal.3.NA, ])
biplot(soilTypeBal.3.pca)


########
# MISC #
########
# Rename blocs where two or more sources are added
data_N$NoBlocF <- factor(data_N$NoBloc)
data_N$NoBlocC <- as.character(data_N$NoBloc)

# Length of growing season
data_N$growing.season <- data_N$DOYrecolte - data_N$DOYplantation
