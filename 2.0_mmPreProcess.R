##########################################
## PREPARATION OF THE MIXED MODEL TABLE ##
##########################################
# To be executed previously:
## source('0_init.R')
## source('1_climate.R')

# Select columns to include
col_Nmm <- c(# Response
               "RendVendable",
               "RendVendableMaxParEssai",

               # Random effect
               "NoEssai", "NoBlocF", "NoTraitement", #  "NoBlocC"

               # Divers
               #'SourceProjet',

               # Common effects
               ## Cultivar
               "Cultivar", "Maturity.Order", "growing.season",

               ## Rotation
               "PrecCrop.Five.Classes",

               ## Soil
               "ClasseTexture",

               ## Water
               "Irrigation", "GroupeHydrologique",

               ## climate
               "pc1_.45", "pc2_.45", "pc1_46.", "pc2_46.",
               "pcpCum_0_45", "pcpCum_46_end", "DJCb_0_45", "DJCb_46_end", "SDI_0_45", "SDI_46_end",

               # Doses
               ## Nitrogen
               "NtotDose",

               # balances soil type
               "soilTypeIlr1.3", "soilTypeIlr2.3",

               # balances soil texture
               "soilTextIlr1", "soilTextIlr2", "soilTextIlr3",

               # soil
               "Csol_imp", "pH1dans1_H2O_imp",

               # Extra climate variables
               "pcpCum_0_45", "DJCb_0_45", "pcpCum_46_end", "DJCb_46_end"
)
data.frame(col_Nmm, col_Nmm %in% colnames(data_N))

## Filter to exclude tests where the minimum yield is too low
data_mm <- data_N %>%
  select(one_of(col_Nmm))

colnames(data_mm)[colnames(data_mm) == 'NoBlocF'] <- 'NoBloc' # Or NoBlocC

## Cultivar column, if exists, should be a factor
if("Cultivar" %in% names(data_mm)) {
  data_mm$Cultivar <- factor(data_mm$Cultivar)
  data_mm$Cultivar <- relevel(data_mm$Cultivar, ref = "Superior")
}

## Filter to remove blocks containing only one value (some exists with NoBlocC)
NoEssai_NoBloc <- factor(paste0(data_mm$NoEssai, '_', data_mm$NoBloc))
isUniqueInBlock <- levels(NoEssai_NoBloc)[table(NoEssai_NoBloc) == 1]
data_mm <- data_mm[!NoEssai_NoBloc %in% isUniqueInBlock, ]

## Replace 0 doses by an arbitrary low value to be able to test log of dose instead of dose
data_mm$RendVendable[data_mm$RendVendable == 0] <- 0.0001

## Remove dots from column names (dots are used by nlme to separate meta-parameters from parameters)
### note single elements (e.g. soil_Al) must be listed after balances to avoid replacing it befor
### parsing the balance (if soil_Al is parsed before, soil_Al.P is translated to soil Al.P)
for (i in 1:nrow(translate_col)) {
  index <- which(colnames(data_mm) == translate_col$from_name_mm[i])
  colnames(data_mm)[index] <- translate_col$to_name_mm[i]
}

data_mm$NoBloc <- factor(data_mm$NoBloc)

## Hydrologic group is well represented by soilTypeIlr1_3, soilTypeIlr1_5 and soilTypeIlr3_5
boxplot(soilTypeIlr1_3 ~ GroupeHydrologique, data_mm, xlab = "Hydrologic group", ylab = "Soil type ilr1_3")
boxplot(soilTypeIlr2_3 ~ GroupeHydrologique, data_mm, xlab = "Hydrologic group", ylab = "Soil type ilr2_3")

## Max dose per trial (to filter out later trials tested with low doses)
N_max.dose.trial <- data_mm %>%
  group_by(NoEssai) %>%
  summarize(N_maxDoseTrial = max(NtotDose))
N_min.dose.trial <- data_mm %>%
  group_by(NoEssai) %>%
  summarize(N_minDoseTrial = min(NtotDose))
data_mm <- left_join(data_mm, N_max.dose.trial, 'NoEssai')
data_mm <- left_join(data_mm, N_min.dose.trial, 'NoEssai')

# Number of dose levels per bloc
nDose <- data_mm %>%
  select(one_of(c("NoEssai", "NoBloc", "NtotDose"))) %>%
  group_by(NoEssai, NoBloc) %>%
  summarise(N_numDoseLevels = length(unique(NtotDose)))
data_mm <- left_join(data_mm, nDose, by = c('NoEssai', 'NoBloc'))
