#####################################
## Climate indexes: axis reduction ##
#####################################
# To be executed previously:
## source('0_init.R')

# Selected columns for climate variables according to their importance and their redundancy
clicol_.45 <- c("pcpCum_0_45", "SDI_0_45", "etpKcCum_0_45", "DJCb_0_45", "totEv3po3j_0_45") # , "nb_ppt10_0_45"
clicol_46. <- c("pcpCum_46_end", "SDI_46_end", "etpKcCum_46_end", "DJCb_46_end", "totEv3po3j_46_end") # , "nb_ppt10_46_end"

# The objective is to create a single procedure to reduce climate variables.
# Here, I do compute weights for two periods (.45 and 46.) and take their
# medians to compute climate scores

indcli_.45 <- data_N[, clicol_.45]
colnames(indcli_.45) <-  c("pcpCum", "SDI", "etpKcCum", "DJCb", "totEv3po3j") # , "nb_ppt10"
indcli_46. <- data_N[, clicol_46.]
colnames(indcli_46.) <- colnames(indcli_.45)

# remove duplicated lines (climate indexes are computed on a trial basis)
indcli_.45 <- unique(indcli_.45)
indcli_46. <- unique(indcli_46.)

# melt both periods and standardize
indcli_ <- rbind(indcli_.45, indcli_46.)
indcli_$period <- factor(c(rep('.45', times = nrow(indcli_.45)),
                           rep('46.', times = nrow(indcli_46.))))
st_ <-  data.frame(decostand(indcli_[, -6], method = "standardize", na.rm = TRUE),
                   period = indcli_$period)

apply(st_[, -6], 2, function(x) mean(x, na.rm = TRUE))
apply(st_[, -6], 2, function(x) sd(x, na.rm = TRUE))

# before axis reduction, remove NAs and outliers
indcli_NA <- indcli_[!apply(indcli_, 1, anyNA), ]
indcli_01 <- sign1(indcli_NA[, -6], qcrit=0.975)$wfinal01 # check: sum(st_01 == 0)/length(st_01)
indcli_ar <- indcli_NA[indcli_01 == 1, ]

# standarize
st_ar <-  data.frame(scale(indcli_NA[indcli_01 == 1, -6]),
                     period = indcli_NA$period[indcli_01 == 1])
apply(st_ar[, -6], 2, function(x) mean(x, na.rm = TRUE))
apply(st_ar[, -6], 2, function(x) sd(x, na.rm = TRUE))

# correlations between indexes
st_ar.cor <- corr.test(x = st_ar[, -6],
                       use = "pairwise",method="pearson",
                       adjust="holm", alpha = .05)
st_ar.cor$r
st_ar.cor$p
st_ar.cor.m <- data.frame(melt(st_ar.cor$r), p.value = melt(st_ar.cor$p)[, 3])
colnames(st_ar.cor.m)[3] <- 'correlation'
st_ar.cor.m$test <- ifelse(st_ar.cor.m$p.value <= 0.05, 'significant', 'not significant')

for (i in 1:nlevels(st_ar.cor.m$Var1)) {
  levels(st_ar.cor.m$Var1)[i] <-  translate_col$to_name_mm_long[
    translate_col$from_name_mm == levels(st_ar.cor.m$Var1)[i]
    ]
  levels(st_ar.cor.m$Var2)[i] <-  translate_col$to_name_mm_long[
    translate_col$from_name_mm == levels(st_ar.cor.m$Var2)[i]
    ]
}

dict <- c()
st_ar.cor.m$is.diag <- NA
i=1
for (i in 1:nrow(st_ar.cor.m)) {
  curr <- paste(st_ar.cor.m[i, 1], st_ar.cor.m[i, 2])
  if(curr %in% dict) {
    st_ar.cor.m$is.diag[i] <- TRUE
  } else {
    st_ar.cor.m$is.diag[i] <- FALSE
  }
  dict[i] <- paste(st_ar.cor.m[i, 2], st_ar.cor.m[i, 1])
}

png('images/corcli.png', height = 3000, width = 4000, res = 600)
ggplot(st_ar.cor.m[!st_ar.cor.m$is.diag, ], aes(x = Var1, y = Var2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  xlab('') +
  ylab('') +
  geom_tile(aes(fill = correlation)) +
  geom_text(aes(label = round(correlation, 2))) +
  scale_fill_gradient2(name = "Correlation",
                       low = "#ece2f0",  high = "#1c9099") + # PuBuGn, http://colorbrewer2.org/
  coord_fixed()
dev.off()

# Axis reduction
pca_st <- rda(st_ar[, -6], scale = FALSE)
scaling <- 2
pca_st_s <- scores(pca_st, display = "sites", scaling = scaling)
pca_st_l <- scores(pca_st, display = "species", scaling = scaling)

st_all_.45 <- data.frame(matrix(ncol = length(clicol_.45), nrow = nrow(data_N)))
names(st_all_.45) <- names(st_ar[, -6])
st_all_46. <- data.frame(matrix(ncol = length(clicol_46.), nrow = nrow(data_N)))
names(st_all_46.) <- names(st_ar[, -6])
for(i in 1:ncol(st_all_.45)) {
  st_all_.45[, i] <- (data_N[, clicol_.45[i]][[1]] - mean(indcli_[, i][[1]], na.rm = T)) / sd(indcli_[, i][[1]], na.rm = T)
  st_all_46.[, i] <- (data_N[, clicol_46.[i]][[1]] - mean(indcli_[, i][[1]], na.rm = T)) / sd(indcli_[, i][[1]], na.rm = T)
}

apply(rbind(st_all_.45, st_all_46.), 2, mean) # close to 0 (unequal number of observations per site)
apply(rbind(st_all_.45, st_all_46.), 2, sd) # close to 1 (unequal number of observations per site)

scores_.45 <- as.matrix(st_all_.45) %*% pca_st_l
scores_46. <- as.matrix(st_all_46.) %*% pca_st_l

head(predict(pca_st, newdata = st_all_.45, scaling = scaling, type = "wa"))

scores_.45 <- predict(pca_st, newdata = st_all_.45, scaling = scaling, type = "wa")[, 1:2]
scores_46. <- predict(pca_st, newdata = st_all_46., scaling = scaling, type = "wa")[, 1:2]

apply(rbind(scores_.45, scores_46.), 2, quantile) # unequal number of observations per site

## ggplot
period <- factor(c(rep(".45", nrow(data_N)), rep("46.", nrow(data_N))))
classe_cli <- as.factor(c(as.character(data_N$ClasseClimatique_.45),
                          as.character(data_N$ClasseClimatique_46.)))
ggloadings <- 0.2 * data.frame(pca_st_l)#4 * data.frame(PC1=ra_weights[,1], PC2=ra_weights[,2])
rownames(ggloadings) <- c('PPT', 'SDI', 'ETC', 'GDD', 'NLE')
ggdf_cli <- unique(data.frame(period,
                              rbind(scores_.45, scores_46.),
                              pcpCum = c(data_N$pcpCum_0_45, data_N$pcpCum_46_end),
                              DJC = c(data_N$DJCb_0_45, data_N$DJCb_46_end),
                              SDI = c(data_N$SDI_0_45, data_N$SDI_46_end),
                              classe_cli))

levels(ggdf_cli$period) <- c('Beginning of the season to 45 DAP',
                             '46 DAP to harvest')

TextFrame <- data.frame(PC1 = ggloadings$PC1, PC2 = ggloadings$PC2,
                        LAB = rownames(ggloadings))
TextFrame <- transform(TextFrame,
                       w = strwidth(LAB, 'inches') + 0.05,
                       h = strheight(LAB, 'inches') + 0.05)

climate_gg <- ggplot(data=ggdf_cli, aes(x=PC1, y=PC2)) +
  theme_bw() +
  geom_point(aes(fill=classe_cli, shape=period), size=2.5, colour = 'gray30', alpha = 0.75) +
  geom_segment(data=ggloadings,
               aes(x=0, y=0, xend=PC1, yend=PC2), size=1.5, colour = 'white',
               arrow=arrow(length=unit(0.3,"cm"))) +
  geom_segment(data=ggloadings,
               aes(x=0, y=0, xend=PC1, yend=PC2), size=0.5,
               arrow=arrow(length=unit(0.2,"cm"))) +
  geom_label(data = ggloadings, aes(PC1, PC2, label = rownames(ggloadings)),
             colour = "white", fill = 'grey15', hjust = 0, nudge_x = 0.005,
             nudge_y = c(0.02, -0.075, 0.08, -0.07, -0.07)) +
  xlab("PC1") +
  ylab("PC2") +
  scale_fill_manual(name = "Thornthwaite class",
                    values = c("#d8b365", "#f5f5f5", "#5ab4ac"),
                    guide=guide_legend(override.aes = aes(shape = 21))) +
  scale_shape_manual(name = "Period",
                     values = c(21, 24))
climate_gg

#svg("climateBiplot.svg", width = 7.5, height = 6, pointsize = 12)
png("images/climateBiplot.png", width = 4000, height = 3500 * 0.7, res = 600)
climate_gg
dev.off()

# Thornthwait classes means
cliClassMean <- ggdf_cli %>%
  group_by(period, classe_cli) %>%
  #summarise(mean_ppt = mean(PC1), meanPC2 = mean(PC2))
  summarise(ppt = mean(pcpCum), temp = mean(DJC), sdi = mean(SDI))

# Correlations between indices
pairs(data_N[, clicol_.45])
pairs(data_N[, clicol_46.])

# insert scores into the main table
data_N$pc1_.45 <- scores_.45[, 1]
data_N$pc2_.45 <- scores_.45[, 2]
data_N$pc1_46. <- scores_46.[, 1]
data_N$pc2_46. <- scores_46.[, 2]
