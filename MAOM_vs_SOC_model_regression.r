###E. Lugato

#measured SOM fractions on 352 LUCAS data('LCS.csv') were already made available at:
#https://esdac.jrc.ec.europa.eu/content/soil-organic-matter-som-fractions

#Here, we reported the comparison of linear vs polynomial fitting of MAOM vs total SOC.

#for further information please contact: emanuele.lugato@ec.europa.eu

############################################################################### 
library(ggplot2)


###read LUCAS and fraction data
LCS<-read.csv('LCS.csv', header=T)

###LUCAS selection in cropland, grassland and forest 
LCS<-subset(LCS, LU %in% c( "GR", "FR", "CR"))


fig<- ggplot(LCS, aes(OC_tf, OC_sc_g_kg )) +  geom_point(size = 1, alpha = 0.8, color="darkgreen") +
  theme_bw() + 
  facet_wrap(~LU) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = TRUE) +
  geom_smooth(method = "lm", formula = y ~ x+0, color="green") +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1.0) + 
  xlab("SOC g/kg soil") + ylab("MAOM C g/kg soil")


landuse<-"CR"  #select ("CR", "GR", "FR")

###model comparison with with intercept   
lm1 = lm(OC_sc_g_kg ~ OC_tf , data=subset(LCS, LCS$LU==landuse)) 
lm2 = lm(OC_sc_g_kg ~   poly(OC_tf,2, raw=T) , data=subset(LCS, LCS$LU==landuse))
anova(lm1,lm2)


###model comparison without intercept   
lm1 = lm(OC_sc_g_kg ~ OC_tf + 0, data=subset(LCS, LCS$LU==landuse)) 
lm2 = lm(OC_sc_g_kg ~   -1 + poly(OC_tf,2, raw=T) , data=subset(LCS, LCS$LU==landuse))
anova(lm1,lm2)


