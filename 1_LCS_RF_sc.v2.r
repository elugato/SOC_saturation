rm(list=ls(all=TRUE))

cm_prof<-20 ###set depth LUCAS profile

source("E:\\SIT\\caprese\\soil\\LUCAS\\LUCAS-2009-data-as-distributed-by-ESDAC\\lucasDB.R")



#############################################################################################
##addtional parmeters
#############################################################################################
omad<-stack("E:\\SIT\\caprese\\ven\\sripts_L\\exp\\_LINUXSERVER_soil\\VCFR\\SPT_lyr\\omadFAO.tif")
 irr<-stack("E:\\SIT\\caprese\\ven\\sripts_L\\exp\\_LINUXSERVER_soil\\VCFR\\SPT_lyr\\irr_FAO.tif")

NDEP<-raster("E:\\SIT\\fra\\db\\deposition_2006_to_2010\\ndep_grid")

WCLIM<-stack("E:\\SIT\\caprese\\meteo\\worldclim2\\clim_70-00LAEA.tif")


 RNPPin_<-raster("E:\\SIT\\caprese\\NPP\\hanpp_europe\\RNPPin.tif")
  Npp<-raster("E:\\SIT\\caprese\\NPP\\NPP_modis_2000-10prj.tif")		###modis 2000-2009 1 km

NFERT<-raster("E:\\SIT\\fra\\db\\map\\Ntot_fert.tif")

ERO<-raster("E:\\SIT\\caprese\\soil\\erosion\\rusle2015\\eros_EU28_1k.tif")
PRM<-raster("E:\\SIT\\caprese\\soil\\geology\\parmatdom2_prj.tif")

PET<-raster("E:\\SIT\\caprese\\meteo\\PET\\global-et0_annual\\ET0_EU_prj.tif")

WTngb<-raster('E:\\SIT\\caprese\\soil\\moisture_WT\\WTD_ngb.tif')
# SMwi<-raster('E:\\SIT\\caprese\\soil\\moisture_WT\\SMs/Season/SM_mean_DecFeb.tif')
# SMsp<-raster('E:\\SIT\\caprese\\soil\\moisture_WT\\SMs/Season/SM_mean_MarMay.tif')
# SMsu<-raster('E:\\SIT\\caprese\\soil\\moisture_WT\\SMs/Season/SM_mean_JunAug.tif')
# SMau<-raster('E:\\SIT\\caprese\\soil\\moisture_WT\\SMs/Season/SM_mean_SepNov.tif')

###extract#########

Corg<-extract(omad[[1]], lucasSHP,  method='bilinear', na.rm=T)
CNm<-extract(omad[[2]], lucasSHP,  method='bilinear', na.rm=T)
 irri<-extract(irr, lucasSHP,  method='bilinear', na.rm=T)
 
Ndep_WD_tx<-extract(NDEP, lucasSHP,  method='simple', na.rm=T)

MAT<- extract(WCLIM[[1]], lucasSHP, method='bilinear', na.rm=T)
RAIN<- extract(WCLIM[[2]], lucasSHP, method='bilinear', na.rm=T)
# EOBid<- extract(EOBS, lucasSHP, method='simple', na.rm=T)

NPP<-extract(Npp, lucasSHP, method='simple', na.rm=T)/1000			###ngb
NPPin<-extract(RNPPin_, lucasSHP, method='bilinear', na.rm=T)

Nfert<-extract(NFERT, lucasSHP,  method='bilinear', na.rm=T)
 Nfert<-ifelse(is.na(Nfert), 0, Nfert)

EROS<-extract(ERO, lucasSHP, method='bilinear', na.rm=T)
PRMT<-as.factor(extract(PRM, lucasSHP, method='simple', na.rm=T))

ET0<-extract(PET, lucasSHP, method='bilinear', na.rm=T)

WT<-extract(WTngb, lucasSHP, method='bilinear', na.rm=T)
# SMwn<-extract(SMwi, lucasSHP, method='bilinear', na.rm=T)
# SMsm<-extract(SMsu, lucasSHP, method='bilinear', na.rm=T)

###binding#########

lucasDB<-cbind(lucasDB, MAT, RAIN, Corg, CNm, irri, Ndep_WD_tx, NPP, EROS, Nfert, PRMT, WT) 


CN<-lucasDB$OC/lucasDB$N
NP<-lucasDB$N/lucasDB$P

 lucasDB<-cbind(lucasDB, "CN"=CN, "NP"=NP)



###############################################   LU      ########################################################
LU<-ifelse(lucasDB$LC09 %in% c("F00", "A22"), as.character(lucasDB$LC12), as.character(lucasDB$LC09))								##bare soil to LC15 if any crop
LU<-ifelse(substring(LU,1,2) %in% c("B1","B2","B3","B4","B5","B7","B8"), "CR", as.character(LU))
 LU<-ifelse(LU %in% c("E10", "E20", "E30"), "GR", as.character(LU))
 LU<-ifelse(LU %in% c("C10", "C20", "C30"), "FR", as.character(LU))
 LU<-ifelse(LU %in% c("D10", "D20"), "SHR", as.character(LU))

lucasDB<-cbind(lucasDB, "LU"=LU)													                  ##LU is the LC1 of 2009!!!!!!!!!!



OC_ul<-120; N_ul<-OC_ul/5
lucasDB<-subset(lucasDB, lucasDB$OC > 2 & lucasDB$OC <OC_ul)   				########limit detection OC to 120 g/gk
lucasDB<-subset(lucasDB, (lucasDB$clay+lucasDB$sand) %in% (1:99))   			#########removing mssing texture and OS!!!!!
# lucasDB<-subset(lucasDB, LU %in% c("SHR", "GR", "FR", "CR"))		  				#####shape predicted dataset!!!!!!!!!!
# lucasDB<-subset(lucasDB, LU %in% c("SHR", "E10", "E20", "E30", "C10", "C20", "C30", "CR"))				   				

 lucasDB$LU<- droplevels(lucasDB$LU)

#write.dbf(lucasDB, file="E:\\SIT\\fra\\db\\fractions\\lucasDB.dbf")

##################################################################################merging fractions
##################################################################################
#LCS_frac<-read.csv("E:\\SIT\\fra\\db\\fractions\\forR.csv", header = TRUE)
LCS_frac<-read.csv("E:\\SIT\\fra\\db\\fractions\\forR_allLC_CaCO3.csv", header = TRUE)		#!!!!!!!

#c0<-(LCS_frac$s_c_prc*10-LCS_frac$OC_sc_g_kg)/10+((1000-LCS_frac$s_c_prc*10)-LCS_frac$OC_pom_g_kg)/10
#c1<-LCS_frac$s_c_prc/((LCS_frac$s_c_prc*10-LCS_frac$OC_sc_g_kg)/10)
#c2<-(100-LCS_frac$s_c_prc)/(((1000-LCS_frac$s_c_prc*10)-LCS_frac$OC_pom_g_kg)/10)

#LCS_frac$OC_sc_g_kg<-LCS_frac$OC_sc_g_kg*c1
#LCS_frac$OC_pom_g_kg<-LCS_frac$OC_pom_g_kg*c2
#LCS_frac$s_c_prc<-c0

#plot( (LCS_frac$OC_sc_g_kg+LCS_frac$OC_pom_g_kg),LCS_frac$OC_pom_g_kg)
#points( (LCS_frac$OC_sc_g_kg+LCS_frac$OC_pom_g_kg), LCS_frac$OC_pom_g_kg, col='red', ylim=c(0,200))

#plot( (LCS_frac$OC_sc_g_kg+LCS_frac$OC_pom_g_kg), LCS_frac$OC_sc_g_kg, ylim=c(0,200))
#points( (LCS_frac$OC_sc_g_kg+LCS_frac$OC_pom_g_kg), LCS_frac$OC_sc_g_kg, col='red', ylim=c(0,200))



LCS<-merge(lucasDB, LCS_frac, by = "sample_ID", sort=FALSE, all = TRUE)

 LCSpred<-subset(LCS, is.na(LCS$OC_sc_g_kg))
LCS<-subset(LCS, LCS$OC_sc_g_kg>0)
LCS<-subset(LCS, LCS$OC>0)






#RANDOM FOREST
###########################################################################LCS correction
LCS$EROS[176]<-0.32
############################################################################### RF prediction MAOM
#pairs(OC_sc_g_kg ~ s_c_prc + pH_in_H2O + OC_tf + Ndep_WD_tx + MAT + EROS + LU, data= LCS,cex=0.9, col=alpha("black", 0.5))

set.seed(419)
RF<-randomForest(OC_sc_g_kg ~ s_c_prc + pH_in_H2O + OC_tf + Ndep_WD_tx + MAT + EROS + WT, data = LCS, ntree=1000, mtry=4, importance=TRUE, na.action=na.omit)
RF

#importance(RF)
#plot(RF)
varImpPlot(RF)
windows()

par(mfrow=c(3,3))
 partialPlot(RF, LCS, s_c_prc, main=paste0("%IncMSE=", round(importance(RF)[1],1)), xlab="silt+clay")
 partialPlot(RF, LCS, pH_in_H2O, main=paste0("%IncMSE=", round(importance(RF)[2],1)), xlab="pH")
 partialPlot(RF, LCS, OC_tf, main=paste0("%IncMSE=", round(importance(RF)[3],1)), xlab="soil OC")
 partialPlot(RF, LCS, Ndep_WD_tx, main=paste0("%IncMSE=", round(importance(RF)[4],1)), xlab="N deposition")
 partialPlot(RF, LCS, MAT, main=paste0("%IncMSE=", round(importance(RF)[5],1)), xlab="MAT")
 partialPlot(RF, LCS, EROS, main=paste0("%IncMSE=", round(importance(RF)[6],1)), xlab="erosion")
 partialPlot(RF, LCS, WT, main=paste0("%IncMSE=", round(importance(RF)[7],1)), xlab="WT")


predicted<-predict(RF, LCS)

#RMSEt<-paste("RMSE= ", round(rmse(LCS$OC_sc_g_kg, predicted, na.rm=TRUE),2), sep="")
F2<-qplot(LCS$OC_sc_g_kg, predicted,  size = I(2), alpha=I(1/5), color=I("red")) + 
 labs(x = "measured" , y =  "predicted") +
 scale_x_continuous(limits=c(0, 80)) +
 scale_y_continuous(limits=c(0, 80)) +
 geom_abline(intercept = 0, slope = 1)+
theme(text=element_text(size=16))+
geom_smooth(method='lm',formula=y~x) +
  theme_bw()

F2

plot((LCS$OC_tf-predicted), LCS$OC_pom_g_kg)


########resampling LCS for error propagation

set.seed(458)
LCS1<-LCS[sample(nrow(LCS), replace=TRUE),]
set.seed(258)
LCS2<-LCS[sample(nrow(LCS), replace=TRUE),]
set.seed(158)
LCS3<-LCS[sample(nrow(LCS), replace=TRUE),]

RF1<-randomForest(OC_sc_g_kg ~ s_c_prc + pH_in_H2O + OC_tf + Ndep_WD_tx + MAT + EROS + WT, data = LCS1, ntree=1000, mtry=4, importance=TRUE, na.action=na.omit)
RF2<-randomForest(OC_sc_g_kg ~ s_c_prc + pH_in_H2O + OC_tf + Ndep_WD_tx + MAT + EROS + WT, data = LCS2, ntree=1000, mtry=4, importance=TRUE, na.action=na.omit)
RF3<-randomForest(OC_sc_g_kg ~ s_c_prc + pH_in_H2O + OC_tf + Ndep_WD_tx + MAT + EROS + WT, data = LCS3, ntree=1000, mtry=4, importance=TRUE, na.action=na.omit)



##################################################N maom
#pairs(N_sc_g_kg ~ s_c_prc + pH_in_CaCl + N_tf + Ndep_WD_tx + MAT + LU + NP , data= LCS,cex=0.9, col=alpha("black", 0.5))
#RFn<-randomForest(N_sc_g_kg ~ s_c_prc + pH_in_H2O + N_tf + OC_tf + RAIN + Ndep_WD_tx + CEC, data = LCS, ntree=1000, mtry=5, importance=TRUE, na.action=na.omit)

set.seed(432)
RFn<-randomForest(N_sc_g_kg ~ s_c_prc + pH_in_H2O + N_tf + Ndep_WD_tx + MAT + EROS + WT, data = LCS, ntree=1000, mtry=4, importance=TRUE, na.action=na.omit)
RFn

#importance(RFn)
plot(RFn)
varImpPlot(RFn)


par(mfrow=c(3,3))
 partialPlot(RFn, LCS, s_c_prc, main=paste0("%IncMSE=", round(importance(RFn)[1],1)), xlab="silt+clay")
 partialPlot(RFn, LCS, pH_in_H2O, main=paste0("%IncMSE=", round(importance(RFn)[2],1)), xlab="pH")
 partialPlot(RFn, LCS, N_tf, main=paste0("%IncMSE=", round(importance(RFn)[3],1)), xlab="soil N")
 partialPlot(RFn, LCS, Ndep_WD_tx, main=paste0("%IncMSE=", round(importance(RFn)[4],1)), xlab="N deposition")
 partialPlot(RFn, LCS, MAT, main=paste0("%IncMSE=", round(importance(RFn)[5],1)), xlab="MAT")
 partialPlot(RFn, LCS, EROS, main=paste0("%IncMSE=", round(importance(RFn)[6],1)), xlab="erosion")
 partialPlot(RFn, LCS, WT, main=paste0("%IncMSE=", round(importance(RFn)[7],1)), xlab="WT")


predictedn<-predict(RFn, LCS)

#RMSEt<-paste("RMSE= ", round(rmse(LCS$N_sc_g_kg, predictedn, na.rm=TRUE),2), sep="")
F3<-qplot(LCS$N_sc_g_kg, predictedn,  size = I(2), alpha=I(1/5), color=I("red")) + 
 labs(x = "measured" , y =  "predicted") +
 scale_x_continuous(limits=c(0, 8)) +
 scale_y_continuous(limits=c(0, 8)) +
 geom_abline(intercept = 0, slope = 1)+
theme(text=element_text(size=16)) +
geom_smooth(method='lm',formula=y~x)    

F3



RFn1<-randomForest(N_sc_g_kg ~ s_c_prc + pH_in_H2O + N_tf + Ndep_WD_tx + MAT + EROS + WT, data = LCS1, ntree=1000, mtry=4, importance=TRUE, na.action=na.omit)
RFn2<-randomForest(N_sc_g_kg ~ s_c_prc + pH_in_H2O + N_tf + Ndep_WD_tx + MAT + EROS + WT, data = LCS2, ntree=1000, mtry=4, importance=TRUE, na.action=na.omit)
RFn3<-randomForest(N_sc_g_kg ~ s_c_prc + pH_in_H2O + N_tf + Ndep_WD_tx + MAT + EROS + WT, data = LCS3, ntree=1000, mtry=4, importance=TRUE, na.action=na.omit)





