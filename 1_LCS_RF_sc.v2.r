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
WCLpr<-stack("E:\\SIT\\caprese\\meteo\\worldclim2\\HADes8.5_2080_LAEA.tif")

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
EOBid<- extract(EOBS, lucasSHP, method='simple', na.rm=T)

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

LCS<-merge(lucasDB, LCS_frac, by = "sample_ID", sort=FALSE, all = TRUE)

 LCSpred<-subset(LCS, is.na(LCS$OC_sc_g_kg))
LCS<-subset(LCS, LCS$OC_sc_g_kg>0)
LCS<-subset(LCS, LCS$OC>0)

#################################################LCS correction
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




##########################################################################################################CLC12
CLC<-stack("E:\\SIT\\caprese\\ven\\sripts_L\\exp\\_LINUXSERVER_soil\\VCFR\\SPT_lyr\\CLC2.tif")[[4]]

f <- function(x)  ifelse(x %in% c(12,13,seq(15,29)), x, NA)              ###match LUCAS samplng points for LU_LC
CLC <- calc(CLC, fun=f)

LUmask<-CLC*0+1

CLC12class<-read.csv("E:/SIT/caprese/CLC/CLC12/Legend/clc_legend.csv", sep=";")

#############################################################################################################################################
##layer n1
sand<-stack('E:/SIT/caprese/ven/sripts_L/exp/_LINUXSERVER_soil/VCFR/SPT_lyr/soil_LCS.tif')[[1]]
 silt_clay<- 100-sand                                                                            

##layer n2 
pH_H2O<-stack('E:/SIT/caprese/ven/sripts_L/exp/_LINUXSERVER_soil/VCFR/SPT_lyr/soil_LCS.tif')[[8]] /10
    

##layer n3 SOC
SOC20st_1k<-raster('E:/SIT/caprese/soil/SOC_LUCAS/export/OCstk20_c_bl.tif')/1000     ###kg/hh to Mg/ha
BDm_c_1k<-raster('E:/SIT/caprese/soil/SOC_LUCAS/export/BDm_c_bl.tif')  

SOC20st_1k<-resample(SOC20st_1k, silt_clay) 
BDm_c_1k<-resample(BDm_c_1k, silt_clay) 

SOC20c_1k<-(SOC20st_1k)/(BDm_c_1k*2000)*1000                  #SOC concentrario g/kg

excl<-SOC20c_1k; excl[excl>OC_ul]<-NA ; excl<-excl*0+1        #exclusion >120 g/kg

SOC20c_1k<-SOC20c_1k*excl*LUmask                              #exclusion for LC_LU CLC-LUCAS  

##layer n3 N (g/kg)
N_1k<-raster('//ies.jrc.it/H05/SOIL/ESDAC/MOSES data/LUCAS-2009-topsoil-chemical-properties/data/N/N.tif')
 N_1k<-resample(N_1k, silt_clay)  


##layer n4 
NDEP1k<-resample(NDEP, silt_clay)  

##layer n5 
MAT1k<-resample(WCLIM[[1]], silt_clay)  
RAIN1k<-resample(WCLIM[[2]], silt_clay)  
MAT1kpr<-resample(WCLpr[[1]], silt_clay)  
RAIN1kpr<-resample(WCLpr[[2]], silt_clay) 

##layer n6 
EROS1k<-ERO#resample(ERO, silt_clay)  

##layer n7 
WTngb


###SOC###
#RF s_c_prc + pH_in_CaCl + OC_tf + Ndep_WD_tx + MAT + EROS + K


inputR<-silt_clay															
inputR<-addLayer(inputR, pH_H2O, SOC20c_1k, NDEP1k, MAT1k, EROS1k, WTngb)			
names(inputR)<-c('s_c_prc', 'pH_in_H2O',  'OC_tf',  'Ndep_WD_tx',  'MAT', 'EROS', 'WT') 


maom_1k<- predict(inputR, RF, type='response', progress='window', na.rm=T)
pom_1k<- SOC20c_1k - maom_1k
m_t_SOCr<-maom_1k/SOC20c_1k


 # writeRaster(maom_1k, filename="maom_1k.tif", format="GTiff",overwrite=T)
 # writeRaster(pom_1k, filename="pom_1k.tif", format="GTiff",overwrite=T)
 # writeRaster(m_t_SOCr, filename="m_t_SOCr_1k.tif", format="GTiff",overwrite=T)


###soil N### 
 
inputR[[3]]<-N_1k 
names(inputR)<-c('s_c_prc', 'pH_in_H2O',  'N_tf',  'Ndep_WD_tx',  'MAT', 'EROS', 'WT') 
 
 maoN_1k<- predict(inputR, RFn, type='response', progress='window', na.rm=T)
 poN_1k<- N_1k - maoN_1k
 
 zonal(maom_1k/maoN_1k, CLC, 'mean', na.rm=T) 	
 zonal(pom_1k/poN_1k, CLC, 'mean', na.rm=T) 	

                                   

###################################################################################scatter plot MAOM vs SOC
library(reshape2)
 library(gridExtra)
  library(rasterVis)
   library(quantreg)
    library(sf)
#    library(maptools)
 
SOCf<-cbind(as.vector(SOC20c_1k), as.vector(maom_1k))
SOCf<-as.data.frame(SOCf)


FIG_sat<-ggplot(data=SOCf, aes(V1,V2)) + 
  geom_hex(bins=100, binwidth=c(1.5,1.5)) +
  scale_fill_gradientn(colours=c('light gray','blue'),name='Frequency',na.value=NA, trans = "sqrt") +
  theme_bw()+
  scale_y_continuous(limits=c(0, OC_ul))+
  scale_x_continuous(limits=c(0, OC_ul))+
  geom_abline(intercept = 0, slope = 1) +
  xlab("SOC") + ylab("MAOM") +
  geom_point(data=LCS, aes(x=OC_tf , y=OC_sc_g_kg), inherit.aes = FALSE, size=0.1, alpha = 0.4) 

                     
grid.arrange(F2,FIG_sat,nrow =1)


MAOM_LU<-zonal(maom_1k*BDm_c_1k*cm_prof*0.1 , CLC, 'sum', na.rm=T)    ###SOC
 POM_LU<-zonal(pom_1k*BDm_c_1k*cm_prof*0.1 , CLC, 'sum', na.rm=T) 
 
MAON_LU<-zonal(maoN_1k*BDm_c_1k*cm_prof*0.1 , CLC, 'sum', na.rm=T)    ###N	
 PON_LU<-zonal(poN_1k*BDm_c_1k*cm_prof*0.1 , CLC, 'sum', na.rm=T)  
 
 
  count<-zonal(maom_1k*0+1 , CLC, 'sum', na.rm=T)
  fraction<-as.data.frame(cbind(POM_LU, MAOM_LU[,2], count[,2],  PON_LU[,2], MAON_LU[,2]))
  names(fraction)<-c("GRID_CODE", "POM", "MAOM", "area", "PON", "MAON")

  fraction<-merge(CLC12class, fraction, by.x = "GRID_CODE", by.y = "GRID_CODE")
   fraction$LC<-as.factor(c('arable', 'arable', 'per. crop', 'per. crop', 'per. crop', 'grassland', 
                'heterog. agr.','heterog. agr.','heterog. agr.','heterog. agr.', 'broad l.', 'conif.', 'mixed f.',
                'shrub-herb.', 'shrub-herb.', 'shrub-herb.', 'shrub-herb.')) 

#################monetary value########################
#cur_p<-25  	##price in Euro per t of CO2
#p_mrt<-4	##ratio mrt

#tot_e<-cur_p*3.67*(sum(fraction$MAOM)+sum(fraction$POM))/10000000         ###billion euro

#pom_p<-tot_e/((sum(fraction$MAOM)*p_mrt+ sum(fraction$POM))*3.67/10000000)	 ###euro per t of CO2 
#maom_p<-pom_p*p_mrt

######################################################

fract_r<-aggregate(. ~ LC, fraction[,5:10], sum)

   fract_r<- melt(fract_r, id=c("LC", "area"))
   value_N<-fract_r$value[17:32]
   fract_r<-fract_r[1:16,]
   fract_r$value_N<- value_N
   fract_r$area_n<-fract_r[,2]/max(fract_r[,2])
   fract_r$mean<-fract_r$value/fract_r$area
   fract_r$CN<-fract_r$value/fract_r$value_N
   #fract_r$price<-ifelse(fract_r$variable=="POM",fract_r$value*3.67*pom_p/10000000,fract_r$value*3.67*maom_p/10000000)
   #fract_r$priceHa<-ifelse(fract_r$variable=="POM",fract_r$mean*3.67*pom_p,fract_r$mean*3.67*maom_p)
     


   
 
p<-ggplot(data=fract_r, aes(x=reorder(LC,value), y=value/10000, fill=variable, width=area_n)) +
   geom_bar(stat="identity") + 
   ylab("SOC Mt") + xlab("") +
  geom_text(aes(label=round(value_N/10000,0)), position = position_stack(vjust = 0.5), size=2.8) +
   coord_flip() +
  theme(legend.position = c(.85, .20),
        legend.title = element_blank(), 
       panel.background = element_rect(fill = "white", colour = "grey50"))+
annotate("text", x = 1, y = 5000, label = "N Mt")
 
p1<-ggplot(data=fract_r, aes(x=reorder(LC,mean), y=mean, fill=variable, width=area_n)) +
   geom_bar(stat="identity") +
  ylab(expression("SOC Mg ha"^-1)) + xlab("") +
geom_text(aes(label=round(CN,0)), position = position_stack(vjust = 0.5), size=2.8) +
  coord_flip() +
  theme(legend.position = "none", 
  panel.background = element_rect(fill = "white", colour = "grey50"))+
annotate("text", x = 1, y = 90, label = "C:N ratio") 


#p2<-ggplot(data=fract_r, aes(x=reorder(LC,price), y=price, fill=variable, width=area_n)) +
#   geom_bar(stat="identity") +
#  ylab("Billion Euro") + xlab("") +
#geom_text(aes(label=round(priceHa/1000,1)), position = position_stack(vjust = 0.5), size=2.8) +
#  coord_flip() +
#  theme(legend.position = "none", 
#  panel.background = element_rect(fill = "white", colour = "grey50"))+
# annotate("text", x = 1, y = 500, label = expression("keuro ha"^-1))


#grid.arrange(p,p1,p2, nrow =1)
grid.arrange(p1,p, nrow =1)
 
################################################################################# MAOM saturation
 # f <- function(x1, x2)  ifelse(x1 %in% c(18,23,24,25,26), x2, NA)
 # MAOM_sat <- overlay(CLC, maom_1k, fun=f, forcefun=TRUE)
 # 
 # MAOM_sat2<- focal(MAOM_sat, w=matrix(1,nrow=5,ncol=5), fun=max, na.rm=T) 
 # 
 # writeRaster(MAOM_sat2-maom_1k, filename="MAOM_def.tif", format="GTiff",overwrite=T)
 
 ################################################################################# MAOM sequestration by2050
#  base<-stack("E:/model/daycent/sim/lucas/results/1_km/_JEODPP/SOC_EU_DayC_grv.tif")[[6]]
#   CC<-stack("E:/model/daycent/sim/lucas/results/1_km/_JEODPP/SOC_EU_DayC_grv_CC.tif")[[6]]
#    CA<-stack("E:/model/daycent/sim/lucas/results/1_km/_JEODPP/SOC_EU_DayC_grv_ha_CA.tif")[[6]]
#    
#    
#    dCC<- CC-base ; dCC[dCC==0]<-NA ; dCC<-dCC/1.5
#    dCA<- CA-base ; dCA[dCA==0]<-NA ; dCA<-dCA/1.5
#  
#    SEQmax <- overlay(dCC, dCA, fun=max)
#    
#    SOC20c_BM<-((SOC20st_1k+SEQmax)/(BDm_c_1k*2000))*1000             ###Cseq + stock to content g/kg
#    
#    inputR<-silt_clay															
#    inputR<-addLayer(inputR, pH_H2O, SOC20c_BM, NDEP1k, MAT1k, EROS1k, K)			
#    names(inputR)<-c('s_c_prc', 'pH_in_H2O',  'OC_tf',  'Ndep_WD_tx',  'MAT', 'EROS', 'K')   
#    
#    
#    
#    maom_BM<- predict(inputR, RF, type='response', progress='window', na.rm=T)  
#    pom_BM<- SOC20c_BM - maom_BM
#    
#    dmaom_BM<- (maom_BM - maom_1k) * BDm_c_1k*cm_prof*0.1
#    dpom_BM<-  (pom_BM - pom_1k) * BDm_c_1k*cm_prof*0.1
#    
#    cellStats(dmaom_BM, sum, na.rm=T); cellStats(dpom_BM, sum, na.rm=T)
#    
# #   writeRaster(dmaom_BM*BDm_c_1k*cm_prof*0.1, filename="dmaom_BM_t_ha.tif", format="GTiff",overwrite=T)   ###stock BM
# #   writeRaster(dpom_BM*BDm_c_1k*cm_prof*0.1, filename="dpom_BM_t_ha.tif", format="GTiff",overwrite=T)   ###stock BM
#    
#  
# ############################################price
# seq_price<- (dmaom_BM*maom_p + dpom_BM*pom_p)*3.67/30	
# # seq_price_<- (dmaom_BM*25 + dpom_BM*25)*3.67/30	
# 
# dSQP<- seq_price - seq_price_

############################################FIGURES
dem<-raster("E:\\SIT\\caprese\\soil\\dem\\hillshade1x1LAEA.tif") 
   map<-dmaom_BM   
   e <- extent(dmaom_BM)
   dem<-crop(dem,e)
   
mapaSHP <- rgdal::readOGR('E:\\SIT\\caprese\\adm_units\\country_def.shp')

t0<-levelplot(dem*0, par.settings = GrTheme, maxpixels=1e7)

####fig1

soc_lab=expression("g C kg"^-1*"soil")

t1 <- levelplot(maom_1k, par.settings = viridisTheme, margin = F, maxpixels=1e7, 
scales = list(draw = FALSE), 
colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label=soc_lab, cex=0.8))
t1 + as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))


t2 <- levelplot(m_t_SOCr, at= c(seq(0, 1, 0.1)), par.settings = RdBuTheme, margin = F, maxpixels=1e7, 
scales = list(draw = FALSE), 
colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label="MAOM:SOC", cex=0.8))
t2 + as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))




####fig3
# fig3<-stack(seq_price, SEQmax/30)
# 
# t3 <- levelplot(fig3, layers=1, par.settings = magmaTheme, margin = F, maxpixels=1e7, 
# scales = list(draw = FALSE), 
# colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label=expression("Euro ha"^-1*"yr"^-1), cex=0.8))
# t3 + as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))
# 
# 
# 
# t4 <- levelplot(fig3, layers=2,  par.settings = viridisTheme, margin = F, maxpixels=1e7, 
# scales = list(draw = FALSE), 
# colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label=expression("Mg C ha"^-1*"yr"^-1), cex=0.8))
# t4 + as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))



###################point distribution
#####################################################################################################
##Projection LAEA
#####################################################################################################
#europe_laea<-sf::st_read("E:\\SIT\\caprese\\adm_units\\country_def.shp")
WGS<-"+proj=longlat +datum=WGS84 +no_defs"
LAEA<-"+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

world_wgs<-sf::st_read("E:\\SIT\\caprese\\adm_units\\ne_50m_admin_0_countries\\ne_50m_admin_0_countries.shp")
countries_laea<-st_transform(world_wgs, LAEA)


data<-LCS
coordinates(data)<-c("GPS_LONG","GPS_LAT")
proj4string(data) <- CRS(WGS)
data<-spTransform(data, CRS(LAEA))
####supplementary CLC
  
CLCr<-ratify(CLC) 
  
rat <- levels(CLCr)[[1]]
rat$landcover <- fraction$LABEL3
rat$landc <- fraction$LC
rat$CLCnew<-paste0(rat$landcover,"/",rat$landc)
rat$CLCnew[7]<-"Annual crop ass. permanent crops/heterog. agr." 
rat$CLCnew[9]<-"Agricultural and natural vegetation/heterog. agr." 
levels(CLCr) <- rat



t5<-levelplot(CLCr, att='CLCnew', margin = F, maxpixels=1e7,
 scales = list(draw = FALSE), 
 col.regions=c('darkgoldenrod1','darkgoldenrod1','chocolate3','chocolate3','chocolate3',
 'chartreuse3','gold4','gold4','gold4','gold4',
 'forestgreen','darkgreen','darkolivegreen','burlywood4','burlywood4',
 'burlywood4','burlywood4'))
t5 + layer(sp.lines(mapaSHP, lwd=0.1, col='darkgray')) + layer(sp.points(data, lwd=0.1, col='black')) 

#writeRaster(CLCr, filename="CLCr.tif", format="GTiff",overwrite=T)


####supplementary SOC
t8 <- levelplot(SOC20c_1k, par.settings = magmaTheme, margin = F, maxpixels=1e7, 
scales = list(draw = FALSE), 
colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label=soc_lab, cex=0.8))
t8 + as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))



####supplementary POM
t6 <- levelplot(pom_1k, par.settings = viridisTheme, margin = F, maxpixels=1e7, 
scales = list(draw = FALSE), 
colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label=soc_lab, cex=0.8))
t6 + as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))



####supplementary seq CC
#t7 <- levelplot(dCC/30, at=seq( -0.04, 0.38, 0.02), par.settings = magmaTheme, margin = F, maxpixels=1e7, 
#scales = list(draw = FALSE), 
#colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label=expression("Mg C ha"^-1*"yr"^-1), cex=0.8))
#t7 + as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))


####supplementary seq price difference
#t8 <- levelplot(dSQP, par.settings = RdBuTheme, margin = F, maxpixels=1e7, 
#scales = list(draw = FALSE), 
#colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label=expression("Euro ha"^-1*"yr"^-1), cex=0.8))
#t8 + as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))

quantile(dSQP, probs = c(0.1, 0.5, 0.9), type=7, names = T)




####supplementary LCS texture distribution 

LCSall<-rbind(LCS, LCSpred)

##ver.3.5
tex <- data.frame(
"CLAY" = LCSall$clay,
"SILT" = 100-LCSall$clay-LCSall$sand,
"SAND"=LCSall$sand
)


library(ggtern)
ggtern(data=tex,aes(SILT,CLAY, SAND)) + geom_point(alpha=0.2)+
geom_point(data = tex[1:352,], col = 'red', size=1.5)


 hist(LCSall$OC,  main='', xlab="SOC (g/ka)", freq=FALSE, lty=1)
 hist(LCS$OC_tf,  add=T, col=alpha("red", 0.2), freq=FALSE,  lty=0 )


############################################################################### 
############################################################################### 
LCS<-subset(LCS, LU %in% c("SHR", "GR", "FR", "CR"))
LCS$LU<- droplevels(LCS$LU)


lm1<-ggplot(data=LCS, aes(MAT, OC_pom_g_kg, colour = RAIN)) + geom_point(alpha=0.5, size=1.8) + ylab("POM")+ 
facet_wrap(~LU)+
geom_smooth(method='lm',formula= y ~ x)+
scale_colour_gradientn(colours = terrain.colors(10))


lm1<-ggplot(data=LCS, aes(MAT, OC_sc_g_kg, colour = RAIN)) + geom_point(alpha=0.5, size=1.8) + ylab("MAOM")+ 
facet_wrap(~LU)+
geom_smooth(method='lm',formula= y ~ x)+
scale_colour_gradientn(colours = terrain.colors(10))

dT<- MAT1kpr-MAT1k
dP<- RAIN1kpr-RAIN1k


dC_LU<-data.frame('LU'= c('CR', 'FR', 'GR', 'SRH', 'CR', 'FR', 'GR', 'SRH'), 
'FRAC'= c('POM', 'POM', 'POM', 'POM', 'MAOM', 'MAOM', 'MAOM', 'MAOM'), 'dC'=rep(0, times = 8), 'area'=rep(0, times = 8))


#POM:j=0 and MAOM:j=4
j=4

if (j==0){
###POM
md0<-lm(LCS$OC_pom_g_kg~LCS$MAT*LCS$LU)
md1<-lm(LCS$OC_pom_g_kg~(LCS$MAT+LCS$RAIN)*LCS$LU)
}else{
###MAOM
md0<-lm(LCS$OC_sc_g_kg~LCS$MAT*LCS$LU)
md1<-lm(LCS$OC_sc_g_kg~(LCS$MAT+LCS$RAIN)*LCS$LU)
}


STcr<-md1$coefficients[2]
STfr<-md1$coefficients[2]+md1$coefficients[7]
STgr<-md1$coefficients[2]+md1$coefficients[8]
STsh<-md1$coefficients[2]+md1$coefficients[9]

SPcr<-md1$coefficients[3]
SPfr<-md1$coefficients[3]+md1$coefficients[10]
SPgr<-md1$coefficients[3]+md1$coefficients[11]
SPsh<-md1$coefficients[3]+md1$coefficients[12]


f <- function(x1, x2, x3)  ifelse(x1 %in% c(12,13,15,16,17,19,20,21,22), x2*STcr + x3*SPcr, NA)		#CR
 CRdt <- overlay(CLC, dT, dP, fun=f, forcefun=TRUE)
 CRdt <- CRdt*(pom_1k*0+1)

f <- function(x1, x2, x3)  ifelse(x1 %in% c(23,24,25), x2*STfr + x3*SPfr, NA)		#FR
 FRdt <- overlay(CLC, dT, dP, fun=f, forcefun=TRUE)
 FRdt <- FRdt*(pom_1k*0+1)
  
f <- function(x1, x2, x3)  ifelse(x1 %in% c(18), x2*STgr + x3*SPgr, NA)		#GR
 GRdt <- overlay(CLC, dT, dP, fun=f, forcefun=TRUE)
 GRdt <- GRdt*(pom_1k*0+1)

f <- function(x1, x2, x3)  ifelse(x1 %in% c(26,27,28,29), x2*STsh + x3*SPsh, NA)		#SH
 SHdt <- overlay(CLC, dT, dP, fun=f, forcefun=TRUE)
 SHdt <- SHdt*(pom_1k*0+1)


dMT<-merge(CRdt,FRdt,GRdt,SHdt)
dMT<-dMT*BDm_c_1k*2				#t/ha of C
cellStats(dMT, sum, na.mr=T)*100/1000000	
tar<-cellStats(dMT*0+1, sum, na.mr=T)

dC_LU[(1+j),3]<-cellStats(CRdt*BDm_c_1k*2, sum, na.mr=T)*100/1000000	
dC_LU[(2+j),3]<-cellStats(FRdt*BDm_c_1k*2, sum, na.mr=T)*100/1000000
dC_LU[(3+j),3]<-cellStats(GRdt*BDm_c_1k*2, sum, na.mr=T)*100/1000000
dC_LU[(4+j),3]<-cellStats(SHdt*BDm_c_1k*2, sum, na.mr=T)*100/1000000

dC_LU[(1+j),4]<-cellStats(CRdt*0+1, sum, na.mr=T)/tar	
dC_LU[(2+j),4]<-cellStats(FRdt*0+1, sum, na.mr=T)/tar
dC_LU[(3+j),4]<-cellStats(GRdt*0+1, sum, na.mr=T)/tar
dC_LU[(4+j),4]<-cellStats(SHdt*0+1, sum, na.mr=T)/tar




#-20, 20, 5
t4 <- levelplot(dMT, at= c(seq(-20, 20, 5)), par.settings = RdBuTheme, margin = F, maxpixels=1e7, 
scales = list(draw = FALSE), 
colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label="dMAOM (Mt of C)", cex=0.8))
t4 + as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))



p<-ggplot(data=dC_LU, aes(x=reorder(LU,dC), y=dC, fill=FRAC, width=V4)) +
  geom_bar(stat="identity")+
 ylab("OC Mt") + xlab("") +
theme(legend.position = c(.20, .90),
        legend.title = element_blank(), 
       panel.background = element_rect(fill = "white", colour = "grey50"))+
coord_flip()






