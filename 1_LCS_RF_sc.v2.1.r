library(reshape2)

#########################################################################################################
#RANDOM FOREST upscaling
##########################################################################################################CLC12
CLC<-stack("E:\\SIT\\caprese\\ven\\sripts_L\\exp\\_LINUXSERVER_soil\\VCFR\\SPT_lyr\\CLC2.tif")[[4]]

f <- function(x)  ifelse(x %in% c(12,13,seq(15,29)), x, NA)              ###match LUCAS samplng points for LU_LC
CLC <- calc(CLC, fun=f)

LUmask<-CLC*0+1

CLC12class<-read.csv("E:/SIT/caprese/CLC/CLC12/Legend/clc_legend.csv", sep=";")

#############################################################################################################################################
##layer n1
sand<-stack('E:/SIT/caprese/ven/sripts_L/exp/_LINUXSERVER_soil/VCFR/SPT_lyr/soil_LCS.tif')[[1]]  #JRC
 silt_clay<- 100-sand    

siltI<-raster('E:/SIT/caprese/soil/SOC_ISRIC/LEONIDAS/SoilGrids/silt_0-20.tif')   #ISRIC
clayI<-raster('E:/SIT/caprese/soil/SOC_ISRIC/LEONIDAS/SoilGrids/clay_0-20.tif')
 silt_clayI<- (siltI+clayI)/10 		                                                                     
 silt_clayI[silt_clayI==0]<-NA

##layer n2 
pH_H2O<-stack('E:/SIT/caprese/ven/sripts_L/exp/_LINUXSERVER_soil/VCFR/SPT_lyr/soil_LCS.tif')[[8]] /10	#JRC
pH_H2OI<-raster('E:/SIT/caprese/soil/SOC_ISRIC/LEONIDAS/SoilGrids/phh2o_0-20.tif')/10	#ISRIC
 pH_H2OI[pH_H2OI==0]<-NA
    

##layer n3 SOC
## JRC ##
SOC20st_1k<-raster('E:/SIT/caprese/soil/SOC_LUCAS/export/OCstk20_c_bl.tif')/1000     ###kg/hh to Mg/ha
BDm_c_1k<-raster('E:/SIT/caprese/soil/SOC_LUCAS/export/BDm_c_bl.tif')  

SOC20st_1k<-resample(SOC20st_1k, silt_clay) 
BDm_c_1k<-resample(BDm_c_1k, silt_clay) 

SOC20c_1k<-(SOC20st_1k)/(BDm_c_1k*2000)*1000                  #SOC concentrario g/kg
 #excl<-SOC20c_1k; excl[excl>OC_ul]<-NA ; excl<-excl*0+1       #exclusion >120 g/kg
  SOC20c_1k[SOC20c_1k>OC_ul]<-NA
  SOC20c_1k<-SOC20c_1k*LUmask                              #exclusion for LC_LU CLC-LUCAS  

## ISRIC ##
SOC20c_1kI<-raster('E:/SIT/caprese/soil/SOC_ISRIC/LEONIDAS/SoilGrids/soc_0-20.tif')/10
BD_1kI<-raster('E:/SIT/caprese/soil/SOC_ISRIC/LEONIDAS/SoilGrids/bdod_0-20.tif')/100
GRV_1kI<-raster('E:/SIT/caprese/soil/SOC_ISRIC/LEONIDAS/SoilGrids/cfvo_0-20.tif')/10

 SOC20c_1kI[SOC20c_1kI==0]<-NA

## limitation to mineral soil JRC########
  f <- function(x1, x2)  ifelse(x1 > 120, NA, x2)
  SOC20c_1k <- overlay(SOC20c_1kI, SOC20c_1k, fun=f, forcefun=TRUE)
##

 SOC20c_1kI[SOC20c_1kI>OC_ul]<-NA
 SOC20c_1kI<-SOC20c_1kI*LUmask

SOC20st_1kI<-SOC20c_1kI*BD_1kI*(1-GRV_1kI/100)*2



## BD stack for maom pom stock conversions ##
BD_1kE<-stack(BDm_c_1k, BD_1kI*(1-GRV_1kI/100), BDm_c_1k, BD_1kI*(1-GRV_1kI/100))
BD_1kE[BD_1kE==0]<-NA



##layer n3 N (g/kg)
N_1k<-raster('//ies.jrc.it/H05/SOIL/ESDAC/MOSES data/LUCAS-2009-topsoil-chemical-properties/data/N/N.tif')
 N_1k<-resample(N_1k, silt_clay)  
 N_1k<-N_1k*(SOC20c_1k*0+1)

N_1kI<-raster('E:/SIT/caprese/soil/SOC_ISRIC/LEONIDAS/SoilGrids/nitrogen_0-20.tif')/100
 N_1kI[N_1kI==0]<-NA
 N_1kI<-N_1kI*(SOC20c_1kI*0+1)


##layer n4 
NDEP1k<-resample(NDEP, silt_clay)  

##layer n5 
MAT1k<-resample(WCLIM[[1]], silt_clay)  
RAIN1k<-resample(WCLIM[[2]], silt_clay)  
 

##layer n6 
EROS1k<-ERO#resample(ERO, silt_clay)  

##layer n7 
WTngb


########################################## spatial predictions ##########################################

### SOC ###############

inputR<-silt_clay															
inputR<-addLayer(inputR, pH_H2O, SOC20c_1k, NDEP1k, MAT1k, EROS1k, WTngb)			
inputR1<-inputR; inputR1[[1]]<-silt_clayI
inputR2<-inputR; inputR2[[3]]<-SOC20c_1kI
inputR3<-inputR; inputR3[[1]]<-silt_clayI; inputR3[[3]]<-SOC20c_1kI
 
namesR<-function(inputR) {
names(inputR)<-c('s_c_prc', 'pH_in_H2O',  'OC_tf',  'Ndep_WD_tx',  'MAT', 'EROS', 'WT') 
return(inputR)
}

inputR<-namesR(inputR); inputR1<-namesR(inputR1); inputR2<-namesR(inputR2); inputR3<-namesR(inputR3)


														#SOC		silt+clay
maom_1k<-  predict(inputR, RF, type='response', progress='window', na.rm=T)		#JRC		JRC
maom_1k1<- predict(inputR1, RF1, type='response', progress='window', na.rm=T)		#JRC 		ISRIC
maom_1k2<- predict(inputR2, RF2, type='response', progress='window', na.rm=T)		#ISRIC	JRC
maom_1k3<- predict(inputR3, RF3, type='response', progress='window', na.rm=T)		#ISRIC	ISRIC

pom_1k<-  SOC20c_1k -  maom_1k ;  m_t_SOCr<- maom_1k/SOC20c_1k	#POM estimation
pom_1k1<- SOC20c_1k -  maom_1k1 ; m_t_SOCr1<-maom_1k1/SOC20c_1k
pom_1k2<- SOC20c_1kI - maom_1k2 ; m_t_SOCr2<-maom_1k2/SOC20c_1kI
pom_1k3<- SOC20c_1kI - maom_1k3 ; m_t_SOCr3<-maom_1k3/SOC20c_1kI

maom_1kE<-stack(maom_1k, maom_1k1, maom_1k2, maom_1k3)	#stack
pom_1kE<-stack(pom_1k, pom_1k1, pom_1k2, pom_1k3)


 # writeRaster(maom_1k, filename="maom_1k.tif", format="GTiff",overwrite=T)
 # writeRaster(pom_1k, filename="pom_1k.tif", format="GTiff",overwrite=T)
 # writeRaster(SOC20c_1k, filename="SOC20c_1k.tif", format="GTiff",overwrite=T)


### soil N ###############

inputR[[3]]<-N_1k 														
inputR1<-inputR; inputR1[[1]]<-silt_clayI
inputR2<-inputR; inputR2[[3]]<-N_1kI
inputR3<-inputR; inputR3[[1]]<-silt_clayI; inputR3[[3]]<-N_1kI
 
namesR<-function(inputR) {
names(inputR)<-c('s_c_prc', 'pH_in_H2O',  'N_tf',  'Ndep_WD_tx',  'MAT', 'EROS', 'WT')
return(inputR)
}

inputR<-namesR(inputR); inputR1<-namesR(inputR1); inputR2<-namesR(inputR2); inputR3<-namesR(inputR3)

														#Nsoil	silt+clay
maoN_1k<-  predict(inputR, RFn, type='response', progress='window', na.rm=T)		#JRC		JRC
maoN_1k1<- predict(inputR1, RFn1, type='response', progress='window', na.rm=T)	#JRC 		ISRIC
maoN_1k2<- predict(inputR2, RFn2, type='response', progress='window', na.rm=T)	#ISRIC	JRC
maoN_1k3<- predict(inputR3, RFn3, type='response', progress='window', na.rm=T)	#ISRIC	ISRIC

poN_1k<-  N_1k -  maoN_1k 		#PON estimation
poN_1k1<- N_1k -  maoN_1k1 
poN_1k2<- N_1kI - maoN_1k2 
poN_1k3<- N_1kI - maoN_1k3 

maoN_1kE<-stack(maoN_1k, maoN_1k1, maoN_1k2, maoN_1k3)	#stack
poN_1kE<-stack(poN_1k, poN_1k1, poN_1k2, poN_1k3)


#inputR[[3]]<-N_1k 
#names(inputR)<-c('s_c_prc', 'pH_in_H2O',  'N_tf',  'Ndep_WD_tx',  'MAT', 'EROS', 'WT') 
 
# maoN_1k<- predict(inputR, RFn, type='response', progress='window', na.rm=T)
# poN_1k<- N_1k - maoN_1k
 
 #zonal(maom_1k/maoN_1k, CLC, 'mean', na.rm=T) 	
 #zonal(pom_1k/poN_1k, CLC, 'mean', na.rm=T) 	

                                   
MAOM_LU<-zonal(maom_1kE*BD_1kE*cm_prof*0.1 , CLC, 'sum', na.rm=T)	 ###SOC
 POM_LU<-zonal(pom_1kE*BD_1kE*cm_prof*0.1 , CLC, 'sum', na.rm=T) 
 
MAON_LU<-zonal(maoN_1kE*BD_1kE*cm_prof*0.1 , CLC, 'sum', na.rm=T)  ###N	
 PON_LU<-zonal(poN_1kE*BD_1kE*cm_prof*0.1 , CLC, 'sum', na.rm=T)  

MAOM_LUE<-cbind(MAOM_LU[,1], 'mean'=apply(MAOM_LU[,2:5], 1, mean, na.rm = TRUE), 'sd'= apply(MAOM_LU[,2:5], 1, sd, na.rm = TRUE)) 
 POM_LUE<-cbind(POM_LU[,1], 'mean'=apply(POM_LU[,2:5], 1, mean, na.rm = TRUE), 'sd'= apply(POM_LU[,2:5], 1, sd, na.rm = TRUE)) 

MAON_LUE<-cbind(MAON_LU[,1], 'mean'=apply(MAON_LU[,2:5], 1, mean, na.rm = TRUE), 'sd'= apply(MAON_LU[,2:5], 1, sd, na.rm = TRUE)) 
 PON_LUE<-cbind(PON_LU[,1], 'mean'=apply(PON_LU[,2:5], 1, mean, na.rm = TRUE), 'sd'= apply(PON_LU[,2:5], 1, sd, na.rm = TRUE)) 
 
 
  count<-zonal(maom_1kE*0+1 , CLC, 'sum', na.rm=T)
  countE<-apply(count[,2:5], 1, mean, na.rm = TRUE)

  fraction<-as.data.frame(cbind(POM_LUE, MAOM_LUE[,2:3], countE,  PON_LUE[,2:3], MAON_LUE[,2:3]))
  names(fraction)<-c("GRID_CODE", "POM", "POMsd", "MAOM", "MAOMsd", "area", "PON", "PONsd", "MAON", "MAONsd")


fraction<-merge(CLC12class, fraction, by.x = "GRID_CODE", by.y = "GRID_CODE")
   fraction$LC<-as.factor(c('arable', 'arable', 'permanent crops', 'permanent crops', 'permanent crops', 'grassland', 
                'heterogenous agriculture','heterogenous agriculture','heterogenous agriculture','heterogenous agriculture', 
                'broad-leaved forests', 'coniferous forests', 'mixed forests',
                'shrubland', 'shrubland', 'shrubland', 'shrubland')) 

fract_r<-aggregate(. ~ LC, fraction[,5:14], sum)
fract_r<-fract_r[c(1,2,4,6,7,9,3,5,8,10)]

   fract_r<- melt(fract_r, id=c("LC", "area"))
   value_N<-fract_r$value[17:32]
   value_sd<-fract_r$value[33:48]
   value_Nsd<-fract_r$value[49:64]

   fract_r<-fract_r[1:16,]
   fract_r$value_N<- value_N
   fract_r$value_sd<- value_sd
   fract_r$value_Nsd<- value_Nsd

   fract_r$area_n<-fract_r[,2]/max(fract_r[,2])
   fract_r$mean<-fract_r$value/fract_r$area
   fract_r$CN<-fract_r$value/fract_r$value_N







###################################################################################
#MAPS and PLOTS
###################################################################################scatter plot MAOM vs SOC

 library(gridExtra)
  library(rasterVis)
   library(quantreg)
    library(sf)
     library(grid)
 

#########################################
#FIG 1
#########################################
dem<-raster("E:\\SIT\\caprese\\soil\\dem\\hillshade1x1LAEA.tif") 
   e <- extent(sand)
   dem<-crop(dem,e)
   
mapaSHP <- rgdal::readOGR('E:\\SIT\\caprese\\adm_units\\country_def.shp')

t0<-levelplot(dem*0, par.settings = GrTheme, maxpixels=1e7)


### average Ensemble ###

maom_1kE_m<-(maom_1k + maom_1k1 + maom_1k2 + maom_1k3)/4
pom_1kE_m<- (pom_1k + pom_1k1 + pom_1k2 + pom_1k3)/4
m_t_SOCrE_m<-(m_t_SOCr + m_t_SOCr1 + m_t_SOCr2 + m_t_SOCr3)/4
BD_1kE_m<-(BD_1kE[[1]]+ BD_1kE[[2]])/2
sandE<- 100-(silt_clay+silt_clayI)/2

### average Ensemble ###

#maom_1k<-stack('MEMS\\MAOM_MEMS1k.tif')[[1]]
#pom_1k<-stack('MEMS\\POM_MEMS1k.tif')[[1]]
#m_t_SOCr<-maom_1k/(maom_1k+pom_1k)

soc_lab=expression("g C kg"^-1*"soil")

#par.settings = rasterTheme(viridis_pal(direction = -1, option = "D")(255))


FIG1.1 <- levelplot(maom_1kE_m, par.settings = viridisTheme, margin = F, maxpixels=1e7, 
scales = list(draw = FALSE), 
colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label=soc_lab, cex=0.8)) +
as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))


FIG1.2 <- levelplot(m_t_SOCrE_m, at= c(seq(0, 1, 0.1)), par.settings = RdBuTheme, margin = F, maxpixels=1e7, 
scales = list(draw = FALSE), 
colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label="MAOM:SOC", cex=0.8)) +
as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))


#### POM
S1.2 <- levelplot(pom_1kE_m, par.settings = viridisTheme, margin = F, maxpixels=1e7, 
scales = list(draw = FALSE), 
colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label=soc_lab, cex=0.8)) +
as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))

grid.arrange(S1.2, FIG1.1,  FIG1.2, ncol=3)



####supplementary SOC
SOC20c_1kE_m<-maom_1kE_m+pom_1kE_m

S1.1 <- levelplot(SOC20c_1kE_m, par.settings = magmaTheme, margin = F, maxpixels=1e7, 
scales = list(draw = FALSE), 
colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label=soc_lab, cex=0.8)) +
as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))

grid.arrange(S1.1, S1.2, ncol=2)



#########################################
#FIG SAT
#########################################


SOCf<-cbind(as.vector(SOC20c_1kE_m), as.vector(maom_1kE_m))
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



#########################################
#FIG 2
######################################### 
 
value_bar<-fract_r$value
value_bar[1:8]<-fract_r$value[1:8]+fract_r$value[9:16]
value_barr<-value_bar/fract_r$area
value_baer<-value_sd/fract_r$area



FIG2.1<-ggplot(data=fract_r, aes(x=reorder(LC,value), y=value/10000, fill=variable, width=area_n)) +
  geom_bar(stat="identity")+
   ylab(expression("SOC Tg"^'')) + xlab("") +
  geom_text(aes(label=round(value_N/10000,0)), position = position_stack(vjust = 0.5), size=2.8) +
  theme( legend.position = c(.85, .20),
        legend.title = element_blank(), 
       panel.background = element_rect(fill = "white", colour = "grey50"))+
annotate("text", x = 1, y = 5000, label = "N Tg") + scale_fill_brewer(palette='Dark2')+
geom_linerange(aes(ymin=(value_bar-value_sd)/10000, ymax=(value_bar+value_sd)/10000), alpha=0.6, size=0.5, position="identity")+
coord_flip() 

FIG2.2<-ggplot(data=fract_r, aes(x=reorder(LC,mean), y=mean, fill=variable, width=area_n)) +
   geom_bar(stat="identity") +
  ylab(expression("SOC Mg ha"^-1)) + xlab("") +
geom_text(aes(label=round(CN,0)), position = position_stack(vjust = 0.5), size=2.8) +
  theme(legend.position = "none", 
  panel.background = element_rect(fill = "white", colour = "grey50"))+
annotate("text", x = 1, y = 90, label = "C:N ratio") + scale_fill_brewer(palette='Dark2')+
geom_linerange(aes(ymin=(value_barr-value_baer), ymax=(value_barr+value_baer)), alpha=0.6, size=0.5, position="identity")+
coord_flip() 



grid.newpage()
grid.draw(cbind(ggplotGrob(FIG2.2), ggplotGrob(FIG2.1), size = "last"))
 
################################################################################# MAOM saturation
 # f <- function(x1, x2)  ifelse(x1 %in% c(18,23,24,25,26), x2, NA)
 # MAOM_sat <- overlay(CLC, maom_1k, fun=f, forcefun=TRUE)
 # 
 # MAOM_sat2<- focal(MAOM_sat, w=matrix(1,nrow=5,ncol=5), fun=max, na.rm=T) 
 # 
 # writeRaster(MAOM_sat2-maom_1k, filename="MAOM_def.tif", format="GTiff",overwrite=T)
################################################################################# 


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

#########################################
#FIG EXTENDED 1
#########################################

FIG_EXT1<-levelplot(CLCr, att='CLCnew', margin = F, maxpixels=1e7,
 scales = list(draw = FALSE), 
 col.regions=c('darkgoldenrod1','darkgoldenrod1','chocolate3','chocolate3','chocolate3',
 'chartreuse3','gold4','gold4','gold4','gold4',
 'forestgreen','darkgreen','darkolivegreen','burlywood4','burlywood4',
 'burlywood4','burlywood4')) +
as.layer(t0, under = TRUE)+ layer(sp.lines(mapaSHP, lwd=0.1, col='darkgray')) + layer(sp.points(data, lwd=0.1, col='black')) 

#writeRaster(CLCr, filename="CLCr.tif", format="GTiff",overwrite=T)


#########################################
#FIG sd MAOM and POM
#########################################

sd_maom<-calc(maom_1kE, sd, na.rm=T)
sd_pom<-calc(pom_1kE, sd, na.rm=T)

FIGsdPOM <- levelplot(sd_pom, par.settings = viridisTheme, margin = F, maxpixels=1e7, 
scales = list(draw = FALSE), 
colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label=soc_lab, cex=0.8)) +
as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))

FIGsdMAOM <- levelplot(sd_maom, par.settings = viridisTheme, margin = F, maxpixels=1e7, 
scales = list(draw = FALSE), 
colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label=soc_lab, cex=0.8)) +
as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))

grid.arrange(FIGsdPOM , FIGsdMAOM, ncol=2)

