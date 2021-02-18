############################################################################### 
#dSOC PROJECTIONS
############################################################################### 
#########################################
#FIG EXTENDED 1
#########################################

LCS<-subset(LCS, LU %in% c("SHR", "GR", "FR", "CR"))
LCS$LU<- droplevels(LCS$LU)


lm1<-ggplot(data=LCS, aes(MAT, OC_pom_g_kg, colour = RAIN)) + geom_point(alpha=0.5, size=1.8) + ylab("POM")+ 
facet_wrap(~LU)+
geom_smooth(method='lm',formula= y ~ x)+
scale_colour_gradientn(colours = terrain.colors(10))


lm2<-ggplot(data=LCS, aes(MAT, OC_sc_g_kg, colour = RAIN)) + geom_point(alpha=0.5, size=1.8) + ylab("MAOM")+ 
facet_wrap(~LU)+
geom_smooth(method='lm',formula= y ~ x)+
scale_colour_gradientn(colours = terrain.colors(10))



GCM<-c('HADes', 'CNRM', 'IPLS')

for (i in 1:3)  {


###GCMprj
WCLpr<-stack(paste0("E:\\SIT\\caprese\\meteo\\worldclim2\\",GCM[i],"8.5_2080_LAEA.tif"))

MAT1kpr<-resample(WCLpr[[1]], silt_clay)  
RAIN1kpr<-resample(WCLpr[[2]], silt_clay)

dT<- MAT1kpr-MAT1k
dP<- RAIN1kpr-RAIN1k

assign(paste0('dT', '_', GCM[i]), dT)
assign(paste0('dP', '_', GCM[i]), dP)


dC_LU<-data.frame('LU'= c('CR', 'FR', 'GR', 'SH', 'CR', 'FR', 'GR', 'SH'), 
'FRAC'= c('POM', 'POM', 'POM', 'POM', 'MAOM', 'MAOM', 'MAOM', 'MAOM'), 'dC'=rep(0, times = 8), 
'area'=rep(0, times = 8), 'rel_change'=rep(0, times = 8))

RAIN_SAND<-dP*sandE


#POM:j=0 and MAOM:j=4
for(j in c(0,4)){

if (j==0){
###POM
 md0<-lm(LCS$OC_pom_g_kg~(LCS$MAT+LCS$RAIN*LCS$sand)*LCS$LU)
#md1<-lm(LCS$OC_pom_g_kg~(LCS$MAT+LCS$RAIN)*LCS$LU)

STcr<-md0$coefficients[2]
STfr<-md0$coefficients[2]+md0$coefficients[9]
STgr<-md0$coefficients[2]+md0$coefficients[10]
STsh<-md0$coefficients[2]+md0$coefficients[11]

SPcr<-md0$coefficients[3]
SPfr<-md0$coefficients[3]+md0$coefficients[12]
SPgr<-md0$coefficients[3]+md0$coefficients[13]
SPsh<-md0$coefficients[3]+md0$coefficients[14]

TXcr<-md0$coefficients[8]
TXfr<-md0$coefficients[8]+md0$coefficients[18]
TXgr<-md0$coefficients[8]+md0$coefficients[19]
TXsh<-md0$coefficients[8]+md0$coefficients[20]


}else{
###MAOM
 md0<-lm(LCS$OC_sc_g_kg~(LCS$MAT+LCS$RAIN*LCS$sand)*LCS$LU)
#md1<-lm(LCS$OC_sc_g_kg~(LCS$MAT+LCS$RAIN)*LCS$LU)

STcr<-md0$coefficients[2]
STfr<-md0$coefficients[2]+md0$coefficients[9]
STgr<-md0$coefficients[2]+md0$coefficients[10]
STsh<-md0$coefficients[2]+md0$coefficients[11]

SPcr<-md0$coefficients[3]
SPfr<-md0$coefficients[3]+md0$coefficients[12]
SPgr<-md0$coefficients[3]+md0$coefficients[13]
SPsh<-md0$coefficients[3]+md0$coefficients[14]

TXcr<-md0$coefficients[8]
TXfr<-md0$coefficients[8]+md0$coefficients[18]
TXgr<-md0$coefficients[8]+md0$coefficients[19]
TXsh<-md0$coefficients[8]+md0$coefficients[20]

}



f <- function(x1, x2, x3, x4)  ifelse(x1 %in% c(12,13,15,16,17,19,20,21,22), x2*STcr + x3*SPcr + x4*TXcr, NA)		#CR
 CRdt <- overlay(CLC, dT, dP, RAIN_SAND, fun=f, forcefun=TRUE)
 CRdt <- CRdt*(maom_1kE_m*0+1)

f <- function(x1, x2, x3, x4)  ifelse(x1 %in% c(23,24,25), x2*STfr + x3*SPfr + x4*TXfr, NA)		#FR
 FRdt <- overlay(CLC, dT, dP, RAIN_SAND, fun=f, forcefun=TRUE)
 FRdt <- FRdt*(maom_1kE_m*0+1)
  
f <- function(x1, x2, x3, x4)  ifelse(x1 %in% c(18), x2*STgr + x3*SPgr + x4*TXgr, NA)		#GR
 GRdt <- overlay(CLC, dT, dP, RAIN_SAND, fun=f, forcefun=TRUE)
 GRdt <- GRdt*(maom_1kE_m*0+1)

f <- function(x1, x2, x3, x4)  ifelse(x1 %in% c(26,27,28,29), x2*STsh + x3*SPsh + x4*TXsh, NA)		#SH
 SHdt <- overlay(CLC, dT, dP, RAIN_SAND, fun=f, forcefun=TRUE)
 SHdt <- SHdt*(maom_1kE_m*0+1)


dMT<-merge(CRdt,FRdt,GRdt,SHdt)
dMT<-dMT*BD_1kE_m*2				#t/ha of C
assign(paste0('dMT', j, '_', GCM[i]), dMT)


cellStats(dMT, sum, na.mr=T)*100/1000000	
tar<-cellStats(dMT*0+1, sum, na.mr=T)

dC_LU[(1+j),3]<-cellStats(CRdt*BD_1kE_m*2, sum, na.mr=T)*100/1000000	
dC_LU[(2+j),3]<-cellStats(FRdt*BD_1kE_m*2, sum, na.mr=T)*100/1000000
dC_LU[(3+j),3]<-cellStats(GRdt*BD_1kE_m*2, sum, na.mr=T)*100/1000000
dC_LU[(4+j),3]<-cellStats(SHdt*BD_1kE_m*2, sum, na.mr=T)*100/1000000

dC_LU[(1+j),4]<-cellStats(CRdt*0+1, sum, na.mr=T)/tar	
dC_LU[(2+j),4]<-cellStats(FRdt*0+1, sum, na.mr=T)/tar
dC_LU[(3+j),4]<-cellStats(GRdt*0+1, sum, na.mr=T)/tar
dC_LU[(4+j),4]<-cellStats(SHdt*0+1, sum, na.mr=T)/tar
 
}


dC_LU$rel_change[1]<-dC_LU$dC[1]/((fract_r$value[1]+fract_r$value[5]+fract_r$value[7])/10000)*100
dC_LU$rel_change[2]<-dC_LU$dC[2]/((fract_r$value[2]+fract_r$value[3]+fract_r$value[6])/10000)*100
dC_LU$rel_change[3]<-dC_LU$dC[3]/((fract_r$value[4])/10000)*100
dC_LU$rel_change[4]<-dC_LU$dC[4]/((fract_r$value[8])/10000)*100

dC_LU$rel_change[5]<-dC_LU$dC[5]/((fract_r$value[9]+fract_r$value[13]+fract_r$value[15])/10000)*100
dC_LU$rel_change[6]<-dC_LU$dC[6]/((fract_r$value[10]+fract_r$value[11]+fract_r$value[14])/10000)*100
dC_LU$rel_change[7]<-dC_LU$dC[7]/((fract_r$value[12])/10000)*100
dC_LU$rel_change[8]<-dC_LU$dC[8]/((fract_r$value[16])/10000)*100

##
assign(paste0('dC_LU_', GCM[i]), dC_LU)

}


dC_LUavg<- (dC_LU_CNRM[,3:5]+  dC_LU_HADes[,3:5]+ dC_LU_IPLS[,3:5])/3
dC_LUavg<-cbind(dC_LU_HADes[,1:2], dC_LUavg, 'sd'=apply(cbind(dC_LU_CNRM[3], dC_LU_HADes[3], dC_LU_IPLS[3]),1,sd))
dC_LUavg<-cbind(dC_LUavg, 'corr'=dC_LUavg$dC)
dC_LUavg$corr[5:8]<-ifelse((dC_LUavg$dC[1:4]<0 & dC_LUavg$dC[5:8]<0), dC_LUavg$dC[1:4]+dC_LUavg$dC[5:8], dC_LUavg$dC[5:8])


FIG3.2<-ggplot(data=dC_LUavg, aes(x=reorder(LU,-dC), y=dC, fill=FRAC)) +
  geom_bar(stat="identity", position="stack")+
 ylab("Tg of C") + xlab("") + ggtitle("\n \n Cumulative changes ")+
theme(legend.position = c(.25, .10),
        legend.title = element_blank(), 
       panel.background = element_rect(fill = "white", colour = "grey50"))+
geom_hline(yintercept=0, linetype="dashed", colour = "grey50")+
#geom_text(aes(y=-1700, label=round(rel_change,1)), size=3, position = position_dodge(0.3) )+
scale_fill_brewer(palette='Dark2', direction=-1)+
geom_linerange(aes(ymin=corr-sd, ymax=corr+sd), alpha=0.6, size=0.5)+
coord_flip()


dMTpom<-(dMT0_CNRM + dMT0_HADes+ dMT0_IPLS)/3
dMTmaom<-(dMT4_CNRM + dMT4_HADes+ dMT4_IPLS)/3

#-20, 20, 5
FIG3.1 <- levelplot(dMTpom, at= c(seq(-25, 25, 5)), par.settings = RdBuTheme, margin = F, maxpixels=1e7, 
scales = list(draw = FALSE), 
colorkey=list(title = "",  height=0.8, width=0.9, space="top"), main=list(label=expression(Delta~"POM (Mg C ha"^-1*")"), cex=0.8)) +
as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))



FIG3.3 <- levelplot(dMTmaom, at= c(seq(-25, 25, 5)), par.settings = RdBuTheme, margin = F, maxpixels=1e7, 
scales = list(draw = FALSE), 
colorkey=list(title = "",  height=0.8, width=0.9, space="top"), main=list(label=expression(Delta~"MAOM (Mg C ha"^-1*")"), cex=0.8))+
as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))

grid.arrange(FIG3.1, FIG3.2, FIG3.3, ncol=3)



############################################################################### 

#########################################
#FIG EXTENDED 1
#########################################
dTavg<- dT-(dT_CNRM + dT_HADes+ dT_IPLS)/3
dPavg<- dP-(dP_CNRM + dP_HADes+ dP_IPLS)/3

q<-quantile(dTavg, probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), type=7, names = T)

figS_dT<- levelplot(dTavg, at= q, par.settings = viridisTheme, margin = F, maxpixels=1e7, 
scales = list(draw = FALSE), 
colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label=expression(Delta~"T (°C)"), cex=0.8))+
as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))



figS_dP<- levelplot(dPavg, at= c(seq(-250, 250, 25)), par.settings = RdBuTheme, margin = F, maxpixels=1e7, 
scales = list(draw = FALSE), 
colorkey=list(title = "", height=0.8, width=0.9, space="top"), main=list(label=expression(Delta~"P (mm)"), cex=0.8)) + 
as.layer(t0, under = TRUE) + layer(sp.lines(mapaSHP, lwd=0.2, col='darkgray'))


grid.arrange(figS_dT, figS_dP, ncol=2)




####
####supplementary LCS texture distribution 
####

LCSall<-rbind(LCS, LCSpred)

##ver.3.5
tex <- data.frame(
"CLAY" = LCSall$clay,
"SILT" = 100-LCSall$clay-LCSall$sand,
"SAND"=LCSall$sand
"LU"=LCSall$LU
)


library(ggtern)
ggtern(data=tex,aes(SILT,CLAY, SAND)) + geom_point(alpha=0.2)+
geom_point(data = tex[1:352,], col = 'red', size=1.5)

tex<-subset(tex, LU %in% c("SHR", "GR", "FR", "CR"))


ggtern(data=tex,aes(SILT,CLAY, SAND)) + 
stat_density_tern(geom='polygon',
                         n         = 400,
                         aes(fill  = ..level..,
                             alpha = ..level..)) +
geom_point(alpha=0) +facet_wrap(~LU)+
guides(color = "none", fill = "none", alpha = "none")


 hist(LCSall$OC,  main='', xlab="SOC (g/ka)", freq=FALSE, lty=1)
 hist(LCS$OC_tf,  add=T, col=alpha("red", 0.2), freq=FALSE,  lty=0 )



ggplot(data=LCS, aes(pH_in_H2O, OC_pom_g_kg, colour = OC_tf)) + geom_point(alpha=0.5, size=1.8) + 
ylab("POM (g/kg)")+ 
facet_wrap(~LU)+
geom_smooth(method = lm, formula = y ~x)+
scale_colour_gradientn(colours = terrain.colors(10))


############################################################################### 
#dSOC PROJECTIONS mrt!!!!!!!!!!!!!!!!!!!!!!!!
###############################################################################

MRT<-LCS$SOC_lcs/(LCS$NPP)

FIGmrt<-ggplot(data=LCS, aes(MAT, MRT)) + geom_point(alpha=0.5, size=1.8) + 
ylab("MRT")+ ylim(0, 40)+
facet_wrap(~LU)+
geom_smooth(method = lm, formula = y ~x)+
scale_colour_gradientn(colours = terrain.colors(10)) +
scale_y_continuous(sec.axis = sec_axis(~ . /40))+
geom_point(aes(MAT, (OC_sc_g_kg /OC_tf)*40), col='red', alpha=0.5)+
geom_smooth(aes(MAT, (OC_sc_g_kg /OC_tf)*40), method=lm,colour="red")





