# 26/05/2020 E. Lugato - SOC_saturation"

R script 
* **Random Forest (RF) model on LUCAS data fractions** 

`LCS_frac<-read.csv("E:\\SIT\\fra\\db\\fractions\\forR_allLC_CaCO3.csv", header = TRUE)`



* **RF model prediction on 1km spatial layers**

```
# upscaling 1km

'inputR<-addLayer(inputR, pH_H2O, SOC20c_1k, NDEP1k, MAT1k, EROS1k, K)'
'names(inputR)<-c('s_c_prc', 'pH_in_H2O',  'OC_tf',  'Ndep_WD_tx',  'MAT', 'EROS', 'K')'
  'maom_1k<- predict(inputR, RF, type='response', progress='window', na.rm=T)'
  
```


* **economic values to MAOM and POM**
```
# cur_p = current price 1t CO2eq; p_mrt = ratio mean residence time MAOM/POM

total euro  'tot_e<-cur_p*3.67*(sum(fraction$MAOM)+sum(fraction$POM))/10000000 '
pom price 'pom_p<-tot_e/((sum(fraction$MAOM)*p_mrt+ sum(fraction$POM))*3.67/10000000)'

```
