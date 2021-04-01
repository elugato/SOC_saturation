# 01/04/2021 E. Lugato - SOC_saturation"

R script '1_LCS_RF_sc.v2.r'
* **Random Forest (RF) model predicting LUCAS measured MAOM fraction** 

`RF<-randomForest(OC_sc_g_kg ~ s_c_prc + pH_in_H2O + OC_tf + Ndep_WD_tx + MAT + EROS + WT, data = LCS, ntree=1000, mtry=4, importance=TRUE, na.action=na.omit)`

---

R script '1_LCS_RF_sc.v2.1.r'
* **RF model prediction on 1km spatial layers**

```
# upscaling 1km

'inputR<-addLayer(inputR, pH_H2O, SOC20c_1k, NDEP1k, MAT1k, EROS1k, WT)'
'names(inputR)<-c('s_c_prc', 'pH_in_H2O',  'OC_tf',  'Ndep_WD_tx',  'MAT', 'EROS', 'WT')'
  'maom_1k<- predict(inputR, RF, type='response', progress='window', na.rm=T)'
  
```

---

R script '1_LCS_RF_sc.v2.2___.r'
* **sensitivity to MAT and P**
```
# model fitting and upscaling with WordClim projections

POM  'md0<-lm(LCS$OC_pom_g_kg~(LCS$MAT+LCS$RAIN*LCS$sand)*LCS$LU)'
MAOM 'md0<-lm(LCS$OC_maom_g_kg~(LCS$MAT+LCS$RAIN*LCS$sand)*LCS$LU)'

```
