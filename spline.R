data = read.csv("//172.27.137.244/g_boilay$/bureau/data2.csv")

#Library

library(tidyverse)
library(dplyr)
library(stringr)
library(janitor)
library(corrplot)
library(BBmisc)
library(lubridate)
library(arsenal)
library(ggfortify)
library(lmtest)
library(questionr)
library(broom)
library(broom.helpers)
library(ggeffects)
library(jtools)
library("writexl")
library(gtsummary)
library(gt)
library(ggcorrplot)
library(vcd)
library(splines)
library(npreg)

# glm

glm_reg = function(y,x,data){
  glm(formula=reformulate(termlabels = x, response= y), data=data, family = binomial(logit))
}

# lm
  
lm_reg = function(y,x,data){
  lm(formula=reformulate(termlabels = x, response= y), data=data)  
}  

 
#data_chol abeta42 continue

data_chol = filter(data, is.na(abeta42) == F,
                   is.na(age) == F,
                   is.na(sexe) == F,
                   is.na(apoe.reg2) == F,
                   is.na(mmse) == F,
                   is.na(ch) == F,
                   is.na(dia) == F,
                   is.na(hta_med) == F,
                   is.na(imc_reg) == F,
                   is.na(presence_biomarqueur) == F,
                   is.na(niveau_etude_reg) == F,
                   is.na(cholesterol_total) == F
)

data_chol$cholesterol_total = scale(data_chol$cholesterol_total)
data_chol$cholesterol_hdl = scale(data_chol$cholesterol_hdl)
data_chol$cholesterol_ldl = scale(data_chol$cholesterol_ldl)

data_chol$sexe = as.factor(data_chol$sexe)
data_chol$ch = as.factor(data_chol$ch)
data_chol$dia = as.factor(data_chol$dia)
data_chol$hta_med = as.factor(data_chol$hta_med)
data_chol$apoe.reg2 = as.factor(data_chol$apoe.reg2)
data_chol$ma = as.factor(data_chol$ma)
data_chol$ma_mci = as.factor(data_chol$ma_mci)


#write.csv(data_chol,"//172.27.137.244/g_boilay/alternance/export/data_chol.csv")


data_ma_0 = data_chol %>% filter(ma == 0)
data_ma_1 = data_chol %>% filter(ma == 1)

data_mamci_0 = data_chol %>% filter(ma_mci == 0)
data_mamci_1 = data_chol %>% filter(ma_mci == 1)

data_apoe_0 = data_chol %>% filter(apoe.reg2 == 0)
data_apoe_1 = data_chol %>% filter(apoe.reg2 != 0)

data_ma2_0 = data_chol %>% filter(ma2 == 0)
data_ma2_1 = data_chol %>% filter(ma2 == 1)
data_ma2_2 = data_chol %>% filter(ma2 == 2)


data_chol$cholesterol_total_sq = data_chol$cholesterol_total^2

# NS

m1 = c("cholesterol_total","ma")
d = data_chol
y = "abeta42"

summary(lm_reg(y,m1,d))

chol_spline <- ns(data_chol$cholesterol_total, df = 3)

plot(data_chol$abeta42 ~ chol_spline[,3])
plot(data_chol$abeta42 ~ data_chol$cholesterol_total)

plot(data_chol$cholesterol_total,data_chol$abeta42,col="grey")
lines(chol_spline,lwd=2,col="purple")



m_spline = c("chol_spline","ma")
d = data_chol
y = "abeta42"

summary(lm_reg(y,m_spline,d))

## bs

plot(data_chol$abeta42 ~ data_chol$cholesterol_total)

bs(age, knots = c(25, 40, 60))

## smooth spline

fit2<-smooth.spline(data_chol$cholesterol_total,data_chol$abeta42,cv = TRUE)
fit2

plot(data_chol$cholesterol_total,data_chol$abeta42,col="grey")
lines(fit2,lwd=2,col="purple")
legend("topright",("Smoothing Splines with 3.065585 df selected by CV"),col="purple",lwd=2)


### spline 5
fit1<-smooth.spline(data_chol$cholesterol_total,data_chol$abeta42,df = 3)
fit1

plot(data_chol$cholesterol_total,data_chol$abeta42,col="grey")
lines(fit1,lwd=2,col="purple")
legend("topright",("Smoothing Splines with 3 df"),col="purple",lwd=2)


## outlier

data_chol_outlier = filter(data, is.na(abeta42) == F,
                   is.na(age) == F,
                   is.na(sexe) == F,
                   is.na(apoe.reg2) == F,
                   is.na(mmse) == F,
                   is.na(ch) == F,
                   is.na(dia) == F,
                   is.na(hta_med) == F,
                   is.na(imc_reg) == F,
                   is.na(presence_biomarqueur) == F,
                   is.na(niveau_etude_reg) == F,
                   is.na(cholesterol_total) == F,
                   cholesterol_total < 3.5
)

data_chol_outlier$cholesterol_total = scale(data_chol_outlier$cholesterol_total)
data_chol_outlier$cholesterol_hdl = scale(data_chol_outlier$cholesterol_hdl)
data_chol_outlier$cholesterol_ldl = scale(data_chol_outlier$cholesterol_ldl)

data_ma_outlier_0 = data_chol_outlier %>% filter(ma == 0)
data_ma_outlier_1 = data_chol_outlier %>% filter(ma == 1)




## smooth spline

fit2<-smooth.spline(data_chol_outlier$cholesterol_total,data_chol_outlier$abeta42,cv = TRUE)
fit2
plot(data_chol_outlier$cholesterol_total,data_chol_outlier$abeta42,col="grey")
lines(fit2,lwd=2,col="purple")
legend("topright",("Smoothing Splines with 2.003235 df selected by CV"),col="purple",lwd=2)



### spline 5
fit1<-smooth.spline(data_chol_outlier$cholesterol_total,data_chol_outlier$abeta42,df = 4)
fit1

plot(data_chol_outlier$cholesterol_total,data_chol_outlier$abeta42,col="grey")
lines(fit1,lwd=2,col="purple")
legend("topright",("Smoothing Splines with 4 df"),col="purple",lwd=2)





fit = ss(data_chol_outlier$cholesterol_total,data_chol_outlier$abeta42, df = 2)
summary(fit)
fit$bic



plot(data_ma_outlier_0$cholesterol_total,data_ma_outlier_0$abeta42)

plot(data_ma_outlier_1$cholesterol_total,data_ma_outlier_1$abeta42)

data_ma_outlier_0$cholesterol_total_sqrt = data_ma_outlier_0$cholesterol_total^2
data_ma_outlier_1$cholesterol_total_sqrt = data_ma_outlier_1$cholesterol_total^2
data_chol_outlier$cholesterol_total_sqrt = data_chol_outlier$cholesterol_total^2


reg0 = lm(abeta42 ~ cholesterol_total + cholesterol_total_sqrt, data = data_ma_outlier_0)
summary(reg0)

reg1 = lm(abeta42 ~ cholesterol_total + cholesterol_total_sqrt, data = data_ma_outlier_1)
summary(reg1)

plot(data_ma_outlier_0$cholesterol_total,data_ma_outlier_0$abeta42)
abline(reg0)

plot(data_ma_outlier_1$cholesterol_total,data_ma_outlier_1$abeta42)
abline(reg1)

#####

data_ma_outlier_1_sp1 = data_ma_outlier_1 %>% filter(cholesterol_total <= 0)
data_ma_outlier_1_sp2 = data_ma_outlier_1 %>% filter(cholesterol_total > 0)

reg1 = lm(abeta42 ~ cholesterol_total + cholesterol_total_sqrt, data = data_ma_outlier_1_sp1)
plot(data_ma_outlier_1_sp1$cholesterol_total,data_ma_outlier_1_sp1$abeta42)
abline(reg1)

reg2 = lm(abeta42 ~ cholesterol_total + cholesterol_total_sqrt, data = data_ma_outlier_1_sp2)
plot(data_ma_outlier_1_sp2$cholesterol_total,data_ma_outlier_1_sp2$abeta42)
abline(reg2)

reg = lm(abeta42 ~ cholesterol_total + cholesterol_total_sqrt, data = data_ma_outlier_1)
plot(data_ma_outlier_1$cholesterol_total,data_ma_outlier_1$abeta42)
abline(reg)


# chol hdl / ldl

data_chol_hdl = data_chol %>% filter(is.na(cholesterol_hdl)==F)

fit2<-smooth.spline(data_chol_hdl$cholesterol_hdl,data_chol_hdl$abeta42,cv = TRUE)

plot(data_chol_hdl$cholesterol_hdl,data_chol_hdl$abeta42,col="grey")
lines(fit2,lwd=2,col="purple")
legend("topright",(""),col="purple",lwd=2)


  data_chol_ldl = data_chol %>% filter(is.na(cholesterol_ldl)==F)
  
  fit2<-smooth.spline(data_chol_ldl$cholesterol_ldl,data_chol_ldl$abeta42,cv = TRUE)
  fit2
  plot(data_chol_ldl$cholesterol_ldl,data_chol_ldl$abeta42,col="grey")
  lines(fit2,lwd=2,col="purple")
  legend("topright",("Smoothing Splines with 2.002034 df selected by CV"),col="purple",lwd=2)
  
  fit2
######################################

## outlier

data_chol_outlier = filter(data, is.na(abeta42) == F,
                   is.na(age) == F,
                   is.na(sexe) == F,
                   is.na(apoe.reg2) == F,
                   is.na(mmse) == F,
                   is.na(ch) == F,
                   is.na(dia) == F,
                   is.na(hta_med) == F,
                   is.na(imc_reg) == F,
                   is.na(presence_biomarqueur) == F,
                   is.na(niveau_etude_reg) == F,
                   is.na(cholesterol_total) == F,
                   cholesterol_total < 3.5
)

data_chol_outlier$cholesterol_total = scale(data_chol_outlier$cholesterol_total)
data_chol_outlier$cholesterol_hdl = scale(data_chol_outlier$cholesterol_hdl)
data_chol_outlier$cholesterol_ldl = scale(data_chol_outlier$cholesterol_ldl)

## ldl

data_chol_outlier_ldl = filter(data, is.na(abeta42) == F,
                           is.na(age) == F,
                           is.na(sexe) == F,
                           is.na(apoe.reg2) == F,
                           is.na(mmse) == F,
                           is.na(ch) == F,
                           is.na(dia) == F,
                           is.na(hta_med) == F,
                           is.na(imc_reg) == F,
                           is.na(presence_biomarqueur) == F,
                           is.na(niveau_etude_reg) == F,
                           is.na(cholesterol_ldl) == F
)

data_chol_outlier_ldl$cholesterol_total = scale(data_chol_outlier_ldl$cholesterol_total)
data_chol_outlier_ldl$cholesterol_hdl = scale(data_chol_outlier_ldl$cholesterol_hdl)
data_chol_outlier_ldl$cholesterol_ldl = scale(data_chol_outlier_ldl$cholesterol_ldl)


write.csv(data_chol_outlier_ldl,"C:/Users/g_boilay/Documents/data_chol_ldl.csv")

### spline 


# Données d'exemple
x <- data_chol_outlier_ldl$cholesterol_ldl
y <- data_chol_outlier_ldl$abeta42

fit1<-smooth.spline(x,y,df=6)

plot(x,y,col="grey",xlab="Age",ylab="Wages")
lines(fit1,col="red",lwd=2)

## test

fit2<-smooth.spline(x,y,cv = TRUE)
plot(x,y,col="grey")
lines(fit2,lwd=2,col="purple")
legend("topright",("Smoothing Splines with 6.78 df selected by CV"),col="purple",lwd=2)

## smooth spline

fit2<-smooth.spline(data_chol_outlier_ldl$cholesterol_ldl,data_chol_outlier_ldl$abeta42,cv = TRUE)
fit2
plot(data_chol_outlier_ldl$cholesterol_ldl,data_chol_outlier_ldl$abeta42,col="grey")
lines(fit2,lwd=2,col="purple")
legend("topright",("Smoothing Splines with 2.003235 df selected by CV"),col="purple",lwd=2)



### spline 5
fit0 = lm(abeta42 ~ cholesterol_total,data = data_chol_outlier)
summary(fit0)$r.squared

extractAIC(fit0,k = 3)
  
fit1<- ss(data_chol_outlier$cholesterol_total,data_chol_outlier$abeta42,df = 2)
summary(fit1)$adj.r.squared

fit2<- ss(data_chol_outlier$cholesterol_total,data_chol_outlier$abeta42,df = 3)
summary(fit2)$adj.r.squared

fit3<- ss(data_chol_outlier$cholesterol_total,data_chol_outlier$abeta42,df = 4)
summary(fit3)$adj.r.squared

fit5<- ss(data_chol_outlier$cholesterol_total,data_chol_outlier$abeta42,df = 6)
summary(fit5)$adj.r.squared

fit10<- ss(data_chol_outlier$cholesterol_total,data_chol_outlier$abeta42,df = 11)
summary(fit10)$adj.r.squared

  
fit <- gam(abeta42 ~ s(cholesterol_total, df = 4), data = data_chol_outlier)




# Installer et charger la bibliothèque ggplot2
install.packages("ggplot2")
library(ggplot2)

# Créer des données simulées
x <- data_ma_outlier_1$abeta42
y <- data_ma_outlier_1$cholesterol_total

# Créer un dataframe avec les données
data <- data.frame(x, y)

# Tracer le graphique avec une cassure en x = 0
ggplot(data, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x * (x >= 0), se = FALSE) +
  theme_minimal()



## chol 1 et 2

data_chol_outlier$chol1 = ifelse(data_chol_outlier$cholesterol_total <= 0, data_chol_outlier$cholesterol_total, 0)
data_chol_outlier$chol2 = ifelse(data_chol_outlier$cholesterol_total <= 0, 0, data_chol_outlier$cholesterol_total)


data_ma_0$chol1 = ifelse(data_ma_0$cholesterol_total <= 0, data_ma_0$cholesterol_total, 0)
data_ma_0$chol2 = ifelse(data_ma_0$cholesterol_total <= 0, 0, data_ma_0$cholesterol_total)

data_ma_1$chol1 = ifelse(data_ma_1$cholesterol_total <= 0, data_ma_1$cholesterol_total, 0)
data_ma_1$chol2 = ifelse(data_ma_1$cholesterol_total <= 0, 0, data_ma_1$cholesterol_total)

regall = lm(abeta42 ~ chol1 + chol2 + ma + chol1*ma + chol2*ma, data = data_chol_outlier)
summary(regall)

reg0 = lm(abeta42 ~ chol1 + chol2, data = data_ma_0)
summary(reg0)

reg1 = lm(abeta42 ~ chol1 + chol2, data = data_ma_1)
summary(reg1)

regclassic = lm(abeta42 ~ cholesterol_total + ma + cholesterol_total*ma, data = data_chol_outlier)
summary(regclassic)

###########

base = c("age","sexe","niveau_etude_reg","mmse","ch","dia","hta_med","imc_reg")

m1 = c("chol1","chol2","ma","age","sexe","niveau_etude_reg","mmse","ch","dia","hta_med","imc_reg","apoe.reg2","chol1*ma","chol2*ma")
d = data_chol_outlier
y = "abeta42"

summary(lm_reg(y,m1,d))


m2 = c("chol1","chol2","age","sexe","niveau_etude_reg","mmse","ch","dia","hta_med","imc_reg","apoe.reg2")
d = data_ma_0
y = "abeta42"

summary(lm_reg(y,m2,d))


m3 = c("chol1","chol2","age","sexe","niveau_etude_reg","mmse","ch","dia","hta_med","imc_reg","apoe.reg2")
d = data_ma_1
y = "abeta42"

summary(lm_reg(y,m3,d))