library(tidyverse)
library(lubridate)

# データ読み込み・前処理（週次の集計）
data <- data.table::fread("pcr_positive_daily.csv")

data$DATE <- ymd(data$日付)
data$Num_Wks <- week(data$DATE)

Wks <- data %>%
  dplyr::select(c("日付", "Num_Wks")) %>% 
  dplyr::distinct(Num_Wks, .keep_all = TRUE)

data_wks <- data %>% 
  dplyr::group_by(`Num_Wks`) %>% 
  summarise(`Covid_positive` = sum(`PCR 検査陽性者数(単日)`)) %>% 
  dplyr::ungroup() %>% 
  dplyr::inner_join(Wks, by="Num_Wks") 

data_wks$DoW <- ymd(data_wks$日付) 
data_wks <- data_wks %>% 
  dplyr::arrange(DoW)

## 週次推移のプロット
options(scipen=1000)
plot(data_wks$Covid_positive, xaxt="n", xlab="週", ylab="PCR検査数", type = "l", main="COVID-19感染者数推移（週次）")
axis(side=2, at=c(0, 10000, 100000, 200000, 300000, 400000))
axis(side=1,at=c(1,8, 15, 22, 29, 36, 43),labels=c(data_wks$DoW[1], data_wks$DoW[8], data_wks$DoW[15], data_wks$DoW[22], data_wks$DoW[29], data_wks$DoW[36], data_wks$DoW[43]))

# モデリング
library(KFAS)

label_list <- c()
## ローカルレベルモデル+固定季節変動
mod1 <- SSModel(data_wks$Covid_positive ? SSMtrend(1, Q=NA) + SSMseasonal(12, Q=0), H=NA)
fit1 <- fitSSM(mod1, numeric(2), method="SANN")
kfs1 <- KFS(fit1$model, smoothing=c("state", "mean", "disturbance"))
logLik1 <- kfs1$logLik - sum(kfs1$Finf>0) * log(2*pi)/2
AIC1 <- -2*logLik1 + 2*(2+12)
#label_list <- "ローカル_固定"

## ローカルレベルモデル+可変季節変動
mod2 <- SSModel(data_wks$Covid_positive ? SSMtrend(1, Q=NA) + SSMseasonal(12, Q=NA), H=NA)
fit2 <- fitSSM(mod2, numeric(3), method="SANN")
kfs2 <- KFS(fit2$model, smoothing=c("state", "mean", "disturbance"))
logLik2 <- kfs2$logLik - sum(kfs2$Finf>0) * log(2*pi)/2
AIC2 <- -2*logLik2 + 2*(3+12)
#label_list <- c(label_list, "ローカル_可変")

## 平滑化トレンドモデル+固定季節変動
mod3 <- SSModel(data_wks$Covid_positive ? SSMtrend(2, Q=list(0, NA)) + SSMseasonal(12, Q=0), H=NA)
fit3 <- fitSSM(mod3, numeric(2), method="SANN")
kfs3 <- KFS(fit3$model, smoothing=c("state", "mean", "disturbance"))
logLik3 <- kfs3$logLik - sum(kfs3$Finf>0) * log(2*pi)/2
AIC3 <- -2*logLik3 + 2*(2+13)
#label_list <- c(label_list, "平滑_固定")

## 平滑化トレンドモデル+可変季節変動
mod4 <- SSModel(data_wks$Covid_positive ? SSMtrend(2, Q=list(0, NA)) + SSMseasonal(12, Q=NA), H=NA)
fit4 <- fitSSM(mod4, numeric(3), method="SANN")
kfs4 <- KFS(fit4$model, smoothing=c("state", "mean", "disturbance"))
logLik4 <- kfs4$logLik - sum(kfs4$Finf>0) * log(2*pi)/2
AIC4 <- -2*logLik4 + 2*(3+13)
#label_list <- c(label_list, "平滑_可変")

## ローカル線形トレンドモデル+固定季節変動
mod5 <- SSModel(data_wks$Covid_positive ? SSMtrend(2, Q=list(NA, NA)) + SSMseasonal(12, Q=0), H=NA)
fit5 <- fitSSM(mod5, numeric(3), method="SANN")
kfs5 <- KFS(fit5$model, smoothing=c("state", "mean", "disturbance"))
logLik5 <- kfs5$logLik - sum(kfs5$Finf>0) * log(2*pi)/2
AIC5 <- -2*logLik5 + 2*(3+13)
#label_list <- c(label_list, "ローカル線形_固定")

## ローカル線形トレンドモデル+可変季節変動
mod6 <- SSModel(data_wks$Covid_positive ? SSMtrend(2, Q=list(NA, NA)) + SSMseasonal(12, Q=NA), H=NA)
fit6 <- fitSSM(mod6, numeric(4), method="SANN")
kfs6 <- KFS(fit6$model, smoothing=c("state", "mean", "disturbance"))
logLik6 <- kfs6$logLik - sum(kfs6$Finf>0) * log(2*pi)/2
AIC6 <- -2*logLik6 + 2*(4+13)
#label_list <- c(label_list, "ローカル線形_可変")

AIC_list <- c(round(AIC1, digits = 2), round(AIC2, digits = 2), round(AIC3, digits = 2), round(AIC4, digits = 2), round(AIC5, digits = 2), round(AIC6, digits = 2))
loglik_list <- c(round(logLik1, digits = 2), round(logLik2, digits = 2), round(logLik3, digits = 2), round(logLik4, digits = 2), round(logLik5, digits = 2), round(logLik6, digits = 2))
AIC_df <- data.frame(AIC_list)
loglik_df <- data.frame(loglik_list)
df <- AIC_df %>% 
  cbind(loglik_df)
colnames(df) <- c("AIC", "MLL")

