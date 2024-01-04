library(dplyr)
library(tableone)
library(data.table)
library(moonBook)
n_seed=2024
set.seed(2024)

##data clean
setwd("~/git/RACING-HTE")
dataPath <- file.path(getwd(), "data", "RACING.xlsx")
data <- readxl::read_xlsx(dataPath)
dataPath <- file.path(getwd(), "data", "RACING_lab.xlsx")
dataLab <- readxl::read_xlsx(dataPath)

data <- as.data.frame(data)
dataLab <- as.data.frame(dataLab)

data <- dataLab
colnames(data)
# 1-2: RCT informations
# 42: intervention 1:combi / 2: Mono ****
# 3-12: Muscle events
# 13-19: Other events
# 20: PP?
# 21-39: FU informations
# 40-41 / 43-45: Demographics ****
# 46-59: comorbidities / base statin ****
# 60: PP?
# 61-105: Outcome informations
# 106-109: Baseline lab ****
# 110-113: 1y lab
# 114-117: 2y lab
# 118-121: 3y lab
# 122-125: Reduced LDLC
# 127-134: L55-L70 informations
# 135-145: baseline medications ****
# 146-156: baseline-2m medications?
# 157-168: 2m-6m medications?
# 169-180: 6m-1y medications?
# 181-192: 1y-2y medications?
# 193-204: 2-3y medications?
# 205-215: Baseline lab ****
# 216-226: 2m lab
# 227-237: 6m lab
# 238-253: 1y lab
# 254-266: 2y lab
# 267-282: 3y lab

# 1-2: RCT informations
# 42: intervention 1:combi / 2: Mono ****
# 40-41 / 43-45: Demographics ****
# 46-59: comorbidities / base statin ****
# 106-109: Baseline lab ****
# 135-145: baseline medications ****
# 205-215: Baseline lab ****

df <- data[,c(1:2, 42, 40:41, 43:45, 46:59, 106:109, 135:145, 206:215)]
features <- colnames(data[,c(1:2, 42, 40:41, 43:45, 46:59, 106:109, 135:145, 206:215)])
df <- df[,-c(22, 38, 47)]
colnames(df) <- c("number", "completed", "intervention", "age", "sex", "height", "weight", "bmi", "mi", "pci", "cabg", "cva", "ckd", "dialysis",
"pad", "hpt", "dm", "insulin", "antidiabetics", "smoking", "currentsmoking", "ldl", "total", "tg", "hdl", "aspirin", "clo", "tica",
"dpat", "noac", "bb", "acei", "arb", "ccb", "diuretics", "aldactone", "ast", "alt", "ck", "hb1ac", "hba1c6.5", "urate", "glucose", "crp")
View(df)

df$ck[df$ck == "ND"] <- NA
df$crp[df$crp == "ND"] <- NA
df$urate[df$urate == "ND"] <- NA
df$ast[df$ast == "ND"] <- NA
df$alt[df$alt == "ND"] <- NA
df$hb1ac[df$hb1ac == "NA"] <- NA
df$hba1c6.5[df$hba1c6.5 == "NA"] <- NA
df$glucose[df$glucose == "NA"] <- NA
df$crp[df$crp == "NA"] <- NA
df$W <- as.numeric(df$intervention)-1
df$W <- ifelse(df$W==1,0,1) # Change to: 0: mono, 1: combi
df <- subset(df, select=-intervention)

intervention_factor <- c("W")
covariates_category <- c("sex", "mi", "pci", "cabg", "cva", "ckd", "dialysis",
                       "pad", "hpt", "dm", "insulin", "antidiabetics", "smoking", "currentsmoking", "aspirin", "clo", "tica",
                       "dpat", "noac", "bb", "acei", "arb", "ccb", "diuretics", "aldactone", "hba1c6.5")
covariates_numeric <- c("age", "height", "weight", "bmi", "ldl", "total", "tg", "hdl", "ast", "alt", "ck", "hb1ac", "urate", "glucose", "crp")


df <- df %>% mutate_at(intervention_factor, as.factor)
df <- df %>% mutate_at(covariates_category, as.factor)
df <- df %>% mutate_at(covariates_numeric, as.numeric)

length(df) # 44
length(intervention_factor) + length(covariates_category) + length(covariates_numeric) # 42

#missRanger
df_imputed <- missRanger::missRanger(subset(df, select=-W))
#
# Missing value imputation by random forests
#
# Variables to impute:		ast, alt, ck, hb1ac, hba1c6.5, urate, glucose, crp
# Variables used to impute:	number, completed, age, sex, height, weight, bmi, mi,
# pci, cabg, cva, ckd, dialysis, pad, hpt, dm, insulin, antidiabetics, smoking,
# currentsmoking, ldl, total, tg, hdl, aspirin, clo, tica, dpat, noac, bb, acei,
# arb, ccb, diuretics, aldactone, ast, alt, ck, hb1ac, hba1c6.5, urate, glucose, crp

#Intervention
df_imputed$W <- as.numeric(df$W)-1

#Primary endpoint
df_imputed$outcome <- as.factor(dataLab[,c("Primary_endpoint")])
df_imputed$outcomedays <- as.numeric(dataLab[,c("PRIMARY_DAYS")])

## 75% of the sample size
smp_size <- floor(0.67 * nrow(df_imputed))
## set the seed to make your partition reproducible
set.seed(2024)
train_ind <- sample(seq_len(nrow(df_imputed)), size = smp_size)

train <- df_imputed[train_ind, ]
test <- df_imputed[-train_ind, ]

table1Features <- c(covariates_category, covariates_numeric, c("outcome", "outcomedays"))
tabUnmatched <- CreateTableOne(vars = table1Features, strata = "W", data = train, test = FALSE)
table1 <- print(tabUnmatched, smd = TRUE)

table1Features <- c(covariates_category, covariates_numeric, c("outcome", "outcomedays"))
tabUnmatched <- CreateTableOne(vars = table1Features, strata = "W", data = test, test = FALSE)
table1 <- print(tabUnmatched, smd = TRUE)


t <- cbind(df, primaryEndpoint = dataLab$Primary_endpoint)
t$primaryEndpoint <- as.factor(t$primaryEndpoint)
table1Features <- c(c("primaryEndpoint"), covariates_category, covariates_numeric)
tabUnmatched <- CreateTableOne(vars = table1Features, strata = "W", data = t, test = FALSE)
table1 <- print(tabUnmatched, smd = TRUE)
write.csv(table1, "baseline.csv")
