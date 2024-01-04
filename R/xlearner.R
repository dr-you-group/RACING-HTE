####### X-leaner
library(ranger)
library(missRanger)

fitXlearner <- train_x_learner(train[,c(covariates_category, covariates_numeric)], as.numeric(train$W), as.numeric(train$outcome), ipcw=0.5)
predTrain <- predict_x_learner(train[,c(covariates_category, covariates_numeric)], as.numeric(train$W), FALSE, TRUE)
predTest <- predict_x_learner(test[,c(covariates_category, covariates_numeric)], as.numeric(test$W), FALSE, TRUE)

plot(predTrain)
plot(predTest)

############ Survival analysis - bleeding
library(survminer)
library(survival)

survdata <- cbind(test, predTest)
survdata$outcome <- as.numeric(survdata$outcome)-1
survdata$benefit <- ifelse(survdata$predTest<=0, 1, 0)

# t <- survdata %>% filter(benefit == 1)
form <- c()
survplot_benefit = c()
dat<-c()
form[[1]] <- formula("Surv(outcomedays, outcome) ~ benefit")
dat[[1]] <- survdata

#Log-Lank Test
km_diff<-c()
km_pval<-c()
tmp<-c()
cox_out<-c()
cox<-c()

for (i in 1 : 1){
  km_diff[[i]] <- survdiff(form[[i]], data =dat[[i]])
  km_pval[i] <- 1 - pchisq(km_diff[[i]]$chisq, length(km_diff[[i]]$n) - 1)
  cox_out[[i]] <- coxph(form[[i]], data =dat[[i]]);
  cox[[i]]<-as.matrix(summary(eval(cox_out[[i]]))$conf.int)
  options(digits=3)
  print(cox_out[[i]])
  i=i+1
}

  title.nm<-c("Primary endpoint")
  y.max<-c(0.1)

  legend.lab<-c()
  legend.lab[[1]]<-c("Non-benefit", "Benefit")

  palette<-c()
  palette[[1]]<-ggsci::pal_lancet("lanonc", alpha = .6)(2)

  tt <- fmsb::riskdifference(table(dat[[1]][dat[[1]][,"W"]==0,"outcome"]>0)[2], table(dat[[1]][dat[[1]][,"W"]==1,"outcome"]>0)[2], table(dat[[1]]$W)[1], table(dat[[1]]$W)[2], CRC=FALSE, conf.level=0.95)

  plotdf <- data.frame(mean = tt$estimate, lower = tt$conf.int[1], upper = tt$conf.int[2], pval = tt$p.value)
  plotdf$mean_p <- paste(round(100*plotdf$mean, 2), "%", sep="")
  plotdf$upper_p <- paste(round(100*plotdf$upper, 2), "%", sep="")
  plotdf$lower_p <- paste(round(100*plotdf$lower, 2), "%", sep="")
  plotdf$pval_p <- round(plotdf$pval, 2)


  #plotting
  survplot_benefit[[1]]<-ggsurvplot(survfit(form[[1]], data =dat[[1]]),
                                    conf.int = F,
                                    pval = F,
                                    censor=F,
                                    break.time.by=90,
                                    break.y.by=0.01,
                                    pval.method=F,
                                    risk.table.col = "strata",
                                    ggtheme=theme_classic2(base_size=12), #, base_family="Arial"),
                                    palette = palette[[1]],
                                    fun = "cumhaz",
                                    xlim = c(0, 1095), ylim = c(0, y.max[1]),
                                    risk.table = TRUE,
                                    #tables.theme = theme_cleantable(),
                                    title = title.nm[1],
                                    fontsize=5.0,
                                    risk.table.pos="out",
                                    risk.table.y.text.col=FALSE,
                                    legend.labs = legend.lab[[1]],
                                    legend.title=""
  )

  survplot_benefit[[1]]$plot <- survplot_benefit[[1]]$plot +
    ggplot2::annotate(
      "text",
      x=180, y=0.095,
      vjust = 1, hjust = 0.5,
      label = paste0("ARD: ",plotdf$mean_p, " (95% CI ", plotdf$lower_p, " to ", plotdf$upper_p, ")"))

  survplot_benefit[[1]]$table <- ggpubr::ggpar(survplot_benefit[[1]]$table,
                                               font.title = list(size = 11, color = "black")
  )

  plotName <- paste0("preliminarySurvivalBenefitVsNonbenefit", i, ".tiff")
  tiff(plotName, width = 8, height = 8, units = 'in', res = 100)
  print(survplot_benefit[[1]])
  dev.off()

  ###########################
  ############ Survival analysis - In benefit group
  library(survminer)
  library(survival)

  survdata_benefit <- survdata[survdata$predTest<=0,]
  survdata_benefit$outcome <- as.numeric(survdata_benefit$outcome)

  survplot_inbenefit = c()
  dat<-c()

  form[[1]] <- formula("Surv(outcomedays, outcome) ~ W")
  dat[[1]] <- survdata_benefit

  #Log-Lank Test
  km_diff<-c()
  km_pval<-c()
  tmp<-c()
  cox_out<-c()
  cox<-c()

  for (i in 1 : 1){
    km_diff[[i]] <-survdiff(form[[i]], data =dat[[i]])
    km_pval[i] <- 1 - pchisq(km_diff[[i]]$chisq, length(km_diff[[i]]$n) - 1)
    cox_out[[i]] <- coxph(form[[i]], data =dat[[i]]);
    cox[[i]]<-as.matrix(summary(eval(cox_out[[i]]))$conf.int)
    options(digits=3)
    print(cox_out[[i]])
    i=i+1
  }

  # km_pval_all<-data.frame(rbind(km_pval[[1]], km_pval[[2]], km_pval[[3]]))
  # names(km_pval_all)<-c("bind.km_pval..1....km_pval..2....km_pval..3...."="log-lank")
  # cox_all<-data.frame(rbind(cox[[1]], cox[[2]], cox[[3]]))
  # Result<-cbind(km_pval_all, cox_all)
  # Result

  for (i in 1:1) {
    title.nm<-c("Primary endpoint")
    y.max<-c(0.15)

    legend.lab<-c()
    legend.lab[[1]] <- c("Monotherapy", "Combination therapy")

    palette<-c()
    palette[[1]]<-ggsci::pal_lancet("lanonc", alpha = .6)(2)


    tt <- fmsb::riskdifference(table(dat[[i]][dat[[i]][,"W"]==0,"outcome"]>0)[2], table(dat[[i]][dat[[i]][,"W"]==1,"outcome"]>0)[2], table(dat[[i]]$W)[1], table(dat[[i]]$W)[2], CRC=FALSE, conf.level=0.95)

    fmsb::riskdifference(table(dat[[1]][dat[[1]][,"W"]==0,"outcome"]>0)[2], table(dat[[1]][dat[[1]][,"W"]==1,"outcome"]>0)[2], table(dat[[1]]$W)[1], table(dat[[1]]$W)[2], CRC=FALSE, conf.level=0.95)

    plotdf <- data.frame(mean = tt$estimate, lower = tt$conf.int[1], upper = tt$conf.int[2], pval = tt$p.value)
    plotdf$mean_p <- paste(round(100*plotdf$mean, 2), "%", sep="")
    plotdf$upper_p <- paste(round(100*plotdf$upper, 2), "%", sep="")
    plotdf$lower_p <- paste(round(100*plotdf$lower, 2), "%", sep="")
    plotdf$pval_p <- round(plotdf$pval, 2)



    #plotting
    survplot_inbenefit[[i]]<-ggsurvplot(survfit(form[[i]], data =dat[[i]]),
                                        conf.int = F,
                                        pval = F,
                                        censor=F,
                                        break.time.by=90,
                                        break.y.by=0.01,
                                        pval.method=F,
                                        risk.table.col = "strata",
                                        ggtheme=theme_classic2(base_size=12), #, base_family="Arial"),
                                        palette = palette[[i]],
                                        fun = "cumhaz",
                                        xlim = c(0, 1095), ylim = c(0, y.max[i]),
                                        risk.table = TRUE,
                                        #tables.theme = theme_cleantable(),
                                        title = title.nm[i],
                                        fontsize=5.0,
                                        risk.table.pos="out",
                                        risk.table.y.text.col=FALSE,
                                        legend.labs = legend.lab[[i]],
                                        legend.title=""
    )

    survplot_inbenefit[[i]]$plot <- survplot_inbenefit[[i]]$plot +
      ggplot2::annotate(
        "text",
        x=180, y=0.095,
        vjust = 1, hjust = 0.5,
        label = paste0("ARD: ",plotdf$mean_p, " (95% CI ", plotdf$lower_p, " to ", plotdf$upper_p, ")"))

    survplot_inbenefit[[i]]$table <- ggpubr::ggpar(survplot_inbenefit[[i]]$table,
                                                   font.title = list(size = 11, color = "black")
    )

    # survplot_inbenefit[[i]]$table <- survplot_inbenefit[[i]]$table +
    #   theme(
    #     plot.title = element_text(hjust = 0, vjust = 0),
    #     plot.margin = unit(c(1, 5.5, 5.5, -8), "points"))
  }


  plotName <- paste0("preliminarySurvivalInNonbenefit", i, ".tiff")
  tiff(plotName, width = 8, height = 8, units = 'in', res = 100)
  print(survplot_inbenefit[[1]])
  dev.off()


  ####### benefit vs nonbenefit table
  survdata$benefit <- as.factor(survdata$benefit)
  survdata$outcome <- as.factor(survdata$outcome)
  table1Features <- c(covariates_category, covariates_numeric, c("outcome", "outcomedays"))
  tabUnmatched <- CreateTableOne(vars = table1Features, strata = "benefit", data = survdata, test = FALSE)
  table1 <- print(tabUnmatched, smd = TRUE)
  write.csv(table1, "nonbenefitVsBenefit.csv")

  survdata_benefit$outcome <- as.factor(survdata_benefit$outcome)
  tabUnmatched <- CreateTableOne(vars = table1Features, strata = "W", data = survdata_benefit, test = FALSE)
  table1 <- print(tabUnmatched, smd = TRUE)



library(kernelshap)
library(shapviz)

#kernelshap_xf0 <- kernelshap(train_x_learner_fit$xf0, X=survdata[,c(covariates_category, covariates_numeric)][survdata$W==0,], bg_X = survdata[,c(covariates_category, covariates_numeric)][survdata$W==0,])
#kernelshap_xf1 <- kernelshap(train_x_learner_fit$xf1, X=survdata[,c(covariates_category, covariates_numeric)][survdata$W==1,], bg_X = survdata[,c(covariates_category, covariates_numeric)][survdata$W==1,])

sv_kernelshap_xf0 <- shapviz(kernelshap_xf0)
sv_kernelshap_xf1 <- shapviz(kernelshap_xf1)

sv_importance(sv_kernelshap_xf0, kind = "bee", max_display = 30)
sv_importance(sv_kernelshap_xf1, kind = "bee", max_display = 30)

plotName <- paste0("SHAP_xf0", ".tiff")
tiff(plotName, width = 8, height = 12, units = 'in', res = 100)
t <- sv_importance(sv_kernelshap_xf0, kind = "bee", max_display = 30)
print(t)
dev.off()

plotName <- paste0("SHAP_xf1", ".tiff")
tiff(plotName, width = 8, height = 12, units = 'in', res = 100)
t <- sv_importance(sv_kernelshap_xf1, kind = "bee", max_display = 30)
print(t)
dev.off()


