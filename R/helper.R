## Helper
get_mapping_to_full_dataset = function(X, W, assignment) {

  mapping = vector("numeric", length = dim(X)[1])
  count = 1
  for (i in 1:dim(X)[1]) {
    if (W[i] == assignment) {
      mapping[i] = count
      count = count + 1
    }
  }
  return(mapping)
}

get_oob_predictions = function(X, forest, mapping) {

  raw_preds = predict(forest, X, predict.all = TRUE)$predictions
  final_preds = vector("numeric", length = dim(X)[1])
  inbag_counts = forest$inbag.counts

  for (i in 1:dim(X)[1]) {
    if (mapping[i] == 0 || i > length(mapping)) {
      final_preds[i] = mean(raw_preds[i,])
    } else {
      temp_preds = vector("list", length = forest$num.trees)
      for (j in 1:forest$num.trees) {
        if (inbag_counts[j][[1]][mapping[i]] == 0) {
          temp_preds[[j]] = raw_preds[i,j]
        }
      }
      final_preds[i] = mean(unlist(Filter(is.numeric, temp_preds)))
    }
  }
  return(final_preds)
}

predict_x_learner <- function(X, W, estimate_propensities, predict_oob) {
  if (predict_oob) {
    preds_1 = get_oob_predictions(X, train_x_learner_fit$xf1, train_x_learner_fit$mapping1)
    preds_0 = get_oob_predictions(X, train_x_learner_fit$xf0, train_x_learner_fit$mapping0)
  } else {
    preds_1 = predict(train_x_learner_fit$xf1, X)$predictions
    preds_0 = predict(train_x_learner_fit$xf0, X)$predictions
  }

  if (estimate_propensities) {
    propf = ranger(W ~ ., data = data.frame(X, W = W),
                   min.node.size = 1)
    ehat = propf$predictions
    preds = (1 - ehat) * preds_1 + ehat * preds_0
  } else {
    preds = 0.5 * preds_1 + 0.5 * preds_0
  }
  return(preds=preds)
}


train_x_learner <- function(X, W, Y, ipcw, save = TRUE, fileList = c("xf0", "xf1", "mapping0", "mapping1", "yhat0", "yhat1")) {
  if (all(file.exists(paste0("./model/", fileList, ".rds")))) {
    xf0 <- readRDS("./model/xf0.rds")
    xf1 <- readRDS("./model/xf1.rds")
    mapping0 <- readRDS("./model/mapping0.rds")
    mapping1 <- readRDS("./model/mapping1.rds")
    yhat0 <- readRDS("./model/yhat0.rds")
    yhat1 <- readRDS("./model/yhat1.rds")
    print("Load models")
  } else {
    hyper_grid <- expand.grid(
      mtry       = 1:4,
      node_size  = 1:3,
      num.trees = seq(50,500,50),
      OOB_RMSE   = 0
    )

    for(i in 1:nrow(hyper_grid)) {
      # train model
      rf <- ranger(
        formula        = Y ~ .,
        data           = data.frame(X[W == 0,], Y = Y[W == 0]),
        num.trees      = hyper_grid$num.trees[i],
        mtry           = hyper_grid$mtry[i],
        min.node.size  = hyper_grid$node_size[i],
        case.weights   = rep(ipcw, nrow(X[W == 0,])),
        importance = 'impurity')
      # add OOB error to grid
      hyper_grid$OOB_RMSE[i] <- sqrt(rf$prediction.error)
    }
    position = which.min(hyper_grid$OOB_RMSE)
    # fit best model
    tf0 <- ranger(Y ~ .,data = data.frame(X[W == 0,], Y = Y[W == 0]),
                  num.trees = hyper_grid$num.trees[position],
                  importance = 'impurity',
                  probability = FALSE,
                  min.node.size = hyper_grid$node_size[position],
                  mtry = hyper_grid$mtry[position],
                  case.weights = rep(ipcw, nrow(X[W == 0,])),
                  keep.inbag = T)
    # tf0 = ranger(Y ~ ., data = data.frame(X[W == 0,], Y = Y[W == 0]),
    #              num.trees = 1000, min.node.size = 30,#43
    #              case.weights = rep(ipcw, nrow(X[W == 0,])))
    yhat0 = predict(tf0, X[W == 1,])$predictions

    for(i in 1:nrow(hyper_grid)) {
      rf <- ranger(
        formula        = Y ~ .,
        data           = data.frame(Y = Y[W == 1] - yhat0, X[W == 1,]),
        num.trees      = hyper_grid$num.trees[i],
        mtry           = hyper_grid$mtry[i],
        min.node.size  = hyper_grid$node_size[i],
        case.weights   = rep(ipcw, nrow(X[W == 1,])),
        importance = 'impurity')
      # add OOB error to grid
      hyper_grid$OOB_RMSE[i] <- sqrt(rf$prediction.error)
    }
    position = which.min(hyper_grid$OOB_RMSE)

    xf1 <- ranger(Y ~ .,data = data.frame(Y = Y[W == 1] - yhat0, X[W == 1,]),
                  num.trees = hyper_grid$num.trees[position],
                  importance = 'impurity',
                  probability = FALSE,
                  min.node.size = hyper_grid$node_size[position],
                  mtry = hyper_grid$mtry[position],
                  case.weights = rep(ipcw, nrow(X[W == 1,])),
                  keep.inbag = T)
    # xf1 = ranger(Y ~ ., data = data.frame(Y = Y[W == 1] - yhat0, X[W == 1,]),
    #              keep.inbag = TRUE, num.trees = 500, min.node.size = 24,
    #              case.weights = rep(ipcw, nrow(X[W == 1,])), importance = "impurity")
    mapping1 = get_mapping_to_full_dataset(X, W, 1)


    for(i in 1:nrow(hyper_grid)) {
      rf <- ranger(
        formula        = Y ~ .,
        data           = data.frame(X[W == 1,], Y = Y[W == 1]),
        num.trees      = hyper_grid$num.trees[i],
        mtry           = hyper_grid$mtry[i],
        min.node.size  = hyper_grid$node_size[i],
        case.weights   = rep(ipcw, nrow(X[W == 1,])),
        importance = 'impurity')
      # add OOB error to grid
      hyper_grid$OOB_RMSE[i] <- sqrt(rf$prediction.error)
    }
    position = which.min(hyper_grid$OOB_RMSE)

    tf1 <- ranger(Y ~ .,data = data.frame(X[W == 1,], Y = Y[W == 1]),
                  num.trees = hyper_grid$num.trees[position],
                  importance = 'impurity',
                  probability = FALSE,
                  min.node.size = hyper_grid$node_size[position],
                  mtry = hyper_grid$mtry[position],
                  case.weights = rep(ipcw, nrow(X[W == 1,])),
                  keep.inbag = T)
    # tf1 = ranger(Y ~ ., data = data.frame(X[W == 1,], Y = Y[W == 1]),
    #              num.trees = 500, min.node.size = 24,
    #              case.weights = rep(ipcw, nrow(X[W == 1,])))
    yhat1 = predict(tf1, X[W == 0,])$predictions

    for(i in 1:nrow(hyper_grid)) {
      rf <- ranger(
        formula        = Y ~ .,
        data           = data.frame(Y = yhat1 - Y[W == 0], X[W == 0,]),
        num.trees      = hyper_grid$num.trees[i],
        mtry           = hyper_grid$mtry[i],
        min.node.size  = hyper_grid$node_size[i],
        case.weights   = rep(ipcw, nrow(X[W == 0,])),
        importance = 'impurity')
      # add OOB error to grid
      hyper_grid$OOB_RMSE[i] <- sqrt(rf$prediction.error)
    }
    position = which.min(hyper_grid$OOB_RMSE)

    xf0 <- ranger(Y ~ .,data.frame(Y = yhat1 - Y[W == 0], X[W == 0,]),
                  num.trees = hyper_grid$num.trees[position],
                  importance = 'impurity',
                  probability = FALSE,
                  min.node.size = hyper_grid$node_size[position],
                  mtry = hyper_grid$mtry[position],
                  case.weights = rep(ipcw, nrow(X[W == 0,])),
                  keep.inbag = T)

    # xf0 = ranger(Y ~ ., data= data.frame(Y = yhat1 - Y[W == 0], X[W == 0,]),
    #              keep.inbag = TRUE, num.trees = 1000, min.node.size = 30,
    #              case.weights = rep(ipcw, nrow(X[W == 0,])), importance = "impurity")
    mapping0 = get_mapping_to_full_dataset(X, W, 0)

    if (save) {
      saveRDS(xf0, "./model/xf0.rds")
      saveRDS(xf1, "./model/xf1.rds")
      saveRDS(mapping0, "./model/mapping0.rds")
      saveRDS(mapping1, "./model/mapping1.rds")
      saveRDS(yhat0, "./model/yhat0.rds")
      saveRDS(yhat1, "./model/yhat1.rds")
    }
  }

  return(list(xf0 = xf0, xf1 = xf1, mapping0 = mapping0, mapping1 = mapping1, yhat0=yhat0, yhat1=yhat1))
}



#
# #######
# train_x_learner <- function(X, W, Y, ipcw, save = TRUE, fileList = c("xf0", "xf1", "mapping0", "mapping1", "yhat0", "yhat1")) {
#   if (all(file.exists(paste0("./model/", fileList, ".rds")))) {
#     xf0 <- readRDS("./model/xf0.rds")
#     xf1 <- readRDS("./model/xf1.rds")
#     mapping0 <- readRDS("./model/mapping0.rds")
#     mapping1 <- readRDS("./model/mapping1.rds")
#     yhat0 <- readRDS("./model/yhat0.rds")
#     yhat1 <- readRDS("./model/yhat1.rds")
#     print("Load models")
#   } else {
#
#     tf0 = ranger(Y ~ ., data = data.frame(X[W == 0,], Y = Y[W == 0]),
#                 num.trees = 2000, min.node.size = 1,
#                 case.weights = rep(ipcw, nrow(X[W == 0,])))
#     yhat0 = predict(tf0, X[W == 1,])$predictions
#
#     xf1 = ranger(Y ~ ., data = data.frame(Y = Y[W == 1] - yhat0, X[W == 1,]),
#                 keep.inbag = TRUE, num.trees = 1000, min.node.size = 1,
#                 case.weights = rep(ipcw, nrow(X[W == 1,])), importance = "impurity")
#     mapping1 = get_mapping_to_full_dataset(X, W, 1)
#
#     tf1 = ranger(Y ~ ., data = data.frame(X[W == 1,], Y = Y[W == 1]),
#                  num.trees = 2000, min.node.size = 1,
#                  case.weights = rep(ipcw, nrow(X[W == 0,])))
#     yhat1 = predict(tf1, X[W == 0,])$predictions
#     xf0 = ranger(Y ~ ., data= data.frame(Y = yhat1 - Y[W == 0], X[W == 0,]),
#                 keep.inbag = TRUE, num.trees = 1000, min.node.size = 30,
#                 case.weights = rep(ipcw, nrow(X[W == 0,])), importance = "impurity")
#     mapping0 = get_mapping_to_full_dataset(X, W, 0)
#
#     if (save) {
#       saveRDS(xf0, "./model/xf0.rds")
#       saveRDS(xf1, "./model/xf1.rds")
#       saveRDS(mapping0, "./model/mapping0.rds")
#       saveRDS(mapping1, "./model/mapping1.rds")
#       saveRDS(yhat0, "./model/yhat0.rds")
#       saveRDS(yhat1, "./model/yhat1.rds")
#     }
#   }
#
#   return(list(xf0 = xf0, xf1 = xf1, mapping0 = mapping0, mapping1 = mapping1, yhat0=yhat0, yhat1=yhat1))
# }

#######
Calib_risk_x <- function(data, n.tile){

  m_ps = glm(as.factor(test$Y)~predict_x_learner_fit2+W, family=binomial(), data=data)
  prs_df = data.frame(pr_score = predict(m_ps, type = 'response'))

  test2<-cbind(data, prs_df)

  test.ex <- mutate(test2, quantile_rank = ntile(test2$pr_score, n.tile))

  pred_prr = data.frame()
  obs_prr = data.frame()


  for (i in (1:n.tile)) {

    pdat = test.ex[test.ex$quantile_rank==i,]

    pred_prr[i,1] =mean(pdat$pr_score)

    #Observed ARR for Propensity score

    obs_prr[i,1] = mean(pdat$Y)

  }
  paired_rr = cbind(pred_prr, obs_prr)
  names(paired_rr) = c("pred_prr", "obs_prr")

  return(paired_rr)
}

# calibration plot(for Risk-Probabbility with PRECISE-DAPT model)
Calib_risk_dapt <- function(data, n.tile){

  m_ps = glm(as.factor(dapt.test$Y)~precise_DAPT, family=binomial(), data=dapt.test)
  prs_df = data.frame(pr_score = predict(m_ps, type = 'response'))

  dapt.test2<-cbind(dapt.test, prs_df)

  test.ex <- mutate(dapt.test2, quantile_rank = ntile(dapt.test2$pr_score, n.tile))

  pred_prr = data.frame()
  obs_prr = data.frame()

  for (i in (1:n.tile)) {

    pdat = test.ex[test.ex$quantile_rank==i,]

    pred_prr[i,1] =mean(pdat$pr_score)

    #Observed ARR for Propensity score
    m_ps2 = glm(as.factor(Y)~as.factor(W), family=binomial(), data=pdat)
    summary(m_ps2)

    prs_df2 = data.frame(y_hat = predict(m_ps2, type = 'response'), W=pdat$W)

    #obs_prr[i,1] = mean(prs_df2$y_hat)
    obs_prr[i,1] = mean(pdat$Y)
  }
  paired_rr = cbind(pred_prr, obs_prr)
  names(paired_rr) = c("pred_prr", "obs_prr")

  return(paired_rr)
}


####### Survival plot
surv_func = function(test){

  test$DAPT<-ifelse(test$W==1,"Short", "Long")
  test$DAPT2<-as.factor(test$DAPT)
  test$Time<-ifelse(test$DAYS.TIMI.MAJOR.BLEEDING<365, test$DAYS.TIMI.MAJOR.BLEEDING, 365)
  test$class<-ifelse(test$predict_x_learner_fit<=0, "benefit", "no-benefit")
  test$class2<-ifelse(test$predict_x_learner_fit<=0, 1, 0)

  fit1<-survfit(Surv(Time, Y) ~ class, data = test)
  fit2<-survfit(Surv(Time, Y) ~ DAPT, data = subset(test, class=="benefit"))
  fit3<-survfit(Surv(Time, Y) ~ DAPT, data = subset(test, class=="no-benefit"))

  summary(fit1);summary(fit2);summary(fit3)

  return(list(test=test, fit1=fit1, fit2=fit2, fit3=fit3, test=test))

}


####### C-index

c_indx<-function(data, gp){

  Surv.obj<-Surv(data[,c('Time')], data[,c('Y')])
  cox.obj<-coxph(Surv.obj ~data[,c(gp)], data=data)
  summary(cox.obj)
  tidy(cox.obj)

  data.frame(tidy(cox.obj, exponentiate = T, conf.int = T)[, c("term", "estimate", "conf.low", "conf.high", "p.value")])
  conco_ben=concordance(cox.obj)$concordance


  test.benefit<-subset(data, get(gp)=="benefit")
  Surv.obj2<-Surv(test.benefit$Time, test.benefit$Y)
  cox.obj2<-coxph(Surv.obj2 ~ DAPT , data=test.benefit)
  summary(cox.obj2);summary(cox.obj2)
  tidy(cox.obj2)

  data.frame(tidy(cox.obj2, exponentiate = T, conf.int = T)[, c("term", "estimate", "conf.low", "conf.high", "p.value")])
  conco_dapt=concordance(cox.obj2)$concordance

  return(list(conco_ben=conco_ben, conco_dapt=conco_dapt))

}



######### Aply
aply=function(name) {

  var<-df[,c(name)]

  df<-RCT %>%dplyr::select(one_of(all_variables_names))

  df_ex<-c()
  df_Missi_ex<-c()
  df_train_ex<-c()
  df_test_ex<-c()


  df_ex <- df %>% rename(Y=name, W=DAPT.SHORT.LONG)

  df_ex$W<-ifelse(df_ex$W==1,1,0)
  df_ex$SEX2<-ifelse(df_ex$SEX==1,0,1)
  df_ex$PRESENT_MI<-ifelse(df_ex$CLINICAL.PRESENTATION == 3,1,0)
  df_ex<-df_ex[,!names(df) %in% c("W")]

  # Converting all columns to numerical
  df_ex <- data.frame(lapply(df_ex, function(x) as.numeric(as.character(x))))

  df_Missi_ex<-subset(df_ex,is.na(df_ex$WBC.COUNT)==FALSE & is.na(df_ex$BASELINE.HEMOGLOBIN)==FALSE)

  n <- dim(df_Missi_ex)[1]

  set.seed(1)

  df_train_ex <- df_Missi_ex[df_Missi_ex$STUDY!=3,]
  df_test_ex <- df_Missi_ex[df_Missi_ex$STUDY==3,]#TICO



  class_ex<-surv_d$test[,c("class")]
  test_ex<-cbind(df_test_ex, class_ex)

  test_ex$DAPT<-ifelse(test_ex$W==1,"Short", "Long")
  test_ex$DAPT2<-as.factor(test_ex$DAPT)

  Time<- paste0("DAYS.", paste(name))

  test_ex$Time<-ifelse(test_ex[,c(Time)]<360, test_ex[,c(Time)], 360)
  Time<-ifelse(test_ex[,c(Time)]<360, test_ex[,c(Time)], 360)

  fit1<-survfit(Surv(Time, Y) ~ class_ex, data = test_ex)
  fit3<-survfit(Surv(Time, Y) ~ DAPT, data = subset(test_ex, class_ex=="benefit"))
  fit4<-survfit(Surv(Time, Y) ~ DAPT, data = subset(test_ex, class_ex=="no-benefit"))

  return(list(fit1=fit1, fit3=fit3, fit4=fit4, test_ex=test_ex, name=name))
}


# calibration plot(for Risk-Probabbility with PRECISE-DAPT model)
Calib_risk_dapt <- function(data, n.tile){

  m_ps = glm(as.factor(dapt.test$Y)~precise_DAPT+W, family=binomial(), data=dapt.test)
  prs_df = data.frame(pr_score = predict(m_ps, type = 'response'))

  dapt.test2<-cbind(dapt.test, prs_df)

  test.ex <- mutate(dapt.test2, quantile_rank = ntile(dapt.test2$pr_score, n.tile))

  pred_prr = data.frame()
  obs_prr = data.frame()

  for (i in (1:n.tile)) {

    pdat = test.ex[test.ex$quantile_rank==i,]

    pred_prr[i,1] =mean(pdat$pr_score)

    #Observed ARR for Propensity score
    m_ps2 = glm(as.factor(Y)~as.factor(W), family=binomial(), data=pdat)
    summary(m_ps2)

    prs_df2 = data.frame(y_hat = predict(m_ps2, type = 'response'), W=pdat$W)

    #obs_prr[i,1] = mean(prs_df2$y_hat)
    obs_prr[i,1] = mean(pdat$Y)
  }
  paired_rr = cbind(pred_prr, obs_prr)
  names(paired_rr) = c("pred_prr", "obs_prr")

  return(paired_rr)
}


theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))

}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

cal_plot <- function(model, model_name, pred_var, ...){

  require(tidyverse)
  require(viridis)
  require(gridExtra)

  # The calibration plot
  g1 <- mutate(data, bin = ntile(get(pred_var), 5)) %>%
    # Bin prediction into 10ths
    group_by(bin) %>%
    mutate(n = n(), # Get ests and CIs
           bin_pred = mean(get(pred_var)),
           bin_prob = mean(as.numeric(Y) - 1),
           se = sqrt((bin_prob * (1 - bin_prob)) / n),
           ul = bin_prob + 1.96 * se,
           ll = bin_prob - 1.96 * se) %>%
    ungroup() %>%
    ggplot(aes(x = bin_pred, y = bin_prob, ymin = ll, ymax = ul)) +
    geom_pointrange(size = 0.5, color = "black") +
    scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by = 0.01)) +
    scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by = 0.01)) +
    geom_abline() + # 45 degree line indicating perfect calibration
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed",
                color = "black", formula = y~-1 + x) +
    # straight line fit through estimates
    # geom_smooth(aes(x = get(pred_var), y = as.numeric(Y) - 1),
    #             color = "red", se = FALSE, method = "loess") +
    # loess fit through estimates
    xlab("") +
    ylab("Observed Probability") +
    theme_minimal() +
    ggtitle(model_name)

  # The distribution plot
  g2 <- ggplot(data, aes(x = get(pred_var))) +
    geom_histogram(fill = "black", bins = 200) +
    scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by = 0.01)) +
    xlab("Predicted Probability") +
    ylab("") +
    theme_minimal() +
    #scale_y_continuous(breaks = c(0, 40)) +
    theme(panel.grid.minor = element_blank())

  # Combine them
  g <- arrangeGrob(g1, g2, respect = TRUE, heights = c(1, 0.25), ncol = 1)
  grid.newpage()
  grid.draw(g)
  return(g[[3]])

}
