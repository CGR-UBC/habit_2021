
sessionDesc <- function(data, var, sessionVar) {
  all <- c((mean(data[,var]) + sd(data[,var])), mean(data[,var]), (mean(data[,var]) - sd(data[,var])))
  means <- tapply(data[,var], data[,sessionVar], mean)
  sds <- tapply(data[,var], data[,sessionVar], sd)
  means_plus <- means + sds
  means_minus <- means - sds
  labels <- c("Session 1", "\n\nSession 2", "\n\nSession 3")
  result <- paste(labels, "\n+1SD:", means_plus, "\nmean:", means, "\n-1SD:", means_minus)
  all_labels <- c("+1SD:", "\nMean", "\n-1SD:")
  all_result <- paste(all_labels, all)
  cat("Overall\n",all_result, "\n\n", result)
  all_hist <- ggplot(data, aes_string(x = var)) +
    ggtitle("All Sessions") +
    geom_histogram(fill = "white", colour = "black", binwidth = 1) 
  session_hist <- ggplot(data, aes_string(x = var)) +
    ggtitle("Per Session") +
    geom_histogram(fill = "white", colour = "black") +
    facet_grid(sessionVar)
  list(all_hist, session_hist)
}

term_replacer <- function(orig, new, model) {
  formula <- (gsub(orig, new, formula(model)))
  formula <- paste(formula[2], formula[1], formula[3])
  call <- model["call"]
  call[["call"]]["formula"] <- formula
  eval(parse(text = call))
}

diagnostics <- function(model, wide = FALSE) {
  name <- deparse(substitute(model))
  if (!file.exists(file.path("diagnostics", "OLS", name))) {dir.create(file.path("diagnostics", "OLS", name), recursive = TRUE)}
  jpeg_width <- ifelse(wide, 1200, 600)

  resid <- sjPlot::plot_model(model, type = "resid", show.data = TRUE)
  sjPlot::save_plot(file.path("diagnostics", "OLS", name, paste(name, "_resid.jpg")), fig = resid, width = 40, height = 20)
  
  if (!(model$family["family"] == "binomial")) {
    slope <- sjPlot::plot_model(model, type = "slope", show.data = TRUE)
    sjPlot::save_plot(file.path("diagnostics", "OLS", name, paste(name, "_slope.jpg")), fig = slope, width = 40, height = 20)
  }
  
  jpeg(file.path("diagnostics", "OLS", name, paste(name, "_performance.jpg")), height = 1500, width = jpeg_width)
  print(performance::check_model(model, panel = TRUE, verbose = FALSE, 
                                 check = c("qq", "normality", "linearity", "ncv", "homogeneity")))
  dev.off()
  
  outliers <- plot(performance::check_outliers(model))
  ggplot2::ggsave(filename = file.path("diagnostics", "OLS", name, paste(name, "_outliers.jpg")),
                  plot = outliers, dpi = 600, width = 10, height = 8)
}

# robust plots
# plot showing reduction in deviance in robust model
robust_diagnostics<-function(model,robustmodel){
  name <- deparse(substitute(robustmodel))
  if (!file.exists(file.path("diagnostics", "robust", name))) {dir.create(file.path("diagnostics", "robust", name), recursive = TRUE)}
  
  dev1<-abs(robustmodel$residuals)
  dev2<-abs(model$residuals)
  
  n <- length(dev1)
  ord1 <- order(dev1)
  sdev1 <- sort(dev1) 
  sdev2 <- sort(dev2) 
  
  jpeg(file.path("diagnostics", "robust", name, 'robust_residuals.jpg'),width = 1200, height = 500, units = "px",pointsize = 20)
  
  par(mfrow=c(1,3)) 
  plot(ppoints(n), sdev1, type="b",pch=1,xlab="quantiles", ylab= " residuals")
  lines(ppoints(n), sdev2, type="b",pch=2)
  xuu <- ppoints(n)[n]
  text(xuu - .03, max(sdev1) + .1, ord1[n])
  text(xuu, max(sdev2) + .3, ord1[n])
  legend(x="topleft",legend=c("Robust","Conventional"), pch=c(1,2))
  plot(sdev1,ylab= " Robust model residuals")
  plot(sdev2,ylab= " Conventional model residuals")
  dev.off()
  
  jpeg(file.path("diagnostics", "robust", name, 'robust_qq.jpg'),width = 1200, height = 500, units = "px",pointsize = 20)
  par(mfrow=c(1,2)) 
  
  qqnorm(robustmodel$residuals, ylab="Robust model Residuals",  xlab="Normal Scores") 
  qqline(robustmodel$residuals)
  qqnorm(model$residuals, ylab="Conventional model residuals",  xlab="Normal Scores") 
  qqline(model$residuals)
  dev.off()
  
}

robmodel_weights<-function(data4diagnostics,model4diagnostics,DV){
  name <- deparse(substitute(model4diagnostics))
  if (!file.exists(file.path("diagnostics", "robust", name))) {dir.create(file.path("diagnostics", "robust", name), recursive = TRUE)}
  
  if (DV=='ISI'){data4diagnostics$rweights<-model4diagnostics$rweights
  } else {data4diagnostics$rweights<-model4diagnostics$w.r
  }
  
  jpeg(file.path("diagnostics", "robust", name, 'weights_by_participant.jpg'),width = 1200, height = 500, units = "px",pointsize = 20)
  par(mfrow=c(1,1)) 
  cat("Participant weights               ", plot(data4diagnostics$Participant.ID,data4diagnostics$rweights),"\n\n")
  dev.off()
  
  data4diagnostics$rweights_binary<-cut(data4diagnostics$rweights, breaks = c(-Inf, 0.5, Inf), 
                                        labels = c("< 0.75","> 0.75"))
  weights_binary_per_p<- data4diagnostics %>% group_by(Participant.f,rweights_binary) %>% tally() %>%  ungroup() %>%
    tidyr::complete(Participant.f,rweights_binary, fill = list(N = 0, freq = 0))
  weight_trials_per_p<- data4diagnostics %>% group_by(Participant.f) %>% tally() %>%  ungroup() %>%
    tidyr::complete(Participant.f, fill = list(N = 0, freq = 0))
  small_weights<-subset(weights_binary_per_p,rweights_binary=='< 0.75')
  large_weights<-subset(weights_binary_per_p,rweights_binary=='> 0.75')
  weight_trials_per_p$p_small<-small_weights$n/(small_weights$n+large_weights$n)
  write.table(weight_trials_per_p,file.path("diagnostics", "robust", name, "weight_less_than_point75_trials__proportion_per_p.csv"),sep=",",row.names=FALSE)
  
  
  if(DV=='BET'){
    weights_binary_per_p_by_BET<- data4diagnostics %>% group_by(Participant.f,rweights_binary,Next.Bet.Changed) %>% tally() %>%  ungroup() %>%
      tidyr::complete(Participant.f,rweights_binary,Next.Bet.Changed, fill = list(N = 0, freq = 0))
    write.table(weights_binary_per_p_by_BET, file.path("diagnostics", "robust", name, "by_BET_weight_less_than_point75_trials__proportion_per_p.csv"),sep=",",row.names=FALSE)
    next_bet_1<-subset(data4diagnostics,Next.Bet.Changed== 1)
    next_bet_0<-subset(data4diagnostics,Next.Bet.Changed == 0)
    cat("Percent weights small next BET 0            ", nrow(subset(next_bet_0,rweights_binary=='< 0.75'))/nrow(next_bet_0),"\n")
    cat("Percent weights small next BET 1            ", nrow(subset(next_bet_1,rweights_binary=='< 0.75'))/nrow(next_bet_1),"\n\n")
    
  }
  
  modelVars <- all.vars(formula(model4diagnostics))[3:length(all.vars(formula(model4diagnostics)))]
  rows=ceiling(length(modelVars)/2)
  
  data4diagnostics_small<-subset(data4diagnostics,rweights_binary=='< 0.75')
  data4diagnostics_large<-subset(data4diagnostics,rweights_binary=='> 0.75')
  jpeg(file.path("diagnostics", "robust", name, "robust_weights.jpg"),width = 2400, height = 2400, units = "px")
  par(mfrow=c(rows,2)) 
  for (var in modelVars) {
    plot(data4diagnostics[,var], data4diagnostics$rweights , ylab="weights", xlab=var,col = rgb(red = 0, green = 0, blue = 0, alpha = 0.3))
    if (!is.factor(data4diagnostics[,var])){
      cat(var," weights small (median)            ", mean(data4diagnostics_small[,var]),"\n")
      cat(var," weights large (median)            ", mean(data4diagnostics_large[,var]),"\n")
      cat(var," weights correlation               ", cor(data4diagnostics[,var],data4diagnostics$rweights),"\n\n")
    }
  }
  dev.off()
  
  if (DV == "ISI"){
    cat(var," weights correlation               ", cor(data4diagnostics$Time.Until.Next.Bet.log,data4diagnostics$rweights),"\n\n")
    jpeg(file.path("diagnostics", "robust", name, "DV_and_weights.jpg"),width = 600, height = 600, units = "px")
    plot(data4diagnostics$Time.Until.Next.Bet.log,data4diagnostics$rweights,col = rgb(red = 0, green = 0, blue = 0, alpha = 0.3))
    abline(lm(data4diagnostics$Time.Until.Next.Bet.log~data4diagnostics$rweights), col="red") 
    dev.off()
  } 
  if (DV == "BET"){
    cat(DV," weights correlation               ", cor(data4diagnostics$Next.Bet.Changed,data4diagnostics$rweights),"\n\n")
    jpeg(file.path("diagnostics", "robust", name, "DV_and_weights.jpg"),width = 600, height = 600, units = "px")
    plot(data4diagnostics$Next.Bet.Changed,data4diagnostics$rweights,col = rgb(red = 0, green = 0, blue = 0, alpha = 0.3))
    abline(lm(data4diagnostics$Next.Bet.Changed~data4diagnostics$rweights), col="red") 
    dev.off()
  }
  #return(as.data.frame(data4diagnostics)$rweights)
}

tab_formatter <- function(model, preds, sf = 2) {
  table <- tail(broom::tidy(model, conf.int = TRUE), length(preds))[,c(1,2,6,7,5)]
  table[1] <- preds
  colnames(table) <- c("Predictors", "Estimate", "95% CI (lower)", "95% CI (upper)", "p")
  
  if (!is.null(model$family) && model$family$family == "binomial") {
    colnames(table)[2] <- "Odds Ratio"
    for (i in 2:4) {table[i] <- exp(table[i])}
  }
  
  for (i in 2:4) {
    if (data.class(table[,i]) == "numeric") {
      table[,i] <- ifelse((table[,i] < 1), 
                          formatC(signif(table[,i], digits = sf), digits = sf, format = "fg", flag = "#"), 
                          formatC(round(table[,i], digits = sf), digits = sf + 1, flag = "#"))
    }
  }  
  
  table[,"p"] <- ifelse(table[,"p"] < 0.001, "< .001", 
                        formatC(round(table[,"p"], digits = 3), digits = 3, format = "f", flag = "#"))
  
  table[, "95% CI"] <- paste(table[,3], table[,4], sep = ", ")
  table <- table[,c(1,2,6,5)]
  table <- rbind(table, list("Observations", paste(nrow(model$model)),"","",""))
  
  return(table)
}