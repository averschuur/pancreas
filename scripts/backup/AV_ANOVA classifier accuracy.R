


> perf_acc1_rf <- perf_acc1[,1]
> perf_acc1_rf <- as.data.frame(perf_acc1[,1])
> View(perf_acc1_rf)
> perf_acc1_rf <- perf_acc1_rf %>%
  +     mutate(method = "RF")
> colnames(perf_acc1_rf) <- c("Accuracy", "method")
> perf_acc1_rf <- perf_acc1_rf[,c(2,1)]


perf_acc1_xgb <- as.data.frame(perf_acc1[,2])

perf_acc1_xgb <- perf_acc1_xgb %>%
     mutate(method = "xgb")
colnames(perf_acc1_xgb) <- c("Accuracy", "method")
perf_acc1_xgb <- perf_acc1_xgb[,c(2,1)]


perf_acc1_nn <- as.data.frame(perf_acc1[,2])

perf_acc1_nn <- perf_acc1_nn %>%
  mutate(method = "nn")
colnames(perf_acc1_nn) <- c("Accuracy", "method")
perf_acc1_nn <- perf_acc1_nn[,c(2,1)]

perf_acc1_all <- rbind(perf_acc1_rf, perf_acc1_xgb, perf_acc1_nn)

boxplot(Accuracy~method , data=perf_acc1_all)

pg<-lm(Accuracy~method , data=perf_acc1_all)
anova(pg)

