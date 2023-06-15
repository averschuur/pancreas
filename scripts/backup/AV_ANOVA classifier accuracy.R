

a<-conf_mat[[1]][["byClass"]]
a<-a[,11]
method = "NN"
a <- as.data.frame(cbind(a,method))
rownames(a) = NULL
colnames(a) <- c("accuracy", "method")

b<-conf_mat[[2]][["byClass"]]
b<-b[,11]
method = "RF"
b <- as.data.frame(cbind(b,method))
rownames(b) = NULL
colnames(b) <- c("accuracy", "method")

c<-conf_mat[[3]][["byClass"]]
c<-c[,11]
method = "XGB"
c <- as.data.frame(cbind(c,method))
rownames(c) = NULL
colnames(c) <- c("accuracy", "method")

all <- rbind(a,b)
all <- rbind(all,c)


attach(all)
boxplot(accuracy~method)


pg<-lm(accuracy~method)
anova(pg)
