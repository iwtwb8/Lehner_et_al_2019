
library(ggplot2)
df_old = data.frame(log2_fold_change=c(-0.489797937,
                      0.510977961,
                      0.16315903,
                      1.562636622,
                      0.719269668,
                      0.398024938,
                      0.931474791,
                      1.332001248,
                      1.251857165,
                      1.49107416,
                      0.547994663,
                      0.785012712,
                      -0.506751667,
                      -0.416046869,
                      0.402595583,
                      -0.705561784,
                      -0.854565633,
                      -0.847802825,
                      -0.87445892,
                      -1.330698763,
                      -1.022396797,
                      -1.018015333,
                      -0.881538437,
                      -1.148443578
), section=c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4), rep=c(1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6))

df=df_old
df$log2_fold_change=df_old$log2_fold_change-mean(df_old[c(19,20,21,22,23,24),1])


df$section=as.factor(df$section)
df$rep=as.factor(df$rep)

model2=aov(log2_fold_change~section, data=df)
model3=lm(log2_fold_change~section, data=df)

p<-ggplot(df, aes(x=as.factor(section), y=log2_fold_change))  

p<- p + geom_boxplot()

#p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
#                 geom = "crossbar", width = 0.4)

p

#pairwise ttests
pairwise.t.test(df$log2_fold_change, df$section, p.adj = "bonf")

