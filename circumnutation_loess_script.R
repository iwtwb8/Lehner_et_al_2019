library(ggplot2)
library(splines)


#read data in and initial plotting
xy_FN1=read.csv("1507_1503_04_26_2018_FN_179_1.csv", header = TRUE, row.names=1)
xy_FN2=read.csv("1507_1503_04_26_2018_FN_179_2.csv", header = TRUE, row.names=1)
xy_WT1=read.csv("1507_1503_04_26_2018_x_Kitaake_1.csv", header = TRUE, row.names=1)
xy_WT2=read.csv("1507_1503_04_26_2018_x_Kitaake_3.csv", header = TRUE, row.names=1)
xy_FN3=read.csv("1507_1503_05_04_2018_FN_179_1.csv", header = TRUE, row.names=1)
xy_FN4=read.csv("1507_1503_05_04_2018_FN_179_2.csv", header = TRUE, row.names=1)
xy_WT3=read.csv("1507_1503_05_04_2018_x_Kitaake_1.csv", header = TRUE, row.names=1)
xy_WT4=read.csv("1507_1503_05_04_2018_x_Kitaake_2.csv", header = TRUE, row.names=1)
xy_FN5=read.csv("1507_1503_05_07_2018_FN_179_2.csv", header = TRUE, row.names=1)
xy_WT5=read.csv("1507_1503_05_07_2018_x_Kitaake_1.csv", header = TRUE, row.names=1)
xy_WT6=read.csv("1507_1503_05_07_2018_x_Kitaake_2.csv", header = TRUE, row.names=1)
xy_WT7=read.csv("1507_1503_05_07_2018_x_Kitaake_3.csv", header = TRUE, row.names=1)
xy_WT8=read.csv("1507_1503_03_13_2018_x_Kitaake_3.csv", header = TRUE, row.names=1)
xy_WT9=read.csv("1507_1503_03_13_2018_x_Kitaake_2.csv", header = TRUE, row.names=1)
xy_WT10=read.csv("1507_1503_03_13_2018_x_Kitaake_1.csv", header = TRUE, row.names=1)
xy_FN6=read.csv("1507_1503_03_13_2018_FN_179_2.csv", header = TRUE, row.names=1)
xy_FN7=read.csv("1507_1503_03_13_2018_FN_179_1.csv", header = TRUE, row.names=1)


#create single data.frame
xys=list(xy_FN1,xy_FN2, xy_FN3, xy_FN4, xy_FN5, xy_FN6, xy_FN7, xy_WT1, xy_WT2,xy_WT3,xy_WT4,xy_WT5,xy_WT6,xy_WT7, xy_WT8, xy_WT9, xy_WT10)
#4_26, 5_04, 5_07, 3_13

#estimate sequence specific distance factor 
#601.37 is the apparent mm of the back of the box in our test image. The denominator is the number of pixels across the back of the box in each set.
mpp=c(60.137/1208,60.137/1208,60.137/1220,60.137/1220, 60.137/1218, 60.137/1158, 60.137/1158,60.137/1208,60.137/1208,60.137/1220,60.137/1220,60.137/1218,60.137/1218,60.137/1218,60.137/1158,60.137/1158,60.137/1158) #mm per pixel for each set

#create data.frames for period and amplitudes
df_amps=data.frame(matrix(ncol=17))
df_pers=data.frame(matrix(ncol=17))
colnames(df_amps)=c("xy_FN1","xy_FN2", "xy_FN3", "xy_FN4", "xy_FN5", "xy_FN6","xy_FN7","xy_WT1", "xy_WT2","xy_WT3","xy_WT4","xy_WT5","xy_WT6","xy_WT7", "xy_WT8","xy_WT9","xy_10")
colnames(df_pers)=c("xy_FN1","xy_FN2", "xy_FN3", "xy_FN4", "xy_FN5", "xy_FN6","xy_FN7", "xy_WT1", "xy_WT2","xy_WT3","xy_WT4","xy_WT5","xy_WT6","xy_WT7", "xy_WT8","xy_WT9","xy_10")

#c is a count variable
c=0 
for (i in xys) { #loop through all data.frames in list
  c=c+1 #counts which element of the list we are on for later use.
  xy_df=i[,2:3] #columns 2 and 3 are the coordinates

  #loess regression to estimate center line of root
  loessMod10 <- loess(BX~BY, data=xy_df, span=0.20) #use Loess to estimate centerline. Span can be tuned. BX and BY are from imageJ
  plot(BX~BY, data=xy_df, main=colnames(df_amps[c]))
  lines(xy_df$BY, predict(loessMod10, data.frame(BY = xy_df$BY)))
  
  
  #correct for drift
  BX_corrected=0 #will be the BX adjusted to the point along the loess prediction surface nearest to the observed point
  BX_min=0 #keeps track of shortest distance
  dist=0 #placeholder variable for distance
  
  #count is another counter variable
  count=0
  for (obs in xy_df$BY){
    count=count+1 #track which observation is being analyzed
    BX_min=xy_df[count,1]-predict(loessMod10, xy_df[count,])
    for (offset in c(-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6)) { #the pixel offset range we're scanning over. Might want to make larger.
      if(!is.na(predict(loessMod10, obs+offset))) #deal with predicting at the ends where offsets will produce NAs
        dist=euc.dist(xy_df[count,],c(predict(loessMod10, obs+offset),obs+offset)) #calculate distance of each point to point along curve with offset
        #print(dist)
        if(dist<abs(BX_min)) {
          BX_min=BX_min/abs(BX_min)*dist #save offset to minimum if it's smaller than the previously found minimum. Preserve sign.
        }
    }
    BX_corrected[count]=BX_min
  }
  
  #plot(BX_corrected, ylim=c(-25,25), main=colnames(df_amps)[c])
  plot(BX_corrected*mpp[c]~c(seq(0,54,.25)), ylim=c(-1.5,1.5), main=colnames(df_amps)[c], xlab="hours",ylab="distance from center line (millimeters)")
  
  #plot a spline over the corrected values
  spline.lm <- lm((BX_corrected*mpp[c]) ~ ns(seq(0,54,.25), df=100), data=xy_df)
  lines(seq(0,54,.25), predict(spline.lm), lwd=2, col='red')
  
  #find indices of local maxima
  loc_max=find_peaks(BX_corrected+rnorm(length(BX_corrected), mean=0, sd= .0001))
  
  #find indices of local minima
  loc_min=find_peaks(-BX_corrected+rnorm(length(BX_corrected), mean=0, sd= .0001))
  
  #plot basedon average
  len=min(length(loc_max), length(loc_min))
  if(length(loc_max)<length(loc_min)) {
    points=loc_max
  } else {
    points=loc_min
  }
  amplitudes=abs((BX_corrected[loc_max[1:len]]-BX_corrected[loc_min[1:len]])/2)
  periods=diff(points)

  df_amps[1:18,c]=tapply(amplitudes, cut(points, seq(1,217,12)), mean)
  df_amps[1:18,c]= df_amps[1:18,c] * mpp[c]
  df_pers[1:6,c]=tapply(periods, cut((points[1:length(points)-1]), seq(1,217,36)), mean)
}

#plot mutants AMPLITUDES
boxplot(t(df_amps[,1:7]), ylim=c(0,max(df_amps, na.rm=T)), ylab="average amplitude (millimeters)")

#plot wildtypes AMPLITUDES
boxplot(t(df_amps[,c(8,9,10,12,13,14,15,16,17)]), ylim=c(0,max(df_amps, na.rm=T)),ylab="average amplitude (millimeters)")

#plot mutants PERIODS
boxplot(na.exclude(t(df_pers[,1:7]))/4, ylab="average period (hours)")

#plot wildtypes PERIODS
boxplot(na.exclude(t(df_pers[,c(8,9,10,12,13,14,15,16,17)]))/4, ylab="average period (hours)")



#euclidean distance function. From user Shambho. https://stackoverflow.com/questions/24746892/how-to-calculate-euclidian-distance-between-two-points-defined-by-matrix-contain
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))


#local maximum function. From user stas g. https://stackoverflow.com/questions/34205515/finding-local-maxima-and-minima-in-r
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}


