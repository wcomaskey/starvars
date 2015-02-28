#####################################################################
##### run this script to perform AL from start to finish on ASAS data
## and save all relevant output, classifications + performance metrics
######################################################################

####################################################
# load needed packages and files

set.seed(1)

source("utils_classify.R")
## source(paste(path,"R/class_cv.R",sep=""))
## source("/Users/jwrichar/Documents/CDI/R/missForest.R")
library(randomForest)
library(nnet)

####################################################
####################################################
# load data

library(foreign)

####################################################
# Load TRAINING data
cat("Loading ASAS + Hipp + OGLE Training data\n")
feat.train = read.table("data/training_set_features.dat",header=TRUE,sep=',')
ID.train = feat.train$ID.train
exclude = c("ID.train","qso_log_chi2_qsonu","qso_log_chi2nuNULL_chi2nu","freq_n_alias")
feat.train = feat.train[,-which(names(feat.train) %in% exclude)]
ra.dec.train = read.table("data/training_set_ra_dec.dat",header=TRUE,sep=',')
class.train = read.table("data/training_set_id_class.dat",header=FALSE,sep='\t')[,2]
n = length(ID.train)

cat(paste("Training data loaded: ",n," sources\n"))


####################################################
# load LINEAR data
cat("Loading LINEAR data\n")

#lineardat=read.arff(file="data/200k_combo.arff")
lineardat = read.table("data/200k_final_combo.csv",sep=",",header=TRUE)
lineardat$freq_rrd = rep(0,dim(lineardat)[1])
feat.tmp = data.frame(lineardat)
feat.linear = data.frame(matrix(0,dim(feat.tmp)[1],dim(feat.train)[2]))
  for(ii in 1:dim(feat.train)[2]){
    feat.linear[,ii]= feat.tmp[,names(feat.tmp)==names(feat.train)[ii]]
  }
colnames(feat.linear) = names(feat.train)
feat.linear[is.na(feat.linear)] = 0

# fix the ? problem
fixfeat = function(feature){
  feature = paste(feature)
  feature[feature=="?"] = 0
  feature = as.numeric(feature)
  return(feature)
}

feat.linear$amplitude = fixfeat(feat.linear$amplitude)
feat.linear$freq_frequency_ratio_21 = fixfeat(feat.linear$freq_frequency_ratio_21)
feat.linear$freq_frequency_ratio_31 = fixfeat(feat.linear$freq_frequency_ratio_31)
feat.linear$freq_signif_ratio_21 = fixfeat(feat.linear$freq_signif_ratio_21)
feat.linear$freq_signif_ratio_31 = fixfeat(feat.linear$freq_signif_ratio_31)
feat.linear$freq1_lambda = fixfeat(feat.linear$freq1_lambda)
feat.linear$max_slope = fixfeat(feat.linear$max_slope)
feat.linear$median_absolute_deviation = fixfeat(feat.linear$median_absolute_deviation)
feat.linear$p2p_scatter_2praw = fixfeat(feat.linear$p2p_scatter_2praw)
feat.linear$p2p_scatter_over_mad = fixfeat(feat.linear$p2p_scatter_over_mad)
feat.linear$p2p_scatter_pfold_over_mad = fixfeat(feat.linear$p2p_scatter_pfold_over_mad)
feat.linear[is.na(feat.linear)] = 0

ID.linear = lineardat$source_id
N = dim(feat.linear)[1]

## add RRLd feature
feat.linear$freq_rrd = ifelse(abs(feat.linear$freq_frequency_ratio_21 - 0.746) < 0.0035 | abs(feat.linear$freq_frequency_ratio_31 - 0.746) < 0.0035, 1,0)

##############################################
# read in mastermain file (ra/dec plus colors)
ra.dec.linear = read.table("data/masterMain.dat.txt",header=TRUE)
# align the mastermain table to the feature table
infeature = which(ra.dec.linear$objectID %in% ID.linear)
ra.dec.linear = ra.dec.linear[infeature,]
orderfeat = order(ID.linear)
ordermain = order(ra.dec.linear$objectID)
# sort features and main table by object ID
feat.linear = feat.linear[orderfeat,]
ID.linear = ID.linear[orderfeat]
ra.dec.linear = ra.dec.linear[ordermain,]


# compute colors for the LINEAR data
source("color_convert_ugriz_BVR.R")
linear.bvr = color.transform(ra.dec.linear$uMod,ra.dec.linear$gMod,ra.dec.linear$rMod,
  ra.dec.linear$iMod,ra.dec.linear$zMod)

feat.linear$color_diff_bj = linear.bvr$B - ra.dec.linear$J
feat.linear$color_diff_rj = linear.bvr$R - ra.dec.linear$J
feat.linear$color_diff_vj = linear.bvr$V - ra.dec.linear$J


cat(paste("LINEAR data loaded: ",n," sources\n"))


####################################################
# load LINEAR labels on the 7k
## cat("Loading LINEAR Labels for the 7k data set\n")
## linear.eb = read.table("data/EBs_IDs.txt",header=FALSE)[,1]
## linear.rrc = read.table("data/data_list_RRc.txt",header=TRUE)[,1]
## linear.rrab = read.table("data/RRab.dat",header=FALSE)[,1]
## linear.sxe = read.table("data/2013_02_04_SXPhe.dat",header=TRUE)[,1]

## linear.class = rbind(cbind(linear.eb, rep("EB",length(linear.eb))),
##   cbind(linear.rrab, rep("RRab",length(linear.rrab))),
##   cbind(linear.rrc, rep("RRc",length(linear.rrc))),
##   cbind(linear.sxe, rep("SXPhe",length(linear.sxe))))

### 2013-06-17 USING THE LINEARattributesFinalApr2013.dat file
# 1 = ab RR Lyr (2612/2923); 2 = c RR Lyr
# (864/990); 4 = Algol-like with 2 minima (342/357); 5 = contact binary (2246/2385);
# 6 = delta Scu/SX Phe (82/112);
linear.class = read.table("data/LINEARattributesFinalApr2013.dat")[,c(14,13)]
linear.class[,2] = factor(linear.class[,2],labels=c("RRab","RRc","Algol","Contact","DScu/SXPhe"))

cat(paste("Number of manually labeled LINEAR objects: ",dim(linear.class)[1]," sources\n"))


####################################################
## add LINEAR - trainingset overlap data to training set
cat("Adding LINEAR training objects that overlap with Existing Training Set (ASAS + Debosscher)\n")
counter = 0
for(ii in 1:dim(ra.dec.train)[1]){
  ind = which(round(ra.dec.linear$ra,2)==round(ra.dec.train[ii,2],2) & round(ra.dec.linear$dec,2)==round(ra.dec.train[ii,3],2))
  if(length(ind)>0){
    ind1 = which(ID.linear == ra.dec.linear$objectID[ind])
    ind2 = which(ID.train == ra.dec.train[ii,1])
    if(length(ind2)==1){
      counter = counter+1
      cl.tmp = paste(class.train[ind2])
      cat("Exchanging ID: ",ID.train[ind2]," Class: ",cl.tmp,"\n")
      ## remove old Training (ASAS / Deb) data
      feat.train = feat.train[-ind2,]
      class.train = class.train[-ind2]
      ID.train = ID.train[-ind2]
      ## add new LINEAR data
      feat.train = rbind(feat.train,feat.linear[ind1,])
      class.train = c(paste(class.train),paste(cl.tmp))
      ID.train = c(ID.train,ID.linear[ind1])
    }
  }
}
class.train = factor(class.train)
nn = length(class.train)
# total size: 
print(table(class.train[(nn-counter+1):nn])[table(class.train[(nn-counter+1):nn])>0])
cat("number added:",counter,"\n")
# 1 RRL FM, LINEAR ID 6196048 / 
# Exchanging ID:  230655  Class:  g. RR Lyrae FM 



################################################
################################################
# fit classifier on ASAS + Debosscher training set
################################################
cat("Final training set size:",dim(feat.train)[1],"\n") 
cat("Number of unique sources:",length(unique(ID.train)),"\n")

cat("Fitting Random Forest\n")
rf.tr = randomForest(x=as.matrix(feat.train),y=class.train,mtry=10,ntree=1000)
pred.linear = predict(rf.tr,newdata = feat.linear)
pred.linear.prob = predict(rf.tr,newdata = feat.linear,type='prob')


# save RF object (for later use; use load(file) to fire it up in future sessions)
save(rf.tr,file="data/linear_randomForest.Rdat")

#  performance metrics
p.hat = apply(pred.linear.prob,1,max)
mean.p.hat = mean(p.hat)
perc.conf = sum(p.hat>0.5)/N
cat("Final performance metrics:\n")
cat("Mean max p-hat:",mean.p.hat,"\n")
cat("Percent w/ prob > 0.5:",perc.conf,"\n")
#Mean max p-hat: 0.2065849 
#Percent w/ prob > 0.5: 0.01955851

  pdf("plots/linear_hist_maxprob.pdf",height=7,width=7)
hist(p.hat,col=4, main="Histogram of Maximum RF P(class)",xlab="max P(class)")
  dev.off()


# plot confmat for objects w/ manual labels
# get everything in right order to make table
inval = which(ID.linear %in% as.numeric(linear.class[,1]))
class.true = class.pred = NULL
for(ii in 1:length(inval)){
  ind = which(linear.class[,1] == ID.linear[inval[ii]])[1]
  class.true = c(class.true, paste(linear.class[ind,2]))
  class.pred = c(class.pred, paste(pred.linear[inval[ii]]))
}
linear.tab = table(as.factor(class.pred), as.factor(class.true))
  pdf("plots/linear_comparison_ml_manual.pdf",height=6,width=11)
par(mar=c(2.5,9,9,1.95))
plot.table(linear.tab,title="",cexprob=1,cexaxis=1.06)
mtext("Manually Identified Class",2,padj=-6.05,cex=2)
mtext("RF Predicted Class",3,padj=-5.8,cex=2)
  dev.off()

hist(p.hat[inval],col=4, main="Histogram of Maximum RF P(class)",xlab="max P(class)")



print(table(pred.linear))
 ##               a. Mira         b1. Semireg PV             b2. SARG A 
 ##                  4811                   6338                    121 
 ##            b3. SARG B                b4. LSP            c. RV Tauri 
 ##                  1752                    835                   8421 
 ##  d. Classical Cepheid     e. Pop. II Cepheid f. Multi. Mode Cepheid 
 ##                   855                     93                   2067 
 ##        g. RR Lyrae FM         h. RR Lyrae FO         i. RR Lyrae DM 
 ##                 10680                    932                   2166 
 ##        j. Delta Scuti             j1. SX Phe       k. Lambda Bootis 
 ##                  3605                      0                      0 
 ##        l. Beta Cephei      m. Slowly Puls. B       n. Gamma Doradus 
 ##                    37                      0                      0 
 ##       o. Pulsating Be        p. Per. Var. SG      q. Chem. Peculiar 
 ##                  4996                    970                      7 
 ##         r. Wolf-Rayet                r1. RCB     s1. Class. T Tauri 
 ##                  3085                   1962                    926 
 ## s2. Weak-line T Tauri             s3. RS CVn        t. Herbig AE/BE 
 ##                  3250                   2219                   4190 
 ##          u. S Doradus         v. Ellipsoidal         w. Beta Persei 
 ##                     0                      0                 109243 
 ##         x. Beta Lyrae        y. W Ursae Maj. 
 ##                 11054                   3385 

print(table(pred.linear[p.hat>0.4]))
 ##               a. Mira         b1. Semireg PV             b2. SARG A 
 ##                    55                     13                      0 
 ##            b3. SARG B                b4. LSP            c. RV Tauri 
 ##                     0                      0                      0 
 ##  d. Classical Cepheid     e. Pop. II Cepheid f. Multi. Mode Cepheid 
 ##                     5                      0                      4 
 ##        g. RR Lyrae FM         h. RR Lyrae FO         i. RR Lyrae DM 
 ##                  1951                    122                     18 
 ##        j. Delta Scuti             j1. SX Phe       k. Lambda Bootis 
 ##                     6                      0                      0 
 ##        l. Beta Cephei      m. Slowly Puls. B       n. Gamma Doradus 
 ##                     0                      0                      0 
 ##       o. Pulsating Be        p. Per. Var. SG      q. Chem. Peculiar 
 ##                     0                      0                      0 
 ##         r. Wolf-Rayet                r1. RCB     s1. Class. T Tauri 
 ##                     0                      0                      0 
 ## s2. Weak-line T Tauri             s3. RS CVn        t. Herbig AE/BE 
 ##                     0                      0                      0 
 ##          u. S Doradus         v. Ellipsoidal         w. Beta Persei 
 ##                     0                      0                    328 
 ##         x. Beta Lyrae        y. W Ursae Maj. 
 ##                    41                   1026 

####### MAKE SUMMARY FILES FOR A FEW CLASSES (THRESHOLD AT P>0.4)

# classical cepheids:
ind.cep = which(pred.linear=="d. Classical Cepheid" & p.hat > 0.4)
tab.cep = cbind(ID.linear[ind.cep], p.hat[ind.cep], 1/feat.linear$freq1_harmonics_freq_0[ind.cep], feat.linear$amplitude[ind.cep])
colnames(tab.cep) = c("ID","Probability","Period","Amplitude")
write.table(tab.cep, "tables/linear_classical_cepheid_high_prob.dat", row.names=FALSE,sep=",",quote=FALSE)


# multi-mode cepheids:
ind.mmcep = which(pred.linear=="f. Multi. Mode Cepheid" & p.hat > 0.4)
tab.cep = cbind(ID.linear[ind.mmcep], p.hat[ind.mmcep], 1/feat.linear$freq1_harmonics_freq_0[ind.mmcep], feat.linear$amplitude[ind.mmcep])
colnames(tab.cep) = c("ID","Probability","Period","Amplitude")
write.table(tab.cep, "tables/linear_multimode_cepheid_high_prob.dat", row.names=FALSE,sep=",",quote=FALSE)

# mira variables:
ind.mira = which(pred.linear=="a. Mira" & p.hat > 0.5)
tab.mira = cbind(ID.linear[ind.mira], p.hat[ind.mira], 1/feat.linear$freq1_harmonics_freq_0[ind.mira], feat.linear$amplitude[ind.mira])
colnames(tab.mira) = c("ID","Probability","Period","Amplitude")
write.table(tab.mira, "tables/linear_mira_high_prob.dat", row.names=FALSE,sep=",",quote=FALSE)

# R Cor Bor:
ind.rcb = which(pred.linear=="r1. RCB" & p.hat > 0.4)
tab.rcb = cbind(ID.linear[ind.rcb], p.hat[ind.rcb], 1/feat.linear$freq1_harmonics_freq_0[ind.rcb], feat.linear$amplitude[ind.rcb])
colnames(tab.rcb) = c("ID","Probability","Period","Amplitude")
write.table(tab.rcb, "tables/linear_rcb_high_prob.dat", row.names=FALSE,sep=",",quote=FALSE)

# RRd:
ind.rrd = which(pred.linear=="i. RR Lyrae DM" & p.hat > 0.5)
tab.rrd = cbind(ID.linear[ind.rrd], p.hat[ind.rrd], 1/feat.linear$freq1_harmonics_freq_0[ind.rrd], feat.linear$amplitude[ind.rrd])
colnames(tab.rrd) = c("ID","Probability","Period","Amplitude")
write.table(tab.rrd, "tables/linear_rrd_high_prob.dat", row.names=FALSE,sep=",",quote=FALSE)

## ########
## ## plot feature importance
## plotfi = FALSE
## if(plotfi){        
##   featimp = matrix(0,length(rf.tr$importance),5)
##   rownames(featimp) = rownames(rf.tr$importance)
##   featimp[,1] = rf.tr$importance
##   for(ii in 2:5){
##     rf.tmp =  randomForest(x=as.matrix(feat.train),y=class.train,mtry=17,ntree=1000,nodesize=1)
##     featimp[,ii] = rf.tmp$importance
##   }
##   featimp.mean = apply(featimp,1,mean)
##   featimp.sd = apply(featimp,1,sd)

##   featimp.sort = featimp.mean[sort(featimp.mean,index.return=TRUE,decreasing=T)$ix]
##   featimp.sdsort = featimp.sd[sort(featimp.mean,index.return=TRUE,decreasing=T)$ix]

## #  pdf("/Users/jwrichar/Documents/CDI/ASAS/plots/asas_rf_imp.pdf",height=6,width=9)
##   par(mar=c(5,14,1,1))
##   barplot(featimp.sort[20:1],horiz=TRUE,xlab="Mean Gini Decrease",names.arg=rep("",20),space=0.2,lwd=2,col=2,xlim=c(0,150))
##   axis(2,labels=names(featimp.sort)[20:1],at=0.7+(0:19)*1.2,tick=FALSE,las=2)
##   arrows(featimp.sort[20:1] - featimp.sdsort[20:1], 0.7+(0:19)*1.2, featimp.sort[20:1] + featimp.sdsort[20:1], 0.7+(0:19)*1.2, code=3, angle=90, length=0.05,lwd=2)
##   legend('bottomright',"Feature Importance",cex=2,bty="n",text.col=1)
## #  dev.off()
## }


######################################################
# compute confidence measure
######################################################

## DO DISTANCE-BASED OUTLIER MEASURE
# cat("Computing classifier confidence values\n")
#source(paste(path.linear,"R/compute_confidence1.R",sep=""))
#outlier.score = anomScore(feat.train,class.train, feat.asas, ntree=500, knn=2, metric="rf")
#write(outlier.score,paste(path.linear,"data/outScore_asas_class.dat",sep=""),ncolumns=1)

