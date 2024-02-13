# LaPorte 2024-02-08
# This script can be used for between environment prediction, including the 2- and 13-gene version. Additional notes about this on the accompanying script.
#this version of the script created 06/29/21, is an adaptation so that 1. It can be run on farm dynamically and 2. for loops are parallelized

#####################################################
#
#           Analysis Libraries
#
#####################################################
library(glmnet)
library(rrBLUP)
library('MASS')
library(gplots)
library("compiler")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

######################################################
#
#  Libraries for parallelization, added by MF 6/29/21
#
######################################################

library(future.apply)
library(parallel)
library(rlist)

setwd("/home/mflapo/")

############### Using Alex Lipka's code for 10x CV, modified #################
#                      Assign proper values to each input below
# Set Y equal to single trait of interest and Geno equal to hapmap genotype file
# This is automatically defined if using GLMNet.tenfoldCV() function

#STEP 1:
#Indicate whether or not you want to filter out SNPs, and specify regression type.
#Note: Gene information is read in on Line 21. This particular file is for the list of SNPs within +/- 250 kb of 58 candidate genes.
#To look for a list of other candidate genes (e.g., biosynthetic pathway genes), create similar file formatted exaclty like "Gene_Information.txt"

#######################################################################################


filter.out.SNPs = FALSE
filter.by.MAF = FALSE

method = 'RR-BLUP' #Options: GLMNet or RR-BLUP
alpha = 1 #only need to specify alpha if running GLMNet (1 for Lasso, 0 for ridge, intermediate values for elastic net analysis.)
#######################################################################################
# G is the matrix of SNPs
# Y is the vector of phenotypes
# K is the kinship matrix
# G, Y, and X have the entity ids removed

#Step 1: Set data directory and import files
#######################################################################################
mydataPath = getwd()

#Genotypic Data
myG = read.table(paste(mydataPath,"/CAM_geno_122320.hmp.txt",sep=""),
                 sep='\t',head=FALSE,check.names = FALSE)
#dim(myG)
#myG[1:15,1:15]

#Phenotypic Data
Y1 = read.table(paste("/home/mflapo/projects/REDO/Between-Enviro/Between-Enviro/AF12_averagesPhenotype_PATTAMA_v5_carotenoids_fromKevin.missing_cis_isomers.txt",sep=""),
                 head = TRUE,check.names = FALSE)

Y2 = read.table(paste("/home/mflapo/projects/REDO/Between-Enviro/Between-Enviro/AF13_averagesPhenotype_PATTAMA_v5_carotenoids_fromKevin.missing_cis_isomers.txt",sep=""),
                   head = TRUE,check.names = FALSE)
#added 2022-1-19

interum <- na.omit(merge(Y1, Y2, by = '<Trait>'))


common_1 <- cbind(interum[,1],interum[,2:6])
common_2 <- cbind(interum[,1],interum[,7:11])

colnames(common_1) <- c("<Trait>","lut","zea","bcx","bc","proa")
colnames(common_2) <- c("<Trait>","lut","zea","bcx","bc","proa")

myY_1 <- common_1
#changed ML 2/3/22
myY_2 <- common_2

n.traits_1 = ncol(myY_1)-1
n.traits_2 = ncol(myY_2)-1

### the way that it is written, must be the same size
n.indiv_phenos_1 = nrow(myY_1)
n.indiv_phenos_2 = nrow(myY_2)

#######################################################################################
##For obtaining a subset of SNPs. NOTE that this is done before SNPs with MAF < 0.05 are filtered out.

if(filter.out.SNPs){
  #Read in the information on the bp position of the candidate genes
  gene.information <- read.delim("/Gene_Information", head = TRUE)
  #Make the bp positions numeric so that mathematical logical statments can be assessed
  bp.numeric <- as.numeric(as.vector(myG[,4]))
  #Obtain all SNPs within +/- 250 kb of candidate genes
  logical.statement.250kb <- NULL
  for(j in 1:nrow(gene.information)) logical.statement.250kb <- c(logical.statement.250kb, which((myG[,3] == gene.information[j,2]) & (bp.numeric>gene.information[j,5]) & (bp.numeric<gene.information[j,6])) )
  # Use logical.statement.250kb to filter out SNPs that are not within +/- 250 kb of candidate genes
  myG.250kb <- myG[c(1,logical.statement.250kb),]
}

#######################################################################################
#End code for obtaining a subset of SNPs

#######################################################################################
#Begin for-loop to analyze one trait at a time (each has different numbers of individuals due to missing phenotypic data)
path.for.results <- "/home/mflapo/projects/REDO/Between-Enviro/"

#Geno = myG.250kb # if using filter.out.SNPs
myG_genos = myG
n.indiv = nrow(myG_genos)-1
n.SNPs = ncol(myG_genos)-1

r.gy.all = NULL                   #initialize results vector

#######################################################################################
#
#         Loop to apply each to trait
#
#######################################################################################
#for (i in 1:n.traits) {
GP_master_script <- function(i) {
  current_trait_name <- colnames(myY_1)[i+1]
  Y_1 = myY_1[,c(1,1+i)]                                               #extracts variety names and phenotypic values for trait i
  Y_2 = myY_2[,c(1,1+i)]                                             #extracts variety names and phenotypic values for trait i
  NAlines_1 = which(is.na(Y_1[,2]))
  NAlines_2 = which(is.na(Y_2[,2])) #finds row numbers having missing phenotypic data
  NAlines_total <- unique(append(NAlines_1, NAlines_2))
  ##CHD moved random.sample_adjustedpertrait definition line (from random.sample, removing NA lines) down below filter.by.MAF section so that it is defined after random.sample is refined (otherwise will return an error)
  Geno = myG_genos
  
  ########### Numericalizing Genotype Data ###################
  #use phenos for first trait as "filler" Y--dummy needed as CV for GAPIT.HapMap to work.
  CV_1 =Y_1
  CV_2 =Y_2
  CV_1[,2]=1
  CV_2[,2]=1
  colnames(CV_1)=c("taxa","overall")
  colnames(CV_2)=c("taxa","overall")
  
  hm=GAPIT.HapMap(G = Geno,SNP.effect="Add",SNP.impute="Major") #We'll want to verify that this numericalization function is still working properly; currently returning warnings.
  # View(hm$GD[1:10,1:10])
  GK <- cbind(hm$GT, hm$GD)
  print("debug1")
  #####################################
  if(filter.by.MAF){
    #Obtain the mafs of all SNPs
    #Total number of lines
    ns <- nrow(hm$GD)
    #ns = n.indiv
    #Sum of the allele scores for each SNP
    ss <- apply(hm$GD, 2, sum)
    #ss <- apply(myG_genos_numMatr, 2, sum)
    #Combine two situations: one where the allele coded as "2" is major; one where "0" is coded as major.
    maf.matrix <- rbind((.5*ss/ns), (1-(0.5*ss/ns)))
    #Copy the minor allele frequencies for all SNPs
    maf <- apply(maf.matrix, 2, min)
    #Find out which SNPs have MAF < 0.05
    snps.below.0.05.maf <- which(maf < 0.05)
    # Remove these SNPs from hm$GD
    hm.GD.without.snps.below.0.05.maf <- hm$GD[,-snps.below.0.05.maf] #CHD tried uncommenting this line 12/17/2020
    #myG_genos_numMatr.without.snps.below.0.05.maf <- myG_genos_numMatr[,-snps.below.0.05.maf] #CHD tried commenting this line 12/17/2020
    GK <- cbind(hm$GT, hm.GD.without.snps.below.0.05.maf)
  }
  ###############################
  

  

  
  
  qc_1=GAPIT.QC(Y = Y_1, GT = hm$GT, CV = CV_1, GK = GK)
  qc_2=GAPIT.QC(Y = Y_2, GT = hm$GT, CV = CV_2, GK = GK)
  
  y_1 <- as.matrix(qc_1$Y[-1])
  y_2 <- as.matrix(qc_2$Y[-1])
  
  G <- as.numeric(qc_1$GK[,-1])
  G <- matrix(G, nrow(y_1), ncol(qc_1$GK[,-1]))
  G <- G - 1
  Geno <- G
  
  Y_1 <- as.data.frame(y_1)
  Y_2 <- as.data.frame(y_2)
  
  
  #Randomly sort the number of lines, and subdivide them into ten subgroups
  sample.size <- n.indiv_phenos_1 #number of taxa.kept in QC function
  sequence.sample <- rep(1:sample.size)
  #set.seed(312) #so that the same random sample is drawn each time (remove this line if conducting iterations with different random samples)
  seeder <- 312

  
  ############################################################################################
  #Run five-fold cross validation
  

  
  r.gy <- numeric()

  
  for(m in 0:19){
      
    set.seed(seeder)
    seeder <- seeder + 1
      
    random.sample <- sample(1:sample.size, replace = FALSE)
      
    random.sample_adjustedpertrait_1 = random.sample[!(random.sample %in% NAlines_total)]  #removes those varieties from random sample vector #CHD moved this line here on 8/28/2020
    random.sample_adjustedpertrait_2 = random.sample[!(random.sample %in% NAlines_total)]  #removes those varieties from random sample vector #CHD moved this line here on 8/28/2020
      
    sample.size_adjustedpertrait_1 = length(random.sample_adjustedpertrait_1)
    sample.size_adjustedpertrait_2 = length(random.sample_adjustedpertrait_2)
     
    random.sample_eliminatedYs_1 = setdiff(random.sample,random.sample_adjustedpertrait_1)
    random.sample_eliminatedYs_2 = setdiff(random.sample,random.sample_adjustedpertrait_2)
      
    indices_eliminated_1 = sort(unique(c(NAlines_total,random.sample_eliminatedYs_1)))
    indices_eliminated_2 = sort(unique(c(NAlines_total,random.sample_eliminatedYs_2)))
      
    print(sample.size_adjustedpertrait_1)
      
    increment_1 = ceiling(sample.size_adjustedpertrait_1/5)                                  #defines how many individuals will be in each fold (ceiling finds the nearest integer without going over)
    increment_2 = ceiling(sample.size_adjustedpertrait_2/5)
      
    iteration_val <- m
    #have a "for" loop, start it at 0, and end it at 9 for 10x CV or 4 for 5x CV
    for(j in 0:4){
      cross_val_val <- j
      pred_1 <- random.sample_adjustedpertrait_1[((increment_1*j)+1):min(((increment_1*j)+increment_1) , sample.size_adjustedpertrait_1)]                    #e.g. For Trait 1 (200 lines total) you'll have 160 individuals in training population
      pred_2 <- random.sample_adjustedpertrait_2[((increment_2*j)+1):min(((increment_2*j)+increment_2) , sample.size_adjustedpertrait_2)]
      
      train_1 <- random.sample_adjustedpertrait_1[-(((increment_1*j)+1):min(((increment_1*j)+increment_1) , sample.size_adjustedpertrait_1))]                #e.g. For Trait 1 (200 lines total) you'll have 40 individuals in prediction set
      train_2 <- random.sample_adjustedpertrait_2[-(((increment_2*j)+1):min(((increment_2*j)+increment_2) , sample.size_adjustedpertrait_2))]
      
      ## For loop to get rid of the indices corresponding to missing data
      count <- 0
      for(k in indices_eliminated_1){
        if(count > 0) k <- k-count #if you need to subtract more than 1 for more than two lines missing, replace 1 with count
        
        if(any(train_1 > k)) {
          train.temp.1 <- train_1[which(train_1 > k)]
          train.temp.2 <- train_1[-which(train_1 > k)]
          try(train.temp.1 <- train.temp.1 - 1)
          train_1 <- c(train.temp.1, train.temp.2)
        }
        
        if(any(pred_1 > k)) {
          pred.temp.1 <- pred_1[which(pred_1 > k)]
          pred.temp.2 <- pred_1[-which(pred_1 > k)]
          try(pred.temp.1 <- pred.temp.1 - 1)
          pred_1 <- c(pred.temp.1, pred.temp.2)
        }
        
        count <- count+1
      } #End for-loop of k (for index correction)
      
      count <- 0
      for(k in indices_eliminated_2){
        if(count > 0) k <- k-count #if you need to subtract more than 1 for more than two lines missing, replace 1 with count
        
        if(any(train_2 > k)) {
          train.temp.1 <- train_2[which(train_2 > k)]
          train.temp.2 <- train_2[-which(train_2 > k)]
          try(train.temp.1 <- train.temp.1 - 1)
          train_2 <- c(train.temp.1, train.temp.2)
        }
        
        if(any(pred_2 > k)) {
          pred.temp.1 <- pred_2[which(pred_2 > k)]
          pred.temp.2 <- pred_2[-which(pred_2 > k)]
          try(pred.temp.1 <- pred.temp.1 - 1)
          pred_2 <- c(pred.temp.1, pred.temp.2)
        }
        
        count <- count+1
      } #End for-loop of k (for index correction)
      
      Geno_1.train <- Geno[train_1,]
      Geno_2.train <- Geno[train_2,]
      Pheno_1.train <- Y_1[train_1,]
      Pheno_2.train <- Y_2[train_2,]
      
      Geno_1.pred <- Geno[pred_1,]
      Geno_2.pred <- Geno[pred_2,]
      Pheno_1.pred <- Y_1[pred_1,]
      Pheno_2.pred <- Y_2[pred_2,]
      
      ##################
      
      ## Next 3 lines used for GLMNet. alpha=1 is lasso penalty, alpha=0 is ridge regression penalty, middle ground is elastic net. two commented cv.fit functions are alternatives for model-fitting.
      if(method == 'GLMNet'){
        #cv.fit <- cv.glmnet(Geno.train, Pheno.train, nlambda=250, nfold=20)
        #cv.fit <- cv.glmnet(Geno.train, Pheno.train, lambda=lambda_seq)
        cv.fit <- cv.glmnet(Geno.train, Pheno.train, alpha = alpha, nlambda=250)
        ans <- predict(cv.fit, Geno.pred, s="lambda.min")
        #r.gy <- c(r.gy, cor(ans, Pheno.pred))
        method.specs = paste(method,"_",alpha,sep='')
      }
      
      ## Next 2 lines used for rrBLUP
      if(method == "RR-BLUP"){
        ans_1 <- kinship.BLUP(y=Pheno_1.train,G.train=Geno_1.train,G.pred=Geno_1.pred,K.method="RR")
        ans_2 <- kinship.BLUP(y=Pheno_2.train,G.train=Geno_2.train,G.pred=Geno_2.pred,K.method="RR")
        
        acc_temp_1 <- cor(ans_1$g.pred,Pheno_1.pred)
        acc_temp_2 <- cor(ans_2$g.pred,Pheno_2.pred)
        acc_temp_1v2 <- cor(ans_1$g.pred,Pheno_2.pred)
        acc_temp_2v1 <- cor(ans_2$g.pred,Pheno_1.pred)
        temprow <- c(current_trait_name, iteration_val, cross_val_val, acc_temp_1, acc_temp_2, acc_temp_1v2, acc_temp_2v1)
        r.gy <- rbind(r.gy, temprow)
        method.specs = method
      }
      
      
    } #End for-loop of j (for CV)
  }
  # Once the loop is over, output the values of the correlation coefficients, as well as their means and standard deviations
  #r.gy = c(r.gy, mean(r.gy), sd(r.gy))
  #r.gy.all = rbind(r.gy.all, r.gy)
  # This is for 5x CV
  
  r.gy.allwithnames <- r.gy
  
  #r.gy.allwithnames = cbind(traitnames_abbrevs_nohead, r.gy.all)
  #r.gy.allwithnames = cbind(traitnames_abbrevs, r.gy.all)
  
  r.gy.output <- as.matrix(r.gy.allwithnames)
  #colnames(r.gy.output) <- c("trait","r_CV1", "r_CV2", "r_CV3", "r_CV4", "r_CV5", "r_CV6", "r_CV7", "r_CV8", "r_CV9", "r_CV10","r_avg", "r_std")
  #colnames(r.gy.output) <- c("trait","r_CV1", "r_CV2", "r_CV3", "r_CV4", "r_CV5","r_avg", "r_std")
  colnames(r.gy.output) <- c("trait","iteration", "fold", "accuracy1","accuracy2","accuracy1v2","accuracy2v1")
  
  write.table(r.gy.output, paste(path.for.results,"/GP.Results_Maize.CAM_Between_AF12-AF13_",method.specs,"-",current_trait_name,"_20-iterations.txt",sep=''), quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  
  
  return(r.gy)
  
} #End function "GP_master_script"
#^^Previously: End for-loop of i (for traits)

########################################################################
#
#   Call the new function with future_lapply
#
########################################################################

looper = 1:n.traits_1
r.gy.all <- mclapply(looper, GP_master_script)
r.gy.all.df <- list.rbind(r.gy.all)
summary(r.gy.all.df)

########################################################################
#
#   Write out results
#
########################################################################
r.gy.allwithnames <- r.gy.all.df
#traitnames_abbrevs = colnames(myY)
#traitnames_abbrevs_nohead = traitnames_abbrevs[-1]
#r.gy.allwithnames = cbind(traitnames_abbrevs_nohead, r.gy.all.df)
#r.gy.allwithnames = cbind(traitnames_abbrevs, r.gy.all)
r.gy.output <- as.matrix(r.gy.allwithnames)
#colnames(r.gy.output) <- c("trait","r_CV1", "r_CV2", "r_CV3", "r_CV4", "r_CV5", "r_CV6", "r_CV7", "r_CV8", "r_CV9", "r_CV10","r_avg", "r_std")
#colnames(r.gy.output) <- c("trait","r_CV1", "r_CV2", "r_CV3", "r_CV4", "r_CV5","r_avg", "r_std")
colnames(r.gy.output) <- c("trait","iteration", "fold", "accuracy1","accuracy2","accuracy1v2","accuracy2v1")
write.table(r.gy.output, paste(path.for.results,"Results-GP-RRBet-AF12-AF13-5-20.txt",sep=''), quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)

