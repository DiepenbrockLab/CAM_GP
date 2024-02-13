#LaPorte 2024-02-08
#This version of the script created 06/29/21 (and modified in 02/2024), is an adaptation so that 
#1. It can be run on an HPC cluster dynamically and 
#2. For loops are parallelized
#3. The code for all models are included


#####################################################
#
#           Import Libraries
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

############### Using Alex Lipka's code for 10x CV, modified  ##########################
#                      Assign proper values to each input below
# Set Y equal to single trait of interest and Geno equal to hapmap genotype file
# This is automatically defined if using GLMNet.tenfoldCV() function

#STEP 1:
#Indicate whether or not you want to filter out SNPs, and specify regression type.
#Note: Gene information is read in on Line 21. This particular file is for the list of SNPs within +/- 250 kb of 58 candidate genes.
#To look for a list of other candidate genes (e.g., biosynthetic pathway genes), create similar file formatted exaclty like "Gene_Information.txt"

#######################################################################################

#Set filter.out.SNPs true for 2- or 13- gene analysis!
filter.out.SNPs = FALSE
filter.by.MAF = FALSE

method = 'RR-BLUP' #Options: GLMNet (for LASSO or EN), RR-BLUP, RKHS

#only need to specify alpha if running GLMNet (1 for Lasso, 0 for ridge, intermediate values (0.1-0.9 in the case of this study) for elastic net analysis.)
alpha = 1 

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

#Phenotypic Data 

myY = read.table(paste(mydataPath,"/AF12_averagesPhenotype_PATTAMA_v5_carotenoids_fromKevin.StandardLabels.txt",sep=""), 
                 head = TRUE,check.names = FALSE)

n.traits = ncol(myY)-1
n.indiv_phenos = nrow(myY)

#######################################################################################
## For obtaining a subset of SNPs. NOTE that this is done before SNPs with MAF < 0.05 are filtered out.
## Note that this is used for the 2-gene and 13-gene analysis!

if(filter.out.SNPs){
  #Read in the information on the bp position of the candidate genes
  gene.information <- read.delim("/Gene_Information_8_Cand_Genes_AEL.txt", head = TRUE)
  
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
path.for.results <- "/home/mflapo/projects/REDO/Ridge"

# if using filter.out.SNPs to do a 2- or 13-gene analysis, un-comment the line below, and comment out "myG_genos = myG"
#Geno = myG.250kb 

myG_genos = myG
n.indiv = nrow(myG_genos)-1
n.SNPs = ncol(myG_genos)-1

#initialize results vector
r.gy.all = NULL                   

#######################################################################################
#
# The function that applies GP analysis to each trait (i)!
#
#######################################################################################

#extracts variety names and phenotypic values for trait i

GP_master_script <- function(i) {  
  current_trait_name <- colnames(myY)[i+1]
  Y = myY[,c(1,1+i)]  
  
  #finds row numbers having missing phenotypic data
  NAlines = which(is.na(Y[,2]))     
  
  ##CHD moved random.sample_adjustedpertrait definition line (from random.sample, removing NA lines) down below filter.by.MAF section so that it is defined after random.sample is refined (otherwise will return an error)
  Geno = myG_genos
  
  ########### Numericalizing Genotype Data ################### 
  #use phenos for first trait as "filler" Y--dummy needed as CV for GAPIT.HapMap to work.
  #note that this must be done for each trait. This portion can be time consuming.
  
  CV=Y
  CV[,2]=1
  colnames(CV)=c("taxa","overall")
  hm=GAPIT.HapMap(G = Geno,SNP.effect="Add",SNP.impute="Major") 
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
  
  
  qc=GAPIT.QC(Y = Y, GT = hm$GT, CV = CV, GK = GK)
  
  y <- as.matrix(qc$Y[-1])
  G <- as.numeric(qc$GK[,-1])
  G <- matrix(G, nrow(y), ncol(qc$GK[,-1]))
  G <- G - 1
  Geno <- G
  Y <- as.data.frame(y)
  
  #Randomly sort the number of lines, and subdivide them into ten subgroups
  sample.size <- n.indiv_phenos #number of taxa.kept in QC function
  sequence.sample <- rep(1:sample.size)
  
  #variable to be used for random seed
  seeder <- 312
  

  
  ############################################################################################
  #Run five-fold cross validation
  r.gy <- numeric()                   

  for(m in 0:19){
    #Set a seed so that the same random sample is drawn each time (remove this line if conducting iterations with different random samples). 
    #In this case, a new random seed is set each iteration.
    
    set.seed(seeder) 
    #set value for new seed for the next iteration
    seeder <- seeder + 1
      
    random.sample <- sample(1:sample.size, replace = FALSE)
    #removes those varieties from random sample vector #CHD moved this line here on 8/28/2020
    random.sample_adjustedpertrait = random.sample[!(random.sample %in% NAlines)]  
      
    sample.size_adjustedpertrait = length(random.sample_adjustedpertrait)
    random.sample_eliminatedYs = setdiff(random.sample,random.sample_adjustedpertrait)
    indices_eliminated = sort(unique(c(NAlines,random.sample_eliminatedYs)))
    print(sample.size_adjustedpertrait)
    
    #this needs to be addressed if a different cross-validation strategy is used.
    increment = ceiling(sample.size_adjustedpertrait/5)                              
      
      
    iteration_val <- m
    
    #have a "for" loop, start it at 0, and end it at 4 for 5-fold CV
    for(j in 0:4){
      cross_val_val <- j
      pred <- random.sample_adjustedpertrait[((increment*j)+1):min(((increment*j)+increment) , sample.size_adjustedpertrait)]                    #e.g. For Trait 1 (200 lines total) you'll have 160 individuals in training population
      train <- random.sample_adjustedpertrait[-(((increment*j)+1):min(((increment*j)+increment) , sample.size_adjustedpertrait))]                #e.g. For Trait 1 (200 lines total) you'll have 40 individuals in prediction set
      
      ## For loop to get rid of the indices corresponding to missing data
      count <- 0
      for(k in indices_eliminated){
        if(count > 0) k <- k-count #if you need to subtract more than 1 for more than two lines missing, replace 1 with count
        
        if(any(train > k)) {
          train.temp.1 <- train[which(train > k)]
          train.temp.2 <- train[-which(train > k)] 
          try(train.temp.1 <- train.temp.1 - 1)
          train <- c(train.temp.1, train.temp.2)
        }
        
        if(any(pred > k)) {
          pred.temp.1 <- pred[which(pred > k)]
          pred.temp.2 <- pred[-which(pred > k)] 
          try(pred.temp.1 <- pred.temp.1 - 1)
          pred <- c(pred.temp.1, pred.temp.2)
        }
        
        count <- count+1
      } #End for-loop of k (for index correction)
      
      Geno.train <- Geno[train,]
      Pheno.train <- Y[train,]
      Geno.pred <- Geno[pred,]  
      Pheno.pred <- Y[pred,]
      
      ##################
      
      ## Next 3 lines used for GLMNet. alpha=1 is lasso penalty, alpha=0 is ridge regression penalty, middle ground is elastic net. two commented cv.fit functions are alternatives for model-fitting.
      if(method == 'GLMNet'){
        #cv.fit <- cv.glmnet(Geno.train, Pheno.train, nlambda=250, nfold=20)
        #cv.fit <- cv.glmnet(Geno.train, Pheno.train, lambda=lambda_seq)
        cv.fit <- cv.glmnet(Geno.train, Pheno.train, alpha = alpha, nlambda=250)
        ans <- predict(cv.fit, Geno.pred, s="lambda.min")
        r.gy <- c(r.gy, cor(ans, Pheno.pred))
        method.specs = paste(method,"_",alpha,sep='')
      }
      
      ## Next 2 lines used for rrBLUP
      if(method == "RR-BLUP"){
        ans <- kinship.BLUP(y=Pheno.train,G.train=Geno.train,G.pred=Geno.pred,K.method="RR")
        acc_temp <- cor(ans$g.pred,Pheno.pred)
        temprow <- c(current_trait_name, iteration_val, cross_val_val, acc_temp)
        r.gy <- rbind(r.gy, temprow)
        method.specs = method
      }
      
      #contribution by Akiyoshi Kiode 2022
      if(method == "RKHS" || method == 'RKHS') {
        # Need to check how numericalization 0-2 or -1 to 1 might mess with this implementation
        # gaussian kernel, 1 bandwidth parameter
        X <- Geno
        Y_RKHS <- as.matrix(Y)
        tst<- pred
        yNA<- as.numeric(Y_RKHS)
        yNA[tst] <- rep(NA, length(yNA[tst]))
        
        ###   X<-Geno.train #data called earlier
        ###   Y_RKHS<-Pheno.train
        
        ### DISTANCE MATRIX #############################
        D<-as.matrix(dist(X,method="euclidean"))^2
        D<-D/mean(D)
        h<-1
        ###   yNA<-Y_RKHS
        K<-exp(-h*D)
        ### MODEL FITTING #################################
        
        # y must be a "numeric,n"
        # ETA is a list
        yNA<-(as.numeric(yNA)) # fixed y[train] dim error
        
        fm=BGLR(y=yNA,ETA=list(list(K=K,model='RKHS')),nIter=6000,burnIn=1000) #nIter was 6000 fm = fitmodel?
        print('reached RKHSCor')
        RKHSCor <- cor(fm$yHat[tst],Pheno.pred)
        temprow <- c(current_trait_name, iteration_val, cross_val_val, RKHSCor)
        r.gy <- c(r.gy, temprow)
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
  colnames(r.gy.output) <- c("trait","iteration", "fold", "accuracy") 
  
  write.table(r.gy.output, paste(path.for.results,"/GP.Results_Maize.CAM_AF12_",method.specs,"-",current_trait_name,"_20-iterations.txt",sep=''), quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  
  
  return(r.gy)
  
} #End function "GP_master_script"
#^^Previously: End for-loop of i (for traits)

########################################################################
#
#   Call the new function with future_lapply
#
########################################################################

looper = 1:n.traits
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
colnames(r.gy.output) <- c("trait","iteration", "fold", "accuracy") 
write.table(r.gy.output, paste(path.for.results,"Results-GP-RR-AF12-5-20.txt",sep=''), quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
