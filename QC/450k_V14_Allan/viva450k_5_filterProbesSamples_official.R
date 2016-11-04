####################################################
# Call in required libraries
###################################################
 
library(data.table)
library(Hmisc)
library(sas7bdat)
library(caret) # to check for rank deficiency
library(MASS) # rlm function for robust linear regression
library(sandwich) #Huber’s estimation of the standard error
library(lmtest) # to use coeftest
library(parallel) # to use multicore approach - part of base R
library(methylumi)
library(matrixStats)


# Load in metylumi dataset and fDat (for filtering out bad probes)
##########################################################################
load("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/datbgdbbmiq_official_CB_2015-12-17.Rdata")
load("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/fDat_official.RData")
load("/n/home03/goa306/Annotation_compiled.Rdata")


intersect = intersect(rownames(dat.bg.db.cord.bmiq), Annot$IlmnID)
length(intersect)

# First, filter out samples that need to be dropped
#########################################################################

dim(dat.bg.db.cord.bmiq)

rownames(dat.bg.db.cord.bmiq)[1:10]

class(rownames(dat.bg.db.cord.bmiq))


mdata <- dat.bg.db.cord.bmiq[,pData(dat.bg.db.cord.bmiq)$lowQual==T & pData(dat.bg.db.cord.bmiq)$filter.gender==T & 
                               pData(dat.bg.db.cord.bmiq)$filter.geno==T & pData(dat.bg.db.cord.bmiq)$repUse==T]

dim(mdata)
rm(dat.bg.db.cord.bmiq)


# Next, filter out at probe level
###################################################

# 1. Filter out bad probes

class(fDat)

head(fDat)

rownames(fDat)[1:10]

class(rownames(fDat))

dim(fDat)

summary(fDat$perc.detected < 99)

summary(fDat$perc.detected < 95)


fDatselect <- fDat[fDat$perc.detected > 95,]
dim(fDatselect)

mdata <- mdata[rownames(mdata) %in% rownames(fDatselect),]
dim(mdata )

# this part added in interactively!! - after running this script as a batch job already
# Save mdata at this point - for cell mixture generation!

# in addition, none of the CpGs at this point had detection p-val > 0.05
# checked with:
# fDatselect <- fDatselect[fDatselect$median <> 0.05,] ..and there were non.. checked several different ways.
save(mdata, file = paste0("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/", "datbgdbbmiq_official_CB_nobadprobes_", Sys.Date(), ".Rdata"))


intersect = intersect(rownames(mdata), Annot$IlmnID)
length(intersect)
# this means all CpGs of mdata are now covered by Annot

# 2. restrict to cpg probes

Annot$probetype <- substring(Annot$IlmnID,1,2)
table(Annot$probetype)

Annotfinal <- Annot[grepl("cg", Annot$IlmnID),]
dim(Annotfinal)

# 3.  Drop X/Y (Allosomal) probes
table(Annotfinal$chr)
Annotfinal <- Annot[!Annot$chr %in% c("chrX","chrY"),]
table(Annotfinal$chr)
dim(Annotfinal)

# 4. Drop SNP probes
  # 1st alternative way - based on Illumina450ProbeVariants.db package
  # THIS METHOD DROPS A LOT MORE - I TESTED IT
  # USE OTHER METHOD BELOW

# library(Illumina450ProbeVariants.db)
# data(probe.450K.VCs.af)
# table(probe.450K.VCs.af$probe10VC.all>0) # 81519!
# SNP10 <- probe.450K.VCs.af[probe.450K.VCs.af$probe10VC.all>0,]
# dim(SNP10)
# class(SNP10)
# rownames(SNP10)[1:5]
# class(rownames(SNP10))
# 
# mdata1 <- mdata[!rownames(mdata) %in% rownames(SNP10),]
# dim(mdata1)


  # 2nd alternative way - based on own annotation - created from Package IlluminaHumanMethylation450kanno.ilmn12.hg19
# Removing SNP associated probes with MAF>5%
Annotfinal$Probe_r_exclude <- rep(0,dim(Annotfinal)[1])
sum((Annotfinal$Probe_maf>=0.05),na.rm=T) #1036
Annotfinal$Probe_r_exclude[Annotfinal$Probe_maf>=0.05] <-1
table(Annotfinal$Probe_r_exclude)
Annotfinal = Annotfinal[Annotfinal$Probe_r_exclude==0,] #select MAF<5%
dim(Annotfinal) 

# SNPs at SBE removing remofe if MAF>=5%
Annotfinal$SNPatSBE <- rep(0,dim(Annotfinal)[1])
sum((Annotfinal$SBE_maf>=0.05),na.rm=T)
Annotfinal$SNPatSBE[Annotfinal$SBE_maf>=0.05] <-1
Annotfinal = Annotfinal[Annotfinal$SNPatSBE==0,] #select non polymorphic @ SBE
dim(Annotfinal) 

# SNPs in CpGs with MAF  removing MAF>=5%
Annotfinal$SNPatCPG <- rep(0,dim(Annotfinal)[1])
sum((Annotfinal$CpG_maf>=0.05),na.rm=T) 
Annotfinal$SNPatCPG[Annotfinal$CpG_maf>=0.05] <-1
Annotfinal = Annotfinal[Annotfinal$SNPatCPG==0,]
dim(Annotfinal) 

mdata <- mdata[rownames(mdata) %in% Annotfinal$IlmnID,]
dim(mdata)


# 4. remove x-reactives
#########################################################
xreactives<-read.csv(file="/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/NonSpecProbes.csv")
class(xreactives)
# xreactives <- data.table(xreactives)
dim(xreactives)
head(xreactives)

class(xreactives$TargetID)
xreactives$TargetID[1:10]
xreactives$TargetID <- as.character(xreactives$TargetID)

mdata <- mdata[!rownames(mdata) %in% xreactives$TargetID,]
dim(mdata)


save(mdata, file = paste0("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/", "datbgdbbmiq_official_CB_clean_", Sys.Date(), ".Rdata"))






