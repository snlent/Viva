# Viva 450K
# Generate official raw dataset-starting from iDATs- in coordination with official V14 data from Dawn and Jarrett
# clean failed samples
##############################

library(Hmisc)
library(methylumi)
library(wateRmelon)
library(ggplot2)
library(matrixStats)
library(sas7bdat)
# library(XLConnect)  CURRENTLY GIVING ERRORS
library(data.table)
library(readxl) # USING THIS INSTEAD OF XLConnect


###########################################################
# CALL IN IDAT TRANSLATION KEY SHEET
############################################################

idattrans <- data.table(read_excel("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/samplesheets/IDAT_Translation_Key_B.xlsx", 
                                   sheet = 1),check.names=TRUE)

dim(idattrans) 

# pull out a the unrelated subtable of chip and batch from idat translation key
beadchipbatch <- idattrans[, list(BeadChip = BeadChip.1, Batch)]
beadchipbatch <- beadchipbatch[!is.na(BeadChip)]
# drop some extra columns
idattrans[, c("Col9", "Col10", "Col11", "BeadChip.1", "Batch") := NULL]
head(idattrans)

# pull out sampleid
idattrans[, sampleid := sub("WG.*_", "", Sample.1)]
idattrans[, samplename := paste0(BeadChip, "_", Position)]

head(idattrans)

table(duplicated(idattrans$samplename)) # no duplicates 


######### FAILED SAMPLES ##################

# Create a flag for 'failed samples'
failed <- grep("Failed",idattrans$Comment)
print(table(idattrans$Comment[failed]))

idattrans[, failed := F]
idattrans[grepl("Fail", Comment), failed := T]
table(idattrans$failed)

############## CONTROL SAMPLES ######################

array.Ctrl <- idattrans[sampleid == "Ctrl"]
dim(array.Ctrl)

# create a flag for control samples
idattrans[, controlSamps := F]
idattrans[sampleid == "Ctrl", controlSamps := T]
table(idattrans$controlSamps)

# make the Ctrl samples have a unique sampleid
# (This doesn't have any immediate use in current script - but may be useful as code for later work on controls..)
idattrans[sampleid == "Ctrl", sampleid := paste0("Ctrl", 1:.N)]
head(idattrans[grepl("Ctrl",sampleid)])

####################################################
# CALL IN SHEET TO EXTRACT VIVA ALIAS
#########################################################

# mapping of SAMPLE e.g. "S-001259392" to Viva alias "122972"
vivaalias <- data.table(read_excel("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/samplesheets/VIVA_sample_mapping_with_age.xlsx", 
                                   sheet = 1),check.names=TRUE)

dim(vivaalias)
head(vivaalias)

vivaalias[, alias := as.numeric(gsub("c|C|-", "", RECRUITMENT.ID))]
head(vivaalias)

# check alias values are all in the fullphenodata that Dawn sent
setnames(vivaalias, "SAMPLE", "sampleid")

describe(vivaalias) # only 728 unique for 1125 unique sampleids in a list of 1161 values

vivaalias[duplicated(vivaalias)] # our repeats are c116996 and c107473

# samples that were run multiple times --> 32 + 16 =48
vivaalias[alias %in% c(116996, 107473)]

# is each sampleid occur only once - not always!
vivaalias[, .N, by=sampleid][, describe(N)]
# well, apart from the two aliases above it is:
vivaalias[!alias %in% c(116996, 107473), .N, by = sampleid][, describe(N)]
# does each sampleid map to a unique alias? Yes!
vivaalias[, list(nalias = length(unique(alias))), by=sampleid][, describe(nalias)]
# each non-qc alias can have 1-3 sampleids
vivaalias[, .N, by=alias][, describe(N)]

vivaalias[, .N, by=sampleid][, describe(N)]

# clean up sample type/stage
vivaalias[, stage := factor(COLLECTION, 
                            levels = c("Proband Delivery", "Proband Age 3", "Proband Age 7"),
                            labels = c("CB", "Y3", "Y7"))]
# make a var for qc samples
vivaalias[alias == 116996 & stage == "CB", qcid := "MR32"]
vivaalias[alias == 107473 & stage == "CB", qcid := "MR16"]



#################################################
# MERGE IN ALIAS WITH IDAT TRANSLATION KEY TO HANDLE MULTIPLE REPEATS
####################################################################

idattrans[, c("alias", "stage", "qcid") := NULL] # take it out first in case repeating
idattrans <- merge(idattrans, unique(vivaalias[,list(sampleid, alias, stage, qcid)]), 
                   by = "sampleid", all.x = T)

### WE WANT TO KEEP 1 FROM EACH MULTIPLE REPEAT - ALLAN KEEPS THE SAME ONES THAT JARRETT DID!! ###
# > pDat[sampleNames=='S-000670341','repUse'] <- TRUE
# > pDat[sampleNames=='S-000621296','repUse'] <- TRUE

# need an indicator for which samples are solely qc
idattrans[, isqc := F]
idattrans[alias %in% c(116996, 107473) & stage == "CB", isqc := T]

# but we get back one of each - we need to pick
# but which are the designated qc samps to use? (which of the 16 or 32 are on same chip as other samples?)

# which chips are the non Cord Blood samples on ?
idattrans[alias == 116996 & stage != "CB", list(sampleid, BeadChip, stage)]

# which samples from that person are on those chips?
idattrans[alias == 116996 & BeadChip %in% .Last.value[, BeadChip], list(sampleid, BeadChip, stage)]

# how many of each sampleid were plated? (especially from the previous command)
idattrans[alias == 116996 & stage == "CB", .N, by=sampleid]

# so we will designate the singleton that is on the same chip as other samples - just like everyone else
idattrans[alias == 116996 & stage == "CB" & sampleid == "S-000621296", isqc := F]

# again for the next Multiple-Repeat sample
idattrans[alias == 107473 & stage != "CB", list(sampleid, BeadChip, stage)]

idattrans[alias == 107473 & BeadChip %in% .Last.value[, BeadChip], list(sampleid, BeadChip, stage)]

idattrans[alias == 107473 & stage == "CB", .N, by=sampleid]

idattrans[alias == 107473 & stage == "CB" & sampleid == "S-000670341", isqc := F]



################################################################
# SEE WHAT NUMBERS LOOK LIKE WHEN YOU SEQUENTIALLY DROP ACCORDING TO FLAGS
##############################################################################

dim(idattrans) 

# dropping failed samples
idattrans1 <- idattrans[failed==F,]
dim(idattrans1) #so 226 dropped - as was the case in Jarrett's code

# dropping control samples
idattrans2 <- idattrans1[controlSamps==F,] # so 25 dropped!! This now matches array.Ctrl in Jarrett's code!
## array.Ctrl
## FALSE  TRUE 
##  1151    25

dim(idattrans2) # THIS MATCHES 1151 IN pDat from V14 FROM DAWN!!

print(idattrans2[is.na(Sample)]) # none - these all disappeared with removing 'failed samples'

head(idattrans2)
class(idattrans2) # [1] "data.table" "data.frame"

# NOW LOAD V14 VIVA DATA OFFICIALLY PROVIDED BY DAWN!!
load("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/VivaThousandQC_BMIQ_1449230271_V14.RData")
pDat <- data.table(pDat)
head(pDat)

# Both these have many of the same variables - make sure that they all correspond!
# create a sampleid var in pDat just like that in idattrans2
pDat[, sampleid := sub("WG.*_", "", Sample.1)]
pDat[, samplename := paste0(Slide, "_", Array)]

# also create a similar variable in idattrans2 from beadchip and position
idattrans2[, BeadChip_Position := paste0(BeadChip, "_", Position)]

dim(pDat) 
dim(idattrans2) 

class(pDat$sampleid) 
class(idattrans2$sampleid) 

class(pDat$samplename) 
class(idattrans2$samplename) 

idattrans2_check <- idattrans2[idattrans2$sampleid %in% pDat$sampleid,]
dim(idattrans2_check) 

idattrans2_check2 <- idattrans2[idattrans2$Sample %in% pDat$Sample,]
dim(idattrans2_check2) 

idattrans2_check3 <- idattrans2[idattrans2$samplename %in% pDat$samplename,]
dim(idattrans2_check3) 

idattrans2_check4 <- idattrans2[idattrans2$BeadChip_Position %in% pDat$samplename,]
dim(idattrans2_check4) 

# Does "questionable quality" from our idattrans2 match "lowQual" from dawn's pDat? - YES!
table(pDat$lowQual) 

idattrans2[, lowQual := F]
idattrans2[grepl("Questionable Quality", Comment), lowQual := T]
table(idattrans2$lowQual)

# does our "isqc" match Dawn's "repUse" - YES!
table(pDat$repUse) 
table(idattrans2$isqc)


# # # so everything matches, now use pDat to make a manifest to replace idattrans2 # # # 
length(unique(idattrans2$samplename)) 
describe(idattrans2$samplename)

length(unique(pDat$samplename)) 
describe(pDat$samplename)

# match idattrans2 to pDat - ensure it is fully aligned 
idattransfinal <- idattrans2[match(pDat$samplename, samplename),]

# Check thoroughly
all(idattransfinal$samplename==pDat$samplename) 
identical(idattransfinal$samplename,pDat$samplename) 
all(idattransfinal$BeadChip_Position==pDat$samplename) 
all(idattransfinal$sampleid==pDat$sampleid) 
all(idattransfinal$Sample==pDat$Sample) 

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# NOW YOU CAN USE Dawn's pDat INSTEAD OF 'IDATTRANS' (i.e.IDATTRANS2) TO MERGE IN AS A MANIFEST
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# But we also want to add information on multiple scans
# Call in master sheet created by Allan - with information on all folders/subfolders
files <- read.csv("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/viva_idat_file_table_20141010.csv")

dim(files) # [1] 3728   18
files <- data.table(files)
head(files)

# propagate subfolderversion to both green and red file for each barcode/folder combo (each scan)
# GA: (b/c Allan initially filled in info form a lot of summary variables just for the 'Grn' files')
files[, subfolderversion := .SD[color == "Grn", subfolderversion], by = c("barcode", "folder")]


# We want to use parent-level(1st-level) folder IDATs - all samples coming from the same scanner!
dim(files[subfolder == F]) 
dim(files[subfolderversion == 1,]) 
# so this corresponds with subfolder == F
# this means that if subfolder == F, then subfolderversion == 1 - these are the parent-level IDATs

# Subset to parent-level folders
filesParent <- files[subfolder == F,]
dim(filesParent) 

# So how many where nbarcode > 2 (i.e. scanned more than once)?
dim(filesParent[nbarcode > 2])

table(filesParent$nbarcode, exclude=NULL)


# make a variable to indicate which beadchips were scanned more than once vs. only once
filesParent[,multipleScans := 0] 
filesParent[nbarcode > 2, multipleScans := 1] 
table(filesParent$multipleScans)


# You also want to make a variable to indicate number of scans - based on 'nbarcode'
filesParent[,NumOfScans:=nbarcode/2] 

table(filesParent$nbarcode, exclude=NULL)
table(filesParent$NumOfScans, exclude=NULL)
# worked perfectly

# now divide # of beadchips by 2, because you have 2 of everything right now (1 for green, 1 for red)
filesParentHalf<- filesParent[color == "Grn"]
dim(filesParentHalf) # [1] 1402   19
table(filesParentHalf$multipleScans)
# 384 matches the number indicated in Allan's email: "Illumina rescanned 384 samples..."
table(filesParentHalf$NumOfScans)


# Now use barcode/samplename to merge scan info into PDat
dim(pDat)
dim(filesParentHalf)

# samplename in pDat = barcode
class(pDat$samplename)
class(filesParentHalf$barcode)

filesParentHalf$barcode <- as.character(filesParentHalf$barcode)

filesParentFinal <- filesParentHalf[filesParentHalf$barcode %in% pDat$samplename,]
dim(filesParentFinal)

filesParentFinal <- filesParentFinal[match(pDat$samplename, barcode),]

all(filesParentFinal$barcode==pDat$samplename) 
identical(filesParentFinal$barcode,pDat$samplename) 

filesParentFINAL <- filesParentFinal[,c("barcode","multipleScans","NumOfScans"),with=FALSE]
head(filesParentFINAL)

all(filesParentFINAL$barcode==pDat$samplename) 
identical(filesParentFINAL$barcode,pDat$samplename) 

pDatFINAL <- cbind(pDat,filesParentFINAL)

all(pDatFINAL$barcode==pDat$samplename) 
identical(pDatFINAL$barcode,pDat$samplename) 
identical(pDatFINAL$barcode,pDatFINAL$samplename) 
identical(filesParentFINAL$barcode,pDatFINAL$samplename) 

table(pDatFINAL$multipleScans)
table(filesParentFinal$multipleScans)
table(filesParentFINAL$multipleScans)

table(pDatFINAL$NumOfScans)
table(filesParentFinal$NumOfScans)
table(filesParentFINAL$NumOfScans)

# Rename pDatFINAL back to 'pDat' before merging it in as methylumi manifest
rm(pDat)
pDat <- pDatFINAL
dim(pDat)

# Call in the official iDATs obtained from Dawn and imported in Dec 2015

load("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/imported/raw_idat_all_official.RData")
load("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/imported/raw_RGset_all_official.RData")

dim(raw.idat.all)
dim(raw.RGset.all)

class(sampleNames(raw.idat.all))
class(sampleNames(raw.RGset.all))
class(pDat$samplename)

# Restrict imported iDAT data to what is in pDat from Dawn

raw.idat.official <- raw.idat.all[,sampleNames(raw.idat.all) %in% pDat$samplename]
dim(raw.idat.official)

raw.RGset.official <- raw.RGset.all[,sampleNames(raw.RGset.all) %in% pDat$samplename]
dim(raw.RGset.official)

rm(raw.idat.all)
rm(raw.RGset.all)


# prepare to merge in with datasets
##########################################################
pDat <- pDat[match(sampleNames(raw.idat.official), samplename),]
identical(sampleNames(raw.idat.official), pDat$samplename)
all(sampleNames(raw.idat.official)==pDat$samplename)

# bring them in as pData
manifest1 <- as.data.frame(pDat[match(sampleNames(raw.idat.official), pDat$samplename),])
row.names(manifest1) <- manifest1$samplename


# put them into the methylumi object
if(exists("raw.idat.official") & identical(sampleNames(raw.idat.official), manifest1$samplename)){
  pData(raw.idat.official) <- manifest1
} else { stop("REORDER YOUR DFCOVARS!!!")}

save(raw.idat.official, file = "/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/raw_idat_V14_official.RData")
rm(raw.idat.official)

######

pDat <- pDat[match(sampleNames(raw.RGset.official), samplename),]
identical(sampleNames(raw.RGset.official), pDat$samplename)
all(sampleNames(raw.RGset.official)==pDat$samplename)


manifest2 <- as.data.frame(pDat[match(sampleNames(raw.RGset.official), pDat$samplename), ])
row.names(manifest2) <- manifest2$samplename


# put them into the minfi object
if(exists("raw.RGset.official") & identical(sampleNames(raw.RGset.official), manifest2$samplename)){
  pData(raw.RGset.official) <- manifest2
} else { stop("REORDER YOUR DFCOVARS!!!")}

save(raw.RGset.official, file = "/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/raw_RGset_V14_official.RData")
rm(raw.RGset.official)
