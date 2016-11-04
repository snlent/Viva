#####################################################################
###################### Viva 450k pipeline ###########################
##### Adapted from Allan Just's code by Sam Lent (lent@bu.edu) ######
#####################################################################


datadir<-'/apps/sasr/data/idats/IDATs_final'
path.samplesheets<-'/apps/sasr/data/idats/samplesheets'
#Should have 3 sample sheet files: VIVA_sample_mapping_with_age.xlsx, BWH-Demeo_2014_SampleSummary.xlsx, IDAT_Trnslation_Key_B.xlsx
lib.path<-'/home/slent/R/x86_64-redhat-linux-gnu-library/3.2/' #Where R packages are installed

importnotload<-FALSE
showplots<-FALSE
saveplots<-TRUE

library(Hmisc,lib.loc=lib.path)
library(methylumi)
library(wateRmelon)
library(ggplot2)
library(data.table)
library(tools)
library(minfi)
library(matrixStats)
library(RColorBrewer)
library(plyr)
library(grid) # for viewport on subplot
library(sas7bdat)
library(XLConnect)


##### Import idats and create a data.table (dt) with information on each idat #####
# dt variables:
# path: idat path (not full path, path from datadir)
# idat: idat filename
# barcode: idat filename minus _Red.idat or _Grn.idat
# md5: MD5 checksums, unique files have unique MD5 checksums
# md5index: Number of times this md5 checksum occurs in all idat files (finds duplicate chips)
# nbarcode: Number of times this barcode occurs (non-duplicate samples will have a 2 here, one for red and one for green)
# mtime: File modification time
# subfolder: flag indicating whether this idat is from a subfolder
# folder: folder containing idat [note: path=paste0(folder,idat)]
# color: Grn or Red
#
# dt saved as viva_idat_file_table_20151217.csv

list.files(datadir)

allidat <- list.files(path = file.path(datadir), recursive = T, pattern = "idat$")
idattrans <- data.table(read_excel("/apps/sasr/data/idats/samplesheets/IDAT_Translation_Key_B.xlsx", 
                                   sheet = 1),check.names=TRUE)
dt <- data.table(path = allidat)
dt[, idat := sub(".*/", "", allidat)]
dt[, barcode := substring(text = idat, 1, 17)]
dt[, md5 := md5sum(file.path(datadir, path))]
if(dt[, length(unique(md5))]==dim(dt)[1]) {print('MD5 checksums all unique')}
if(dt[, length(unique(md5))]==dim(dt)[1]) {print(paste0(dim(dt)[1]-dt[, length(unique(md5))],' MD5 checksums not unique'))} #creates md5 checksums
dt[, md5index := .N, by = md5] #creates md5index
dt[md5index > 1,] # pull out non-unique files

# it looks like for two copies of chip 9721367038
# the Red and Grn idat files are identical


describe(dt[, .N, by = barcode])
dt[, .N]
dt[, length(unique(barcode))]
dt[, nbarcode := .N, by = barcode]
dt[nbarcode == 5]

# extract file modification time
dt[, mtime := file.info(file.path(datadir, dt$path))$mtime]
nrow(dt)
dt[, length(unique(mtime))]
dt[, describe(mtime)]

# are redone barcodes always within a subfolder? 
describe(dt[nbarcode > 2, nchar(path)])
describe(dt[nbarcode == 2, nchar(path)])
# create an index for whether they are in a subfolder
dt[, subfolder := nchar(path) > 37]

# can I step through, import each folder (a single chip scan), run some quick descriptives?
# first - for the path
dt[, folder := substr(path, 1, nchar(path) - 26)]
dt[, color := substr(idat, nchar(idat) - 7, nchar(idat) - 5)]
folders <- unique(dt[, folder])
# and we store it only for Green (to not have duplicates)
Sys.time(); length(folders)
for(i in 1:length(folders)){
  # allowing that some files might fail
  methyfolder <- try(methylumIDAT(barcodes = dt[color == "Grn" & folder == folders[i], barcode], 
               idatPath = paste0(datadir, '/',folders[i])))
  if(class(methyfolder)[1] != "try-error"){
    # store out median methylated & unmethylated intensity and # of failed probes
    dt[color == "Grn" & folder == folders[i], `:=`(mMed = log2(colMedians(methylated(methyfolder))),
                                                 uMed = log2(colMedians(unmethylated(methyfolder))), 
                                                 nbadprobes = colSums(pvals(methyfolder) > 0.05))] ***Figure out why colSums gives evidence of bad probes***
  } else print(paste("error thrown for folder", folders[i]))
  print(i)
}
describe(dt[color == "Grn"])
# order files then number the versions
setkey(dt, barcode, folder)
dt[grep("9422492047", barcode),]
dt[color == "Grn", subfolderversion := 1:.N, by = c("barcode")]
# make a plot
ggplot(dt[color == "Grn"], aes(subfolderversion, nbadprobes, group = barcode)) + geom_line()
# need to remove some obvious outliers - a terrible chip
dt[which.max(nbadprobes)]
dt[, nbadprobesub := nbadprobes]
dt[folder == "9721367038/119537129/", nbadprobesub := NA]
# plot of how bad each sample is, joining copies in a line
ggplot(dt[color == "Grn" & nbarcode > 2], aes(subfolderversion, nbadprobesub, group = barcode)) + geom_line()
# check that the number of failed probes isn't the same for too many samples
dt[order(nbadprobesub, decreasing = T)][subfolderversion == 2 & !is.na(nbadprobesub)][1:10]
# now compute with intensity
dt[color == "Grn", meanintensity := (mMed + uMed)/2]
# some bad values - look at these
dt[meanintensity == -Inf]
# these are from the chip with identical green and red idat files
dt[meanintensity == -Inf, nbadprobesub := NA]
dt[meanintensity == -Inf, meanintensity := NA]
dt[folder == "9721367038/119537129/", meanintensity := NA]
# review - how many samples were rescanned?
dt[nbarcode > 2, length(unique(barcode))]
# how often is the second version better than the first?
# by number of failed probes
dt[color == "Grn" & nbarcode > 2, bestversion := .SD[which.min(nbadprobesub), subfolderversion], by = "barcode"]
dt[color == "Grn" & nbarcode > 2 & subfolder == F, list(nbarcode, bestversion), by = "barcode"][, addmargins(table(nbarcode, bestversion))]
# now by intensity (here higher is better)
ggplot(dt[color == "Grn" & nbarcode > 2], aes(subfolderversion, meanintensity, group = barcode)) + geom_line()
dt[color == "Grn" & nbarcode > 2, bestversionintensity := .SD[which.max(meanintensity), subfolderversion], by = "barcode"]
dt[color == "Grn" & nbarcode > 2 & subfolder == F, list(nbarcode, bestversionintensity), by = "barcode"][, addmargins(table(nbarcode, bestversionintensity))]
# when is the best version in terms of failed probes not the same as for intensity?
dt[bestversion != bestversionintensity, list(folder, 
                                             nbadprobesub, meanintensity, 
                                             bestversion, bestversionintensity),by=barcode][, length(unique(barcode))]
dt[color == "Grn" & nbarcode > 2 & subfolder == F, list(nbarcode, bestversionintensity), by = "barcode"][, addmargins(table(nbarcode, bestversionintensity))]

# how do the first subfolder versions compare to the non-rescanned?
newbest <- dt[color == "Grn" & ((nbarcode == 2 & subfolder == F) | (nbarcode > 2 & subfolder == T & subfolderversion == 2))]

ggplot(dt[color == "Grn"], )

# merge in which samples failed
newbest[, c("keep", "flaggedbad") := NULL]
newbest <- merge(newbest, idattrans[, list(barcode = samplename, keep, flaggedbad)], by = "barcode", all.x = T)
ggplot(newbest[keep == T & is.na(flaggedbad)], aes(subfolder, y = nbadprobesub)) + 
  geom_point(position = position_jitter(height = 0, width = 0.3)) + 
  geom_boxplot(color = "red", outlier.size = 0, alpha = 0.3) + 
  ylab("Number of bad probes per sample\n(detection p-val > 0.05)") + 
  theme_bw(14) + theme(axis.title.x = element_blank())
with(newbest[keep == T & is.na(flaggedbad), ], t.test(nbadprobesub~subfolder))

# export this for Dawn
write.csv(dt, file = "viva_idat_file_table_20141010.csv", row.names = F)
