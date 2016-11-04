 # Viva 450K

# Preprocessing and normalization of methylation data objects

##############################

library(Hmisc)
library(methylumi)
library(wateRmelon)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(matrixStats)
library(data.table)

# showplots <- F
# saveplots <- F

##############################
# Background correction
##############################

load("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/raw_idat_V14_official.RData")

# if we want to work on beta scale (rather than M)
# Triche paper suggests goob is better but not yet working

# go with "noob": normal-exponential out of band background correction
system.time(dat.bg <- methylumi.bgcorr(raw.idat.official, method = "noob", parallel = F))


dat.bg #info

# save this data out - takes a while!
# if(0){
  save(dat.bg, file = file.path("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/", paste0("datbg_official_", Sys.Date(), ".Rdata")))
# }


##############################
# Dye bias and normalization
##############################

# use dye bias function from lumi
#system.time(dat.db <- adjColorBias.quantile(dat.bg))
#gc()
# takes 136 minutes
# note that the resulting betas are identical
# not sure this is working... expects M values?

# use dye bias function from methylumi
# normalizeMethyLumiSet()
system.time(dat.bg.db <- normalizeMethyLumiSet(dat.bg))
# Normalizing via Illumina controls...
# Using sample number 400 as reference level...
# took 1.5 minutes

betas(dat.bg.db[1:3,1:3])
betas(dat.bg.db[1:3,1:3])

# but what is this really doing?
# methylumi:::normalizeViaControls
# normalizeMethyLumiSet

# normalization schemes for comparison
#system.time(dat.db.lumibg <- normalizeMethyLumiSet(dat.bg.lumi))
#system.time(dat.db.nobg <- normalizeMethyLumiSet(raw.pf))

# replot
# if(showplots){
#   boxplotColorBias(dat.bg[,1:30])
#   boxplotColorBias(dat.db[,1:30])
#   plotColorBias1D(dat.bg[,1:30])
#   plotColorBias1D(dat.db[,1:30])
#   plotColorBias2D(dat.bg[,1:30])
#   plotColorBias2D(dat.db[,1:30])
#   # not really appropriate for bimodal distributions
#   boxplot(log(unmethylated(raw.pf[,1:16])), subset = T)
# }


# wateRmelon metrics of agreement
# SNP probe standard error for each genotype (lower is better)
mean(genki(betas(raw.idat.official))) # 4.028e-05
mean(genki(betas(dat.bg))) # 4.531e-05

# variation at imprinting DMRs between samples
dmrse(betas(raw.idat.official))
dmrse(betas(dat.bg))
# just between probes
dmrse_col(betas(raw.idat.official))
dmrse_col(betas(dat.bg))
# just between samples
dmrse_row(betas(raw.idat.official))
dmrse_row(betas(dat.bg))
dmrse_row(betas(dat.bg.db))

rm(dat.bg)
#
## take out data no longer needed (bead # and variances - wait on stripOOB for now)
print(object.size(raw.idat.official), units = "MB")
raw.pf <- stripBeadNs(raw.idat.official)
raw.pf <- stripBeadSDs(raw.idat.official)
gc()


save(dat.bg.db, file = file.path("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/", paste0("datbgdb_official_", Sys.Date(), ".Rdata")))

