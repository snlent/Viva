# Dec 2015 
# Generate a feature data with percent detection p-vals - with same variable as that in V14 data from Dawn
# This is largely adapted from Jarrett's code - I've followed in his methods for generating an fDat

library(minfi)
load("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/raw_RGset_V14_official.RData")

# Obtain detection p-values
detectionPVals <- detectionP(raw.RGset.official)

# extract the matrix of raw beta values
raw.RGset.official.raw <- preprocessRaw(raw.RGset.official)
betas.raw <- getBeta(raw.RGset.official.raw, type = "Illumina")


detected <- (detectionPVals<0.05)

# print(table(detected)) # this takes way too long

# find the percentage of high confidence probes at each location
detected.probes <- apply(detected,1,sum)/dim(detected)[2]*100

print(table(detected.probes<75))
print(table(detected.probes<99))
print(table(detected.probes<95))


# examine the median p-value for each location
detected.median <- apply(detectionPVals,1,median)
# how many have a median value greater than 0.05
print(table(detected.median>0.05))

fDat <- data.frame(site=rownames(raw.RGset.official.raw), perc.detected=detected.probes
                                       ,median=detected.median
                                       ,stringsAsFactors=F)

rownames(fDat) <- rownames(raw.RGset.official.raw)
dim(fDat)
head(fDat)

save(fDat, file = file.path("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/fDat_official.RData"))
