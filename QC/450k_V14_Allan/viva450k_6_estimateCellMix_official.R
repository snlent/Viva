# Estiamte cell mixture
# Code adapted from Andres - doing cell mixture estimation on cleaned betas as opposed to using minfi on the RGset (uncleaned)

library(methylumi)

load("/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/intermediate/datbgdbbmiq_official_CB_nobadprobes_2015-12-18.Rdata")
dim(mdata)

mbetas <- betas(mdata)
rm(mdata)

dim(mbetas)
mbetas[1:5,1:5]


require("minfi")
require("FlowSorted.Blood.450k")

## load betas from cord-blood
coefs = FlowSorted.Blood.450k.JaffeModelPars # cell type means at discriminatory probes
dim(coefs)

#Intersect
intersect = intersect(rownames(mbetas), rownames(coefs))
length(intersect)

coefs = coefs[intersect,]
dim(coefs)

cellProps = minfi:::projectCellType(mbetas[rownames(coefs),], coefs) # project
cellProps = as.data.frame(cellProps)


# Print the result
head(round(cellProps,3))


cellProps <- data.frame(cellProps)

class(rownames(cellProps)) #character
head(cellProps)

cellProps$samplename <- rownames(cellProps) 

class(cellProps$samplename) # "character"
head(cellProps)

# now save without rownames..

write.csv(cellProps, file = "/net/rcnfs05/srv/export/baccarelli_lab/share_root/viva_450k/data/output/estimatedcellcounts_vivacords_official.csv", row.names = F)
