#BiocManager::install('MSnbase')
library(MSnbase)

# DATA STRUCTURE ------------------
# read example mzXML raw data
file <- dir(system.file(package = "MSnbase", dir = "extdata"),
            full.names = TRUE, pattern = "mzXML$")
rawdata <- readMSData(file, msLevel = 2, verbose = T)

# peak list in Mgf format can be read by readMgfData()

# MS data can be exported. here as mzML format. copy=T means metadata is also written
writeMSData(rawdata, file = paste0(tempfile(), ".mzML"), copy = TRUE)

# MSnExp objects, that stores all the spectra of an experiment, as defined by one or multiple raw data files
itraqdata

# Biobase::fData() return a dataframe with features as rows, variables as columns from a eSet-based object
head(fData(itraqdata))

# extract single spectrum data from MSnExp object
sp <- itraqdata[["X1"]]
sp
peaksCount(sp)
head(peaksCount(itraqdata))
rtime(sp) # retention time

# label-based proteome MS use reporter ion (iTRAQ or TMT etc.) to mark sample from different groups
# these common isobaric tags instances are provided in this package
iTRAQ4 # their MZ value (mass/charge) are recorded
TMT16

# BiocManager::install('msdata')
f <- c(system.file("microtofq/MM14.mzML", package = "msdata"))
mtof <- readMSData(f, mode = "onDisk")

# extract the (MS level 1) chromatogram.
# Without providing an m/z and a retention time range the function returns the total ion chromatogram (TIC) for each file within
mtof_tic <- chromatogram(mtof)
mtof_tic

# The chromatogram method returns a Chromatograms object (note the s)
# which holds multiple Chromatogram objects and arranges them in a two-dimensional grid with columns representing files/samples of the MSnExp /OnDiskMSnExp object and rows m/z-retention time ranges
mtof_tic[1, 1] # a single chromatogram object
head(intensity(mtof_tic[1, 1]))
head(rtime(mtof_tic[1, 1]))
mz(mtof_tic[1, 1]) # mz ratio

# extract the base peak chromatogram (the largest peak along the m/z dimension for each retention time/spectrum) we set the aggregationFun argument to "max"
mtof_bpc <- chromatogram(mtof, aggregationFun = "max")

# PLOTTING RAW DATA ----------------
# The MSmap class can be used to isolate specific slices of interest from a complete MS acquisition by specifying m/z and retention time ranges. 
# Below we first download a raw data file from the PRIDE repository 
# and create an MSmap containing all the MS1 spectra between acquired between 30 and 35 minutes and peaks between 521 and 523 m/z
# BiocManager::install('rpx')
library("rpx")
px1 <- PXDataset("PXD000001")
mzf <- pxget(px1, 7)

## reads the data
ms <- openMSfile(mzf)
hd <- header(ms)

## a set of spectra of interest: MS1 spectra eluted
## between 30 and 35 minutes retention time
ms1 <- which(hd$msLevel == 1)
rtsel <- hd$retentionTime[ms1] / 60 > 30 &
  hd$retentionTime[ms1] / 60 < 35

## the map
M <- MSmap(ms, ms1[rtsel], 521, 523, .005, hd, zeroIsNA = TRUE)

# The M map object can be rendered as a heatmap with plot()
plot(M, aspect = 1, allTicks = FALSE)

# One can also render the data in 3 dimension with the plot3D()
plot3D(M)

# plot MS spectra
plot(sp, reporters = iTRAQ4, full = TRUE)

# extract BSA data
sel <- fData(itraqdata)$ProteinAccession == "BSA"
bsa <- itraqdata[sel]
bsa
plot(bsa, reporters = iTRAQ4, full = FALSE) + theme_gray(8)
