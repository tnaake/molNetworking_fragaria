library(xcms)

## abacus
## negative
setwd("~/winhome/Documents/05_collaborations_small_projects/Jose_fragaria/mzML/neg")
xset_neg <- xcmsSet(file = "./", method = "centWave", ppm = 10, snthres = 20, 
                    peakwidth = c(5, 20), prefilter = c(3, 5000))
classes <- c(
    rep("156_G", 3), rep("156_R", 3), rep("156_W", 3), 
    rep("191_G", 3), rep("191_R", 3), rep("191_W", 3), 
    rep("196_G", 3), rep("196_R", 3), rep("196_W", 3),
    rep("282_G", 3), rep("282_R", 3), rep("282_W", 3),
    rep("591_G", 3), rep("591_R", 3), rep("591_W", 3),
    rep("595_G", 3), rep("595_R", 3), rep("595_W", 3),
    rep("660_G", 3), rep("660_R", 3), rep("660_W", 3),
    rep("Amiga_G", 3), rep("Amiga_R", 3), rep("Amiga_W", 3),
    rep("Benicia_G", 3), rep("Benicia_R", 3), rep("Benicia_W", 3),
    rep("Candonga_G", 2), rep("Candonga_R", 3), rep("Candonga_W", 3),
    rep("Fontanilla_G", 3), rep("Fontanilla_R", 3), rep("Fontanilla_W", 3),
    rep("Fuentepina_G", 3), rep("Fuentepina_R", 3), rep("Fuentepina_W", 3),
    rep("SantaClara_G", 3), rep("SantaClara_R", 3), rep("SantaClara_W", 3)
) ## 
sampclass(xset_neg) <- classes
xset2_neg <- group(xset_neg, method = "density", minfrac = 0.5, minsamp = 1, bw = 2, mzwid = 0.015)
xset3_neg <- retcor(xset2_neg, family = "s", plottype = "m", missing = 1, extra = 1, span = 1)
xset4_neg <- group(xset3_neg, method = "density", mzwid = 0.015, minfrac = 0.25, minsamp = 1, bw = 2)
xset5_neg <- fillPeaks(xset4_neg, method = "chrom")

save(xset_neg, xset2_neg, xset3_neg, xset4_neg, xset5_neg, file = "../../xcms_strawberry_neg.RData")

library(CAMERA)
an_neg <- xsAnnotate(xset5_neg)
anF_neg <- groupFWHM(an_neg, perfwhm = 0.6)
anI_neg <- findIsotopes(anF_neg, mzabs = 0.01)
anIC_neg <- groupCorr(anI_neg, cor_eic_th = 0.75, graphMethod = "lpc")
anFA_neg <- findAdducts(anIC_neg, polarity = "negative")
peaklist_neg <- getPeaklist(anFA_neg)

save(an_neg, anF_neg, anI_neg, anIC_neg, anFA_neg, peaklist_neg, file = "../../CAMERA_strawberry_neg.RData")

## positive
setwd("~/winhome/Documents/05_collaborations_small_projects/Jose_fragaria/mzML/pos")
xset_pos <- xcmsSet(file = "./", method = "centWave", ppm = 10, snthres = 20, 
                    peakwidth = c(5, 20), prefilter = c(3, 5000))
classes <- c(
    rep("156_G", 3), rep("156_R", 3), rep("156_W", 3), 
    rep("191_G", 3), rep("191_R", 3), rep("191_W", 3), 
    rep("196_G", 3), rep("196_R", 3), rep("196_W", 3),
    rep("282_G", 3), rep("282_R", 3), rep("282_W", 3),
    rep("591_G", 3), rep("591_R", 3), rep("591_W", 3),
    rep("595_G", 2), rep("595_R", 3), rep("595_W", 3),
    rep("660_G", 3), rep("660_R", 3), rep("660_W", 3),
    rep("Amiga_G", 3), rep("Amiga_R", 3), rep("Amiga_W", 3),
    rep("Benicia_G", 3), rep("Benicia_R", 3), rep("Benicia_W", 3),
    rep("Camarosa_G", 3), rep("Camarosa_R", 3), rep("Camarosa_W", 3),
    rep("Candonga_G", 3), rep("Candonga_R", 3), rep("Candonga_W", 3),
    rep("Fontanilla_G", 3), rep("Fontanilla_R", 3), rep("Fontanilla_W", 3),
    rep("Fuentepina_G", 3), rep("Fuentepina_R", 3), rep("Fuentepina_W", 3),
    rep("SantaClara_G", 3), rep("SantaClara_R", 3), rep("SantaClara_W", 3)
) ## 
sampclass(xset_pos) <- classes
xset2_pos <- group(xset_pos, method = "density", minfrac = 0.5, minsamp = 1, bw = 2, mzwid = 0.015)
xset3_pos <- retcor(xset2_pos, family = "s", plottype = "m", missing = 1, extra = 1, span = 1)
xset4_pos <- group(xset3_pos, method = "density", mzwid = 0.015, minfrac = 0.5, minsamp = 1, bw = 2)
xset5_pos <- fillPeaks(xset4_pos, method = "chrom")

save(xset_pos, xset2_pos, xset3_pos, xset4_pos, xset5_pos, file = "../../xcms_strawberry_pos.RData")

library(CAMERA)
an_pos <- xsAnnotate(xset5_pos)
anF_pos <- groupFWHM(an_pos, perfwhm = 0.6)
anI_pos <- findIsotopes(anF_pos, mzabs = 0.01)
anIC_pos <- groupCorr(anI_pos, cor_eic_th = 0.75, graphMethod = "lpc")
anFA_pos <- findAdducts(anIC_pos, polarity = "positive")
peaklist_pos <- getPeaklist(anFA_pos)

save(an_pos, anF_pos, anI_pos, anIC_pos, anFA_pos, peaklist_pos, file = "../../CAMERA_strawberry_pos.RData")


setwd("~/winhome/Documents/05_collaborations_small_projects/Jose_fragaria/")
load("CAMERA_strawberry_pos.RData"); load("CAMERA_strawberry_neg.RData")

## keep only those features that are present in all biological replicates in 
## at least one of the classes (3 is max)
inds_neg <- which(colnames(peaklist_neg) == "X156_G"):which(colnames(peaklist_neg) == "SantaClara_W")
inds_neg_max <- apply(peaklist_neg[, inds_neg], 2, max)
peaklist_neg <- peaklist_neg[apply(peaklist_neg[, inds_neg], 1, function(x) any(x >= inds_neg_max)), ]

inds_pos <- which(colnames(peaklist_pos) == "X156_G"):which(colnames(peaklist_pos) == "SantaClara_W")
inds_pos_max <- apply(peaklist_pos[, inds_pos], 2, max)
peaklist_pos <- peaklist_pos[apply(peaklist_pos[, inds_pos], 1, function(x) any(x >= inds_pos_max)), ]

##peaklist_pos <- peaklist_pos[apply(peaklist_pos[, inds_pos], 1, max) >= 3, ]
##inds_neg <- which(colnames(peaklist_neg) == "X156_G"):which(colnames(peaklist_neg) == "SantaClara_W")
##peaklist_neg <- peaklist_neg[apply(peaklist_neg[, inds_neg], 1, max) >= 3, ]

## filter: use only features between 1.2-15.0 (rt is in seconds)
peaklist_pos <- peaklist_pos[peaklist_pos[, "rt"] > 1.2*60 & peaklist_pos[, "rt"] < 15.0*60, ]
peaklist_neg <- peaklist_neg[peaklist_neg[, "rt"] > 1.2*60 & peaklist_neg[, "rt"] < 15.0*60, ]

## check re-occurence of features (could be noise)
feat <- round(peaklist_neg[, "mz"], 3)
table(feat)[table(feat) > 10]
feat <- round(peaklist_pos[, "mz"], 3)
table(feat)[table(feat) > 10]

## exclude features that are mz=[101.23,112.984,114.933,116.927,161.045,117.928,118.926,137.008,144.923,158.990,
## 159.029,161.045,177.040,180.973,184.915,192.928,194.925,235.927,253.056,303.914,352.853,420.842] 
## (seems to be contamination, checked in chromatograms)
## mz=[116.927,192.928,235.927,253.056,352.853]
peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 101.022 & peaklist_neg[, "mz"] < 101.024),]
#peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 112.983 & peaklist_neg[, "mz"] < 112.985),]
#peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 114.932 & peaklist_neg[, "mz"] < 114.934),]
peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 115.001 & peaklist_neg[, "mz"] < 115.003),]
peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 116.926 & peaklist_neg[, "mz"] < 116.928),]
#peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 117.927 & peaklist_neg[, "mz"] < 117.929),]
#peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 118.925 & peaklist_neg[, "mz"] < 118.927),]
#peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 137.007 & peaklist_neg[, "mz"] < 137.009),]
#peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 144.922 & peaklist_neg[, "mz"] < 144.924),]
#peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 158.989 & peaklist_neg[, "mz"] < 158.991),]
#peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 159.028 & peaklist_neg[, "mz"] < 159.030),]
#peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 161.044 & peaklist_neg[, "mz"] < 161.046),]
#peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 177.039 & peaklist_neg[, "mz"] < 177.041),]
#peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 180.972 & peaklist_neg[, "mz"] < 180.974),]
#peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 184.914 & peaklist_neg[, "mz"] < 184.916),]
peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 192.927 & peaklist_neg[, "mz"] < 192.929),]
#peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 194.924 & peaklist_neg[, "mz"] < 194.926),]
peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 235.926 & peaklist_neg[, "mz"] < 235.928),]
peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 253.055 & peaklist_neg[, "mz"] < 253.057),]
#peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 303.913 & peaklist_neg[, "mz"] < 303.915),]
peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 352.852 & peaklist_neg[, "mz"] < 352.854),]
#peaklist_neg <- peaklist_neg[!(peaklist_neg[, "mz"] > 420.841 & peaklist_neg[, "mz"] < 420.843),]

## exclude features that are mz=[102.056,109.03,110.02,115.076,118.94,125.097,126.056,127.04,130.051,130.534,131.535,
##132.532,136.022,136.953,141.96,151.112,153.127,158.026-158.027,159.97,164.95,182.99,186.96,199.030,199.988,214.090,
## 215.093, 216.085,227.983
## mz=[130.051,186.96,214.90,215.093,216.085]
peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 102.055 & peaklist_pos[, "mz"] < 102.057), ] ##102.056
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 109.02 & peaklist_pos[, "mz"] < 109.04), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 110.01 & peaklist_pos[, "mz"] < 110.03), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 115.074 & peaklist_pos[, "mz"] < 115.078), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 118.93 & peaklist_pos[, "mz"] < 118.95), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 125.095 & peaklist_pos[, "mz"] < 125.099), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 126.054 & peaklist_pos[, "mz"] < 126.058), ]
peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 127.038 & peaklist_pos[, "mz"] < 127.042), ]
peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 130.050 & peaklist_pos[, "mz"] < 130.052), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 130.532 & peaklist_pos[, "mz"] < 130.536), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 131.533 & peaklist_pos[, "mz"] < 131.537), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 132.530 & peaklist_pos[, "mz"] < 132.534), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 136.020 & peaklist_pos[, "mz"] < 136.024), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 136.951 & peaklist_pos[, "mz"] < 136.955), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 141.95 & peaklist_pos[, "mz"] < 141.97), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 151.110 & peaklist_pos[, "mz"] < 151.114), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 153.125 & peaklist_pos[, "mz"] < 153.129), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 158.024 & peaklist_pos[, "mz"] < 158.029), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 159.96 & peaklist_pos[, "mz"] < 159.98), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 164.94 & peaklist_pos[, "mz"] < 164.96), ]
peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 186.95 & peaklist_pos[, "mz"] < 186.96), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 199.029 & peaklist_pos[, "mz"] < 199.031), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 199.987 & peaklist_pos[, "mz"] < 199.989), ]
peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 214.089 & peaklist_pos[, "mz"] < 214.091), ]
peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 215.092 & peaklist_pos[, "mz"] < 215.094), ]
peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 216.084 & peaklist_pos[, "mz"] < 216.086), ]
#peaklist_pos <- peaklist_pos[!(peaklist_pos[, "mz"] > 227.981 & peaklist_pos[, "mz"] < 227.985), ]


## divide each column by weight 
setwd("~/winhome/Documents/05_collaborations_small_projects/Jose_fragaria/")
weights <- read.table("strawberry_weight_list.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## for positive
inds_pos <- which(colnames(peaklist_pos) == "X156_G1_pos"):which(colnames(peaklist_pos) == "Santa_Clara_W3_pos")
col_pos <- unlist(strsplit(colnames(peaklist_pos[, inds_pos]), split = "_pos"))
weights_pos <- weights[match(col_pos, weights[, "sample"]), "weight_mg"] ## get weights vector for positive samples
peaklist_pos[, inds_pos] <- sweep(peaklist_pos[, inds_pos], 2, weights_pos, FUN = "/")

## for negative
inds_neg <- which(colnames(peaklist_neg) == "X156_G1_neg"):which(colnames(peaklist_neg) == "Santa_Clara_W3_neg")
col_neg <- unlist(strsplit(colnames(peaklist_neg[, inds_neg]), split = "_neg"))
weights_neg <- weights[match(col_neg, weights[, "sample"]), "weight_mg"] ## get weights vector for negative samples
peaklist_neg[, inds_neg] <- sweep(peaklist_neg[, inds_neg], 2, weights_neg, FUN = "/")

## divide each column by 75% quantile (to adjust for differences in sensitivity of each run)
peaklist_pos[, inds_pos] <- apply(peaklist_pos[, inds_pos], 2, function(x) x / quantile(x, 0.75))
peaklist_neg[, inds_neg] <- apply(peaklist_neg[, inds_neg], 2, function(x) x / quantile(x, 0.75))

## log2 transform intensity values
peaklist_pos[, inds_pos] <- log2(peaklist_pos[, inds_pos] + 1)
peaklist_neg[, inds_neg] <- log2(peaklist_neg[, inds_neg] + 1)

rownames(peaklist_pos) <- paste(round(peaklist_pos[, "mz"], 3), "/", round(peaklist_pos[, "rt"], 2), "_", peaklist_pos[, "pcgroup"], sep = "")
rownames(peaklist_neg) <- paste(round(peaklist_neg[, "mz"], 3), "/", round(peaklist_neg[, "rt"], 2), "_", peaklist_neg[, "pcgroup"], sep = "")

## filter for sd 
table(apply(peaklist_pos[, inds_pos], 1, sd) > 0.1) ## removes 279, keeps 15827 features
table(apply(peaklist_neg[, inds_neg], 1, sd) > 0.1) ## removes 655, keeps 17074 features
peaklist_pos <- peaklist_pos[apply(peaklist_pos[, inds_pos], 1, sd) > 0.1, ]
peaklist_neg <- peaklist_neg[apply(peaklist_neg[, inds_neg], 1, sd) > 0.1, ]

## filter low intensity features
table(apply(peaklist_pos[, inds_pos], 1, max) > 1) ## removes 2856, keeps 12971 features
table(apply(peaklist_neg[, inds_neg], 1, max) > 1) ## removes 4039, keeps 13035 features
peaklist_pos <- peaklist_pos[apply(peaklist_pos[, inds_pos], 1, max) > 1, ]
peaklist_neg <- peaklist_neg[apply(peaklist_neg[, inds_neg], 1, max) > 1, ]


pca <- prcomp(peaklist_pos[, inds_pos])
pdf("pca_fragaria_pos.pdf")
plot(pca$rotation[, 1], pca$rotation[, 2])
text(x = pca$rotation[, 1], y = pca$rotation[, 2], labels = rownames(pca$rotation))
dev.off()

pca <- prcomp(peaklist_neg[, inds_neg])
pdf("pca_fragaria_neg.pdf")
plot(pca$rotation[, 1], pca$rotation[, 2])
text(x = pca$rotation[, 1], y = pca$rotation[, 2], labels = rownames(pca$rotation))
dev.off()


write.table(peaklist_pos, file = "strawberry_peaklist_pos.csv", quote = FALSE, sep = "\t")
write.table(peaklist_neg, file = "strawberry_peaklist_neg.csv", quote = FALSE, sep = "\t")


## switch to golem
setwd("/home/naake/05_fragaria/")
## load function from github.com/tnaake/MetNet

## load peaktables
peaklist_pos <- read.table("strawberry_peaklist_pos.csv", header = TRUE, sep = "\t")
peaklist_neg <- read.table("strawberry_peaklist_neg.csv", header = TRUE, sep = "\t")

peaklist_neg[peaklist_neg[, "mz"] < 449.115 & peaklist_neg[, "mz"] > 449.10, c(1,4, 164:165)] ## pelargonidin-3-O gluycose
peaklist_neg[peaklist_neg[, "mz"] < 431.11 & peaklist_neg[, "mz"] > 431.09, c(1,4,164:165)]

peaklist_neg[peaklist_neg[, "mz"] < 465.115 & peaklist_neg[, "mz"] > 465.10, c(1,4, 164:165)] ## cyanidin 3 O glucoside
peaklist_neg[peaklist_neg[, "mz"] < 447.11 & peaklist_neg[, "mz"] > 447.08, c(1,4,164:165)]

peaklist_neg[peaklist_neg[, "mz"] < 509.14 & peaklist_neg[, "mz"] > 509.12, c(1,4, 164:165)] ## cyanidin 3 O glucoside
peaklist_neg[peaklist_neg[, "mz"] < 447.11 & peaklist_neg[, "mz"] > 447.08, c(1,4,164:165)]



## create data.frame with transformations
transformations <- rbind(
    c("transf_Hydrogenation/dehydrogenation", "H2", 2.0156500642, "?"), ## --> +
    c("transf_Acetylation (-H)", "C2H3O2", 59.0133043405, "?"), ## log P=-0.17 (acetic acid, Pubchem) --> +
    c("transf_Acetylation (-H2O)", "C2H2O",  42.0105646863, "?"), ## log P=-0.17 (acetic acid, Pubchem) --> +
    c("transf_benzoyl", "C6H4CO", 104.026213, "?"), ## --> +
    c("transf_hydroxy benzoyl", "C7H4O2", 120.021128, "?"), ## --> +
    c("transf_hydroxy benzyl", "C7H6O", 106.041863, "?"), ## --> +
    c("transf_galloyl", "C7H4O4", 152.010958, "?"), ## --> +
    c("transf_methylation", "CH2", 14.015649, "?"), ## --> +
    c("transf_methoxylation", "CH2O", 30.010564, "?"),
    c("transf_CHO2", "CHO2", 44.9976542763, "?"),
    c("transf_Glyoxylate (-H2O)", "C2O2",  55.9898292442, "?"), ## log Kow=-0.07 (glyoxylic acid, Pubchem)
    c("transf_Biotinyl (-H)", "C10H15N2O3S", 243.0803380482, "?"), ## log Kow=0.39 (biotin, Pubchem) --> +
    c("transf_Biotinyl (-H2O)", "C10H14N2O2S", 226.0775983940, "?"), ## log Kow=0.39 (biotin, Pubchem) --> +
    c("transf_C2H2", "C2H2", 26.0156500642, "?"),
    c("transf_Sulfate (-H2O)", "SO3", 79.9568145563, "?"), 
    c("transf_Isoprene addition (-H)", "C5H7", 67.0547752247, "+"), ## log Kow=2.42 (isoprene, Pubchem)
    c("transf_Ketol group (-H2O)", "C2H2O", 42.0105646863, "?"),
    c("transf_Primary amine", "NH2", 16.0187240694, "?"),
    c("transf_Secondary amine", "NH", 15.0108990373, "?"),
    c("transf_Tertiary amine", "N", 14.0030740052, "?"),
    c("transf_Hydroxylation (-H)", "O", 15.9949146221, "-"),
    c("transf_Malonyl group (-H2O)", "C3H2O3", 86.0003939305, "?"), ## log Kow=-0.81 (malonic acid, Pubchem) --> -
    c("transf_Urea addition (-H)", "CH3N2O", 59.0245377288, "-"), ## log Kow=-2.11 (urea, Pubchem)
    c("transf_D-ribose (-H2O) (ribosylation)", "C5H8O4", 132.0422587452, "-"),
    c("transf_Rhamnose (-H2O)", "C6H10O4", 146.0579088094, "-"),
    c("transf_Disaccharide (-H2O) #1", "C12H20O10", 324.105649, "-"),
    c("transf_Disaccharide (-H2O) #2", "C12H20O11", 340.1005614851, "-"),  
    c("transf_Glucuronic acid (-H2O)", "C6H8O6", 176.0320879894, "-"), ## log P=-2.57 (glucoronic acid, Pubchem)
    c("transf_Monosaccharide (-H2O)", "C6H10O5", 162.0528234315, "-"), ## log Kow=-3.24 (glucose, Pubchem)
    c("transf_Deoxyhexose (-H20)", "C6H10O4", 146.0579088094, "-"),
    c("transf_Pentose (-H2O)", "C5H8O4", 132.042260, "-"),
    c("transf_Trisaccharide (-H2O)", "C18H30O15", 486.1584702945, "-"),
    c("transf_Glucose-O-phosphate (-H2O)", "C6H11O8P", 242.0191538399, "-"), 
    c("transf_coumaroyl (-H2O)", "C9H6O2", 146.0367794368, "+"), ## log P=1.79 (coumaric acid, Pubchem)
    c("transf_coumaroyl hexose", "C15H16O7", 308.089599, "?"), ## --> -
    c("transf_caffeoyl", "C9H6O3", 162.031693, "+"),
    c("transf_feruloyl (-H2O)", "C9H6O2OCH2", 176.0473441231, "+"), ## log P=1.51 (ferulic acid, Pubchem)
    c("transf_sinapoyl (-H2O)", "C9H6O2OCH2OCH2", 206.0579088094, "+"), 
    c("transf_protocatechuoyl", "C7H4O3", 136.016043, "?"), ## --> +
    c("transf_quinic acid (-H2O)", "C7H10O5", 174.052824, "?"), ## logP=-2.007 (quinic acid, chemspicer) --> -
    c("transf_shikimic acid (-H2O)", "C7H8O4", 156.042260, "?"), 
    c("transf_ellagic acid (-H2O)", "	C14H4O7", 283.995705, "?")) ## log Kow=-2.05 (ellagic acid, Pubchem) --> -
    ##c("putrescine to spermidine (+C3H7N)", "C3H7N", 57.0578492299, "?"))

## check within each pcgroup, mass differences with respect to the molecular ion
adducts <- rbind(
    c("adduct_formic acid adduct", "CH2O2",  46.0054792, "?"),
    c("adduct_chlorine adduct", "Cl", 34.968853, "?"),
    c("adduct_sodium formate adduct", "NaHCO2", 67.9874246, "?"),
    c("adduct_Na adduct (+Na-2H)", "Na-H", 21.98146, "?"),
    c("adduct_K adduct (+K-2H)", "K-H", 37.9558834, "?"),
    c("adduct_isotopic+1_H", "isotopic peak 2H", 1.00628, "?"),
    c("adduct_isotopic+2_15N-14N", "isotopic peak 15N-14N", 0.997035, "?"),
    c("adduct_isotopic+2_18O-16O", "isotopic peak 18O-16O", 2.004244, "?"),
    c("adduct_isotopic+2_34S-32S", "isotopic peak 34S-32S", 1.995796, "?"),
    c("adduct_isotopic+1", "isotopic peak 13C1", 1.0033554, "?"),
    c("adduct_isotopic+2", "isotopic peak 13C2", 2.0067108, "?"),
    c("adduct_isotopic+3", "isotopic peak 13C3", 3.0100662, "?"), 
    c("adduct_loss of deoxyhexose (e.g. Rhamnose)", "C6H10O4", -146.0579088094, "?"), 
    c("adduct_loss of hexose", "C5H10O5", -162.0528234315, "?"),
    c("adduct_loss of pentose", "C5H8O4", -132.042260, "?"),
    c("adduct_loss of malonyl group", "C3H2O3", -86.0003939305, "?"),
    c("adduct_decarboxylation (loss from malonyl)", "CO2", -43.989830, "?"),
    c("adduct_loss of water", "H2O", -18.010565, "?"), 
    c("adduct_loss of COCH2", "COCH2", -42.0105646863, "?")
)


transformations <- data.frame(group = c(adducts[, 1], transformations[, 1]),
    formula = c(adducts[, 2], transformations[, 2]), 
    mass = c(as.numeric(adducts[, 3]), as.numeric(transformations[, 3])),
    rt = c(adducts[, 4], transformations[, 4]))

## use function createStructuralAdjacency and remove false positives by function rtCorrection
## pos
struct_adj_pos <- structural(x = peaklist_pos, transformation = transformations,
    ppm = 10, directed = TRUE)
struct_adj_pos <- rtCorrection(structural = struct_adj_pos, x = peaklist_pos, 
    transformation = transformations)
## neg
struct_adj_neg <- structural(x = peaklist_neg, transformation = transformations,
    ppm = 10, directed = TRUE)
struct_adj_neg <- rtCorrection(structural = struct_adj_neg, x = peaklist_neg,
    transformation = transformations)

## save
save(struct_adj_pos, file = "MetNet_strawberry_struct_adj_pos.RData")
save(struct_adj_neg, file = "MetNet_strawberry_struct_adj_neg.RData")

## use function statistical/threshold
inds_pos <- which(colnames(peaklist_pos) == "X156_G1_pos"):which(colnames(peaklist_pos) == "Santa_Clara_W3_pos")
inds_neg <- which(colnames(peaklist_neg) == "X156_G1_neg"):which(colnames(peaklist_neg) == "Santa_Clara_W3_neg")
peaklist_pos_cut <- peaklist_pos[, inds_pos]
peaklist_neg_cut <- peaklist_neg[, inds_neg]

models <- c("pearson", "pearson_partial", "spearman", "spearman_partial", 
    "clr", "aracne", "randomForest")

## apply the function statistical to create weighted adjacency matrices 
## per model
stat_adj_pos <- statistical(as.matrix(peaklist_pos_cut), model = models, 
    correlation_adjust = "BH") 
stat_adj_neg <- statistical(as.matrix(peaklist_neg_cut), model = models, 
    correlation_adjust = "BH")

save(stat_adj_pos, file = "MetNet_strawberry_stat_adj_pos.RData")
save(stat_adj_neg, file = "MetNet_strawberry_stat_adj_neg.RData")

## optional: set links within a pc_group to NA for pearson*/spearman*,
## clr/aracne/bayes/randomForest to reduce unmeaningful links
setLinkWithinPCgroupTo <- function(stat_adj, set_to) {
    
    if (is.list(stat_adj)) {
        pc <- strsplit(rownames(stat_adj[[1]]), split = "_")    
    } else {
        pc <- strsplit(rownames(stat_adj), split = "_")    
    }
    
    pc <- unlist(lapply(pc, "[", 2))
    pc_u <- unique(pc)    
    
    for (i in 1:length(pc_u)) {
        
        ## get inds of colnames (identical to rownames) that are equal to pc_pos_u
        inds <- pc == pc_u[i]
        
        ## iterate through each model and set to NaN
        if (is.list(stat_adj)) {
            for (j in 1:length(stat_adj)) {
                stat_adj[[j]][inds, inds] <- set_to
            }    
        } else {
            stat_adj[inds, inds] <- set_to
        }
        
    }
    
    return(stat_adj)
}

## apply function
stat_adj_pos_NA <- setLinkWithinPCgroupTo(stat_adj_pos, set_to = NaN)
stat_adj_neg_NA <- setLinkWithinPCgroupTo(stat_adj_neg, set_to = NaN)

## apply the function threshold to create unweighted adjacency matrices
## type = "threshold" 
## define thresholds
pdf("hist_pos_pearson.pdf")
hist(stat_adj_pos_NA[["pearson"]])
dev.off()

pdf("hist_pos_pearson_partial.pdf")
hist(stat_adj_pos_NA[["pearson_partial"]])
dev.off()

pdf("hist_pos_spearman.pdf")
hist(stat_adj_pos_NA[["spearman"]])
dev.off()

pdf("hist_pos_spearman_partial.pdf")
hist(stat_adj_pos_NA[["spearman_partial"]])
dev.off()

pdf("hist_pos_clr.pdf")
hist(stat_adj_pos_NA[["clr"]])
dev.off()

pdf("hist_pos_aracne.pdf")
hist(stat_adj_pos_NA[["aracne"]])
dev.off()

pdf("hist_pos_randomForest.pdf")
hist(stat_adj_pos_NA[["randomForest"]])
dev.off()

pdf("hist_neg_pearson.pdf")
hist(stat_adj_neg_NA[["pearson"]])
dev.off()

pdf("hist_neg_pearson_partial.pdf")
hist(stat_adj_neg_NA[["pearson_partial"]])
dev.off()

pdf("hist_neg_spearman.pdf")
hist(stat_adj_neg_NA[["spearman"]])
dev.off()

pdf("hist_neg_spearman_partial.pdf")
hist(stat_adj_neg_NA[["spearman_partial"]])
dev.off()

pdf("hist_neg_clr.pdf")
hist(stat_adj_neg_NA[["clr"]])
dev.off()

pdf("hist_neg_aracne.pdf")
hist(stat_adj_neg_NA[["aracne"]])
dev.off()

pdf("hist_neg_randomForest.pdf")
hist(stat_adj_neg_NA[["randomForest"]])
dev.off()

## check thresholds
table(stat_adj_neg_NA[["pearson"]][upper.tri(stat_adj_neg_NA[["pearson"]])] > 0.75)
table(stat_adj_pos_NA[["pearson"]][upper.tri(stat_adj_pos_NA[["pearson"]])] > 0.75)
table(stat_adj_neg_NA[["pearson_partial"]][upper.tri(stat_adj_neg_NA[["pearson_partial"]])] > 0.5)
table(stat_adj_pos_NA[["pearson_partial"]][upper.tri(stat_adj_pos_NA[["pearson_partial"]])] > 0.5)
table(stat_adj_neg_NA[["spearman"]][upper.tri(stat_adj_neg_NA[["spearman"]])] > 0.75)
table(stat_adj_pos_NA[["spearman"]][upper.tri(stat_adj_pos_NA[["spearman"]])] > 0.75)
table(stat_adj_neg_NA[["spearman_partial"]][upper.tri(stat_adj_neg_NA[["spearman_partial"]])] > 0.5)
table(stat_adj_pos_NA[["spearman_partial"]][upper.tri(stat_adj_pos_NA[["spearman_partial"]])] > 0.5)
table(stat_adj_neg_NA[["clr"]][upper.tri(stat_adj_neg_NA[["clr"]])] > 2)
table(stat_adj_pos_NA[["clr"]][upper.tri(stat_adj_pos_NA[["clr"]])] > 2)
table(stat_adj_neg_NA[["aracne"]][upper.tri(stat_adj_neg_NA[["aracne"]])] > 0.05)
table(stat_adj_pos_NA[["aracne"]][upper.tri(stat_adj_pos_NA[["aracne"]])] > 0.05)
table(stat_adj_neg_NA[["randomForest"]][upper.tri(stat_adj_neg_NA[["randomForest"]])] > 0.001)
table(stat_adj_pos_NA[["randomForest"]][upper.tri(stat_adj_pos_NA[["randomForest"]])] > 0.001)

args <- list("pearson" = 0.75, "pearson_partial" = 0.5, "spearman" = 0.75,
    "spearman_partial" = 0.5, "clr" = 2, "aracne" = 0.05, 
    "randomForest" = 0.001, threshold = 1)
stat_adj_pos_thr <- threshold(statistical = stat_adj_pos_NA, type = "threshold",
    args = args)
stat_adj_neg_thr <- threshold(statistical = stat_adj_neg_NA, type = "threshold",
    args = args)
 
## type = "top1" 
args_top <- list(n = 500000)
stat_adj_pos_top1 <- threshold(statistical = stat_adj_pos, type = "top1",
    args = args_top)
stat_adj_neg_top1 <- threshold(statistical = stat_adj_neg, type = "top1",
    args = args_top)
 
## type = "top2"
stat_adj_pos_top2 <- threshold(statistical = stat_adj_pos, type = "top2",
    args = args_top)
stat_adj_neg_top2 <- threshold(statistical = stat_adj_neg, type = "top2",
    args = args_top)

## type = "mean"
stat_adj_pos_mean <- threshold(statistical = stat_adj_pos, type = "mean",
    args = args_top)
stat_adj_neg_mean <- threshold(statistical = stat_adj_neg, type = "mean",
    args = args_top)


save(stat_adj_pos_thr, stat_adj_pos_top1, stat_adj_pos_top2, stat_adj_pos_mean, 
    file = "MetNet_strawberry_stat_adj_thr_pos.RData")
save(stat_adj_neg_thr, stat_adj_neg_top1, stat_adj_neg_top2, stat_adj_neg_mean, 
     file = "MetNet_strawberry_stat_adj_thr_neg.RData")


## set links between features belonging to different pcgroups in struct_adj to
## 0 for links that contain the "adduct_" pattern (type = "inter")
## set links between features belonging to the same pcgroup in struct_adj to 0
## for links that contain the "transf_" pattern (type = "intra")
setLinkPCgroupTo0 <- function(struct_adj, type = c("inter", "intra"), pattern = "adduct_") {
    pc <- strsplit(rownames(struct_adj[[1]]), split = "_")
    pc <- unlist(lapply(pc, "[", 2))
    pc_u <- unique(pc)    

    inds <- apply(struct_adj[[2]], 2, grepl, pattern = pattern)
    pc_inds <- inds <- which(inds, arr.ind = TRUE)

    pc_inds[, "row"] <- pc[inds[, "row"]]
    pc_inds[, "col"] <- pc[inds[, "col"]]
    
    if (type == "inter") {
        ## remove between different pcgroups (pattern = "adduct_")
        inds_remove <- which(pc_inds[, "row"] != pc_inds[, "col"])
        
        for (i in 1:length(inds_remove)) {
            ind_row <- inds[inds_remove[i], "row"]
            ind_col <- inds[inds_remove[i], "col"]
            struct_adj[[1]][ind_row, ind_col] <- 0
            struct_adj[[2]][ind_row, ind_col] <- ""
        }
    }
        
    if (type == "intra") {
        ## remove within the same pcgroup (pattern = "transf_")
        inds_remove <- which(pc_inds[, "row"] == pc_inds[, "col"])
        
        for (i in 1:length(inds_remove)) {
            ind_row <- inds[inds_remove[i], "row"]
            ind_col <- inds[inds_remove[i], "col"]
            struct_adj[[1]][ind_row, ind_col] <- 0
            struct_adj[[2]][ind_row, ind_col] <- ""
        }
    }
    
    return(struct_adj)
}

struct_adj_neg_mod <- setLinkPCgroupTo0(struct_adj_neg, type = "inter", 
    pattern = "adduct_")
struct_adj_pos_mod <- setLinkPCgroupTo0(struct_adj_pos, type = "inter",
    pattern = "adduct_")
##struct_adj_neg_mod <- setLinkPCgroupTo0(struct_adj_neg_mod, type = "intra", 
##    pattern = "transf_")
##struct_adj_pos_mod <- setLinkPCgroupTo0(struct_adj_pos_mod, type = "intra",
##    pattern = "transf_")

## set all intra pcgroup links in stat_adj to 1 in order to make sure that all 
## links from struct_adj for intra pcgroups links are taken 
stat_adj_pos_thr_mod <- setLinkWithinPCgroupTo(stat_adj_pos_thr, set_to = 1)
stat_adj_neg_thr_mod <- setLinkWithinPCgroupTo(stat_adj_neg_thr, set_to = 1)

## use function combine to combine the structural and statistical information
cons_adj_pos <- combine(structural = struct_adj_pos_mod, 
    statistical = stat_adj_pos_thr_mod)
cons_adj_neg <- combine(structural = struct_adj_neg_mod, 
    statistical = stat_adj_neg_thr_mod)

save(cons_adj_pos, file = "MetNet_strawberry_cons_adj_pos.RData")
save(cons_adj_neg, file = "MetNet_strawberry_cons_adj_neg.RData")

## only retain edges that link to at least two other mass features
g_pos <- graph_from_adjacency_matrix(cons_adj_pos[[1]], mode = "directed")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg[[1]], mode = "directed")
comp_g_pos <- components(g_pos)
comp_g_neg <- components(g_neg)
inds_cut_pos <- comp_g_pos$membership %in% which(comp_g_pos$csize > 2)
inds_cut_neg <- comp_g_neg$membership %in% which(comp_g_neg$csize > 2)


cons_adj_pos_cut <- cons_adj_pos
cons_adj_neg_cut <- cons_adj_neg
cons_adj_pos_cut[[1]] <- cons_adj_pos[[1]][inds_cut_pos, inds_cut_pos]
cons_adj_pos_cut[[2]] <- cons_adj_pos[[2]][inds_cut_pos, inds_cut_pos]
cons_adj_neg_cut[[1]] <- cons_adj_neg[[1]][inds_cut_neg, inds_cut_neg]
cons_adj_neg_cut[[2]] <- cons_adj_neg[[2]][inds_cut_neg, inds_cut_neg]

g_pos <- graph_from_adjacency_matrix(cons_adj_pos_cut[[1]], mode = "directed")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg_cut[[1]], mode = "directed")

plot(g_pos, vertex.label.cex = 0.1, vertex.size = 0.1, edge.arrow.width = 0.1,
    edge.arrow.size = 0.1)

## remove edges between M1+H and M2+Na, M3+H and M4+K, in general when there is 
## a link between an adduct and a molecular ion

## node with outgoing transf link should not have outgoing adduct link when
## the adduct-linking vertex has a ingoing adduct link

## node with ingoing transf link should not have ingoing adduct link when the 
## adduct-linking vertex has an outgoing adduct link

## pos mode: isotopic+1, isotopic+2, isotopic+3, Na adduct, K adduct, 
## formic acid adduct
## neg mode: isotopic+1, isotopic+2, isotopic+3, formic acid, Na adduct,
## Sodium formate 


##          M1+H    M1+Na   M2+H    M2+Na   M3+H    M3+Na
## M1+H     0       1       0       0       0       1       
## M1+Na    0       0       0       0       0       0
## M2+H     1       0       0       1       0       0
## M2+Na    0       1       0       0       1       0
## M3+H     0       1       0       0       0       1
## M3+Na    0       0       0       0       0       0  

# mat_test <- matrix(c(c(0, 1, 1, 0, 0, 1), ## M1+H
#                      c(1, 0, 1, 1, 1, 1), ## M1+Na
#                      c(1, 1, 0, 1, 0, 0), ## M2+H
#                      c(0, 1, 1, 0, 1, 0), ## M2+Na
#                      c(0, 1, 0, 1, 0, 1), ## M3+H
#                      c(0, 1, 0, 0, 1, 0)), ## M3+Na
#                    
mat_test <- matrix(c(c(0, 1, 0, 0, 0, 1), ## M1+H
                     c(0, 0, 0, 0, 0, 0), ## M1+Na
                     c(1, 0, 0, 1, 0, 0), ## M2+H
                     c(0, 1, 0, 0, 0, 0), ## M2+Na
                     c(0, 1, 0, 0, 0, 1), ## M3+H
                     c(0, 0, 0, 0, 0, 0)), ## M3+Na
                   byrow = TRUE, ncol = 6, nrow = 6)
mat_test_transformation <- matrix(c(c("", "adduct_Na", "", "", "", "transf_gluc_F"), ## M1+H
                     c("", "", "", "", "", ""), ## M1+Na
                     c("transf_gluc_T", "", "", "adduct_Na", "", ""), ## M2+H
                     c("", "transf_gluc_T", "", "", "", ""), ## M2+Na
                     c("", "transf_gluc_F", "", "", "", "adduct_Na"), ## M3+H
                     c("", "", "", "", "", "")), ## M3+Na
                   byrow = TRUE, ncol = 6, nrow = 6)
colnames(mat_test) <- c("M1+H", "M1+Na", "M2+H", "M2+Na", "M3+H", "M3+Na")
rownames(mat_test) <- colnames(mat_test)


## node with outgoing transf link should not have outgoing adduct link when
## the adduct-linking vertex has an ingoing adduct link of same type

## do for all components
## iterate through rows

removeFalseLinksAdducts <- function(mat_l) {
    
    mat_num <- mat_l[[1]]
    mat_char <- mat_l[[2]]
    
    for (i in seq_len(nrow(mat_num))) {
        
        rows_i <- mat_char[i, ]
        adduct_i <- rows_i[grep(rows_i, pattern = "adduct_")]
        
        ## if any vertex has outgoing adduct link, check all cols where rows_i
        ## has transf if ingoing adduct link exists for adduct-linking feature
        if (length(adduct_i) > 0) { 
            
            inds <- grep(rows_i, pattern = "transf_")
            
            for (j in inds) {
                cols_j <- mat_char[, j]
                adduct_j <- cols_j[grep(cols_j, pattern = "adduct")]
                print(i)
                print(j)
                if (length(adduct_j) > 0)
                    if (adduct_j %in% adduct_i)
                        mat_num[i, j] <- 0 ## remove in mat_test[i, ]
            }
        }
    }
    
    ## set mat_char to "" where mat_num == 0
    mat_char[which(mat_num == 0)] <- ""
    
    return(list(mat_num, mat_char))
}

cons_adj_pos_cut_rem <- removeFalseLinksAdducts(cons_adj_pos_cut)


g_pos <- graph_from_adjacency_matrix(cons_adj_pos_cut_rem[[1]], mode = "directed")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg_cut[[1]], mode = "directed")

plot(g_pos, vertex.label.cex = 0.1, vertex.size = 0.1, edge.arrow.width = 0.1, edge.arrow.size = 0.1)


tmp <- removeFalseLinksAdducts(list(mat_test, mat_test_transformation))[[1]]
net_test <- graph_from_adjacency_matrix(tmp, mode = "directed")
plot(net_test)




## reduce edges using pc_groups from peaklist_pos and peaklist_neg
## combine edges from one group to one edge representing the pc group
## create a link between a pc_group when there is at least one link between members of the pc_group
reduce_cons <- function (cons_mat) {
    
    ## assign the two entries to cons_mat_bin and cons_mat_type
    cons_mat_bin <- cons_mat[[1]]
    cons_mat_type <- cons_mat[[2]]
    
    ## get pc_groups
    pc_cons <- unlist(lapply(strsplit(colnames(cons_mat_bin), split = "_"), "[", 2))
    pc_cons_u <- unique(pc_cons)
    
    ## create adjacency matrices
    mat_pc_bin <- matrix(0, nrow = length(pc_cons_u), ncol = length(pc_cons_u))
    mat_pc_type <- matrix("", nrow = length(pc_cons_u), ncol = length(pc_cons_u))
    colnames(mat_pc_bin) <- rownames(mat_pc_bin) <- pc_cons_u
    colnames(mat_pc_type) <- rownames(mat_pc_type) <- pc_cons_u
    
    for (i in 1:length(pc_cons_u)) {
        inds <- which(pc_cons == pc_cons_u[i])
        ## get pc_groups to which members of pc_cons_u[i] link to
        pc_col <- NULL
        pc_row <- NULL
        pc_col <- rownames(which(cons_mat[[1]][, inds] == 1, arr.ind = TRUE))
        
        if (!is.null(pc_col)) {
            pc_col <- unlist(lapply(strsplit(pc_col, split = "_"), "[", 2))
        }
        
        pc_row <- rownames(which(t(cons_mat[[1]][inds, ]) == 1, arr.ind = TRUE))
        
        if (!is.null(pc_row)) {
            pc_row <- unlist(lapply(strsplit(pc_row, split = "_"), "[", 2))
        }
        
        if (!is.null(pc_col)) {
            mat_pc_bin[pc_col, pc_cons_u[i]] <- 1
            
            ## collapse-paste all inter feature links
            for (j in 1:length(pc_col)) {
                type_paste <- cons_mat[[2]][pc_cons %in% pc_col[j], pc_cons == pc_cons_u[[i]]]
                type_paste <- names(table(as.vector(type_paste)))
                type_paste <- type_paste[type_paste != ""]
                type_paste <- paste(type_paste, collapse = "/")
                mat_pc_type[pc_col[j], pc_cons_u[i]] <- type_paste
            }
        }
        
        if (!is.null(pc_row)) {
            mat_pc_bin[pc_cons_u[i], pc_row] <- 1
            
            ## collapse-paste all inter feature links
            for (j in 1:length(pc_row)) {
                type_paste <- cons_mat[[2]][pc_cons == pc_cons_u[[i]], pc_cons %in% pc_row[j]]
                type_paste <- names(table(as.vector(type_paste)))
                type_paste <- type_paste[type_paste != ""]
                type_paste <- paste(type_paste, collapse = "/")
                mat_pc_type[pc_cons_u[i], pc_row[j]] <- type_paste
            }
        }
    }
    return(list(mat_pc_bin, mat_pc_type))
}


## reduce the consensus adjacency matrices and save to file
cons_neg_red <- reduce_cons(cons_adj_neg)
cons_pos_red <- reduce_cons(cons_adj_pos)
diag(cons_neg_red[[1]]) <- 0
diag(cons_pos_red[[1]]) <- 0
save(cons_neg_red, file="MetNet_strawberry_cons_neg_red.RData")
save(cons_pos_red, file="MetNet_strawberry_cons_pos_red.RData")

remove_non_connected <- function(cons_red) {
    ## find nodes that do not connect to others and remove
    inds_1 <- apply(cons_red[[1]], 1, function(x) all(x == 0))
    inds_2 <- apply(cons_red[[1]], 2, function(x) all(x == 0))
    inds <- intersect(which(inds_1), which(inds_2))
    cons_red_1 <- cons_red[[1]][-inds, -inds]
    cons_red_2 <- cons_red[[2]][-inds, -inds]
    return(list(cons_red_1, cons_red_2))
}

cons_neg_red_0 <- remove_non_connected(cons_neg_red)
cons_pos_red_0 <- remove_non_connected(cons_pos_red)

g_neg <- graph_from_adjacency_matrix(cons_neg_red_0[[1]], mode="undirected")
g_pos <- graph_from_adjacency_matrix(cons_pos_red_0[[1]], mode="undirected")

pdf("network_cons_neg_red_0.pdf")
plot(g_neg, vertex.size=0.5, vertex.label.cex=0.1)
dev.off()

## 9, 10 only in condition 1
## 7 only in condition 2
## 8 in condition 2 and 3
## 6 only in condition 3
## 1-5 in condition 1-3
peaklist_neg[peaklist_neg[, "pcgroup"] == 2, ]

inds_samples_neg <- which(colnames(peaklist_neg) == "X156_G"):which(colnames(peaklist_neg) == "SantaClara_W")
inds_samples_pos <- which(colnames(peaklist_pos) == "X156_G"):which(colnames(peaklist_pos) == "SantaClara_W")


sample_presence_neg <- rbind(c("X156_G", 3), c("X156_R", 3), c("X156_W", 3), c("X191_G", 3), 
    c("X191_R", 3), c("X191_W", 3), c("X196_G", 3), c("X196_R", 3), c("X196_W", 3), c("X282_G", 3), 
    c("X282_R", 3), c("X282_W", 3), c("X591_G", 3), c("X591_R", 3), c("X591_W", 3), c("X595_G", 3), 
    c("X595_R", 3), c("X595_W", 3), c("X660_G", 3), c("X660_R", 3), c("X660_W", 3), c("Amiga_G", 3),
    c("Amiga_R", 3), c("Amiga_W", 3), c("Benicia_G", 3), c("Benicia_R", 3), c("Benicia_W", 3),
    c("Candonga_G", 2), c("Candonga_R", 3), c("Candonga_W", 3), c("Fontanilla_G", 3), c("Fontanilla_R", 3),
    c("Fontanilla_W", 3), c("Fuentepina_G", 3), c("Fuentepina_R", 3), c("Fuentepina_W", 3), 
    c("SantaClara_G", 3), c("SantaClara_R", 3), c("SantaClara_W", 3))
sample_presence_pos <- rbind(c("X156_G", 3), c("X156_R", 3), c("X156_W", 3), c("X191_G", 3), 
    c("X191_R", 3), c("X191_W", 3), c("X196_G", 3), c("X196_R", 3), c("X196_W", 3), c("X282_G", 3), 
    c("X282_R", 3), c("X282_W", 3), c("X591_G", 3), c("X591_R", 3), c("X591_W", 3), c("X595_G", 3), 
    c("X595_R", 3), c("X595_W", 3), c("X660_G", 3), c("X660_R", 3), c("X660_W", 3), c("Amiga_G", 3),
    c("Amiga_R", 3), c("Amiga_W", 3), c("Benicia_G", 3), c("Benicia_R", 3), c("Benicia_W", 3),
    c("Camarosa_G", 3), c("Camarosa_R", 3), c("Camarosa_W", 3), c("Candonga_G", 3), 
    c("Candonga_R", 3), c("Candonga_W", 3), c("Fontanilla_G", 3), c("Fontanilla_R", 3),
    c("Fontanilla_W", 3), c("Fuentepina_G", 3), c("Fuentepina_R", 3), c("Fuentepina_W", 3), 
    c("SantaClara_G", 3), c("SantaClara_R", 3), c("SantaClara_W", 3))

create_presence_mat <- function(cons_red_0, peaklist, sample_presence) {
    
    sample_presence <- data.frame(geno=as.character(sample_presence[,1]), num=as.numeric(sample_presence[,2]))
    presence_mat <- matrix(0, ncol=length(sample_presence$geno), nrow=0)
    colnames(presence_mat) <- sample_presence$geno
    
    for (i in 1:nrow(cons_red_0)) {
        print(i)
        pcgroup <- rownames(cons_red_0)[i]     
        peak_pc <- peaklist[peaklist[, "pcgroup"] == pcgroup, ]
        ## define a metabolite as present when it is found in more than two thirds of members of a pcgroup
        peak_pc_pres <- peak_pc[, as.character(sample_presence$geno)] >= sample_presence$num
        peak_pc_pres <- apply(peak_pc_pres, 2, function(x) sum(x) / dim(peak_pc_pres)[1] >= 2/3)
        peak_pc_mat <- matrix(peak_pc_pres, nrow=1)
        rownames(peak_pc_mat) <- pcgroup
        presence_mat <- rbind(presence_mat, ifelse(peak_pc_pres, 1, 0))
        rownames(presence_mat)[i] <- pcgroup
    }
    return(presence_mat)
}

## create function that iterates through all pcgroups, it returns the
## row which has highest mean intensity value of this pcgroup
create_reduced_peaklist <- function(cons_red_0, peaklist, inds) {
    peaklist_red <- matrix(0, ncol=ncol(peaklist), nrow=0)
    colnames(peaklist_red) <- colnames(peaklist)
    pcgroup <- rownames(cons_red_0)
    for (i in 1:length(pcgroup)) {
        inds_row <- which(peaklist[, "pcgroup"] == pcgroup[i])
        
        if (length(inds_row) > 1) {
            mean_row <- apply(peaklist[inds_row, inds], 1, mean)
            peaklist_add <- peaklist[inds_row,]##[which.max(mean_row)],]    
            peaklist_red <- rbind(peaklist_red, peaklist_add)
        } else { ## remove singletons
            peaklist_add <- peaklist[inds_row, ]
        }
        
        
    }
    return(peaklist_red)
}
peaklist_neg_red <- create_reduced_peaklist(cons_neg_red_0, peaklist_neg, inds_neg)

source("~/winhome/Documents/R_functions/fuzzy_clustering.R")

estimClustNum_neg <- estimClustNum(peaklist_neg_red[, inds_neg], maxClust=30)
ClustComp_neg <- ClustComp(peaklist_neg_red[, inds_neg], NSs=50, NClust=estimClustNum_neg$numclust, cores=1)
plot_clusters(ClustComp_neg, peaklist_neg_red[, inds_neg], filename="clustering_clustComp_neg.pdf")

peaklist_neg_red_z <- t(apply(peaklist_neg_red[, inds_neg], 1, function(x) (x - mean(x)) / sd(x)))

pheatmap(peaklist_neg_red_z[names(which(ClustComp_neg$Bestcl$cluster == 5)), ], cluster_cols=FALSE)


library(amap)
peaklist_neg_red_z <- t(apply(peaklist_neg_red[, inds_neg], 1, function(x) (x - mean(x)) / sd(x)))
clust_neg <- amap::hcluster(peaklist_neg_red_z, method="correlation")
plot(clust_neg)

pheatmap(peaklist_neg_red_z, clustering_distance_rows="correlation", clustering_distance_cols="correlation")

M.pca <- prcomp(peaklist_neg_red_z, scale = TRUE, center = TRUE)
loadings <- data.frame(M.pca$rotation[, 1:2], names=rownames(M.pca$rotation))
imp <- peaklist_neg_red[names(sort(sqrt((M.pca$x[,1]^2) + (M.pca$x[,2])^2), decreasing=T)[1:100]), "pcgroup"]
cons_neg_red_0_imp <- cons_neg_red_0[imp, imp]
ind <- apply(cons_neg_red_0_imp, 1, function(x) any(x > 0))
cons_neg_red_0_imp <- cons_neg_red_0_imp[ind, ind]

g <- igraph::graph_from_adjacency_matrix(cons_neg_red_0_imp, mode="undirected")
plot(g, vertex.label.cex=0.5, vertex.size=1)


ggplot(loadings, aes(x=PC1, y=PC2)) + 
    ##geom_text(aes(x=PC1, y=PC2+0.008,label=names), size = 2.5) + 
    geom_point(size=2.0) +
    scale_shape_manual("", values=shapes) + 
    scale_fill_manual("", values=fillings) + ## was values=fillings
    xlab("PC1 (82.77%)") + ylab("PC2 (4.23%)")


presence_mat_neg <- create_presence_mat(cons_neg_red_0, peaklist_neg, sample_presence_neg)
presence_mat_pos <- create_presence_mat(cons_pos_red_0, peaklist_pos, sample_presence_pos)

extract_module <- function(reduced_adj, presence_mat, cond) {
    col_c <- which(colnames(presence_mat) %in% cond)
    col_c_n <- which(!colnames(presence_mat) %in% cond)
    ## which cells per row are not present in cond_not
    ## get rows if the condition is true
    feat <- apply(presence_mat, 1, function(x) all(x[col_c_n] != 1))
    feat <- which(feat)
    
    reduced_adj_r <- matrix(reduced_adj[feat, feat], ncol=length(feat), byrow=TRUE)
    rownames(reduced_adj_r) <- colnames(reduced_adj_r) <- names(feat)
    
    return(reduced_adj_r)
}

wild_species <- colnames(peaklist_neg[, inds_samples_neg])[grep(colnames(peaklist_neg[, inds_samples_neg]), pattern="X[0-9]")]

plot(graph_from_adjacency_matrix(extract_module(reduced_mat, presence_mat, cond=c("cond1", "cond2", "cond3")), mode="undirected"))

plot(graph_from_adjacency_matrix(extract_module(cons_neg_red_0, presence_mat_neg, cond=wild_species), mode="undirected"))
##
library(igraph)
g_pos <- graph_from_adjacency_matrix(cons_adj_pos, mode="undirected")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg, mode="undirected")
comp_g_pos <- components(g_pos)
comp_g_neg <- components(g_neg)
plot(g_pos, edge.width=5, vertex.label.cex=0.5, edge.color="grey")
