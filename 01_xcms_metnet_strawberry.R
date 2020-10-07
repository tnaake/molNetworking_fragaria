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

## create peaklist for W, G, R
cols_W <- colnames(peaklist_pos)[grep(colnames(peaklist_pos), pattern = "_W1_|_W2_|_W3_")]
peaklist_pos_W <- peaklist_pos[, c("mz", "rt", cols_W)]
cols_W <- colnames(peaklist_neg)[grep(colnames(peaklist_neg), pattern = "_W1_|_W2_|_W3_")]
peaklist_neg_W <- peaklist_neg[, c("mz", "rt", cols_W)]
cols_G <- colnames(peaklist_pos)[grep(colnames(peaklist_pos), pattern = "_G1_|_G2_|_G3_")]
peaklist_pos_G <- peaklist_pos[, c("mz", "rt", cols_G)]
cols_G <- colnames(peaklist_neg)[grep(colnames(peaklist_neg), pattern = "_G1_|_G2_|_G3_")]
peaklist_neg_G <- peaklist_neg[, c("mz", "rt", cols_G)]
cols_R <- colnames(peaklist_pos)[grep(colnames(peaklist_pos), pattern = "_R1_|_R2_|_R3_")]
peaklist_pos_R <- peaklist_pos[, c("mz", "rt", cols_R)]
cols_R <- colnames(peaklist_neg)[grep(colnames(peaklist_neg), pattern = "_R1_|_R2_|_R3_")]
peaklist_neg_R <- peaklist_neg[, c("mz", "rt", cols_R)]

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

## check within each pcgroup, mass differences with respect to the M+H
adducts_pos <- rbind(
    c("adduct_formic acid adduct", "CH2O2",  46.0054792, "?"), ## only positive
    c("adduct_ammonium", "NH4-H", 18.034374-1.007825, "?"),
    c("adduct_acetonitril H+CH3CN", "CH3CN", 41.026549, "?"),
    c("adduct_acetonitril H+CH3CN+1H2O", "CH3CNH2O", 59.037113, "?"),
    c("adduct_sodium formate adduct", "NaHCO2", 67.9874246, "?"),
    c("adduct_Na adduct (+Na+)", "Na-H", 22.989770-1.007825, "?"), ## only positive
    c("adduct_K adduct (+K+)", "K-H", 38.963708-1.007825, "?"), ## only positive
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

## check within each pcgroup, mass differences with respect to the M-H
adducts_neg <- rbind(
    c("adduct_formate HCOO-", "HCOO+H", 44.997654+1.007825, "?"),
    c("adduct_chlorine adduct", "Cl+H", 34.968853+1.007825, "?"),
    c("adduct_M-2H+Na", "Na-H", 22.989770-1.00782, "?"),
    c("adduct_sodium formate adduct", "NaHCO2", 67.9874246, "?"),
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

## positive mode
transformations_pos <- data.frame(
    group = c(adducts_pos[, 1], transformations[, 1]),
    formula = c(adducts_pos[, 2], transformations[, 2]), 
    mass = c(as.numeric(adducts_pos[, 3]), as.numeric(transformations[, 3])),
    rt = c(adducts_pos[, 4], transformations[, 4]))

## negative mode 
transformations_neg <- data.frame(
    group = c(adducts_neg[, 1], transformations[, 1]),
    formula = c(adducts_neg[, 2], transformations[, 2]), 
    mass = c(as.numeric(adducts_neg[, 3]), as.numeric(transformations[, 3])),
    rt = c(adducts_neg[, 4], transformations[, 4]))

## use function createStructuralAdjacency and remove false positives by 
## function rtCorrection
## pos
struct_adj_pos <- structural(x = peaklist_pos, 
    transformation = transformations_pos, ppm = 10, directed = TRUE)
struct_adj_pos <- rtCorrection(structural = struct_adj_pos, x = peaklist_pos, 
    transformation = transformations_pos)
## neg
struct_adj_neg <- structural(x = peaklist_neg, 
    transformation = transformations_neg, ppm = 10, directed = TRUE)
struct_adj_neg <- rtCorrection(structural = struct_adj_neg, x = peaklist_neg,
    transformation = transformations_neg)

## save
save(struct_adj_pos, file = "MetNet_strawberry_struct_adj_pos.RData")
save(struct_adj_neg, file = "MetNet_strawberry_struct_adj_neg.RData")

## use function statistical/threshold
inds_pos <- which(colnames(peaklist_pos) == "X156_G1_pos"):which(colnames(peaklist_pos) == "Santa_Clara_W3_pos")
inds_neg <- which(colnames(peaklist_neg) == "X156_G1_neg"):which(colnames(peaklist_neg) == "Santa_Clara_W3_neg")
peaklist_pos_cut <- peaklist_pos[, inds_pos]
peaklist_neg_cut <- peaklist_neg[, inds_neg]
inds_pos <- which(colnames(peaklist_pos_G) == "X156_G1_pos"):which(colnames(peaklist_pos_G) == "Santa_Clara_G3_pos")
inds_neg <- which(colnames(peaklist_neg_G) == "X156_G1_neg"):which(colnames(peaklist_neg_G) == "Santa_Clara_G3_neg")
peaklist_pos_G_cut <- peaklist_pos_G[, inds_pos]
peaklist_neg_G_cut <- peaklist_neg_G[, inds_neg]
inds_pos <- which(colnames(peaklist_pos_W) == "X156_W1_pos"):which(colnames(peaklist_pos_W) == "Santa_Clara_W3_pos")
inds_neg <- which(colnames(peaklist_neg_W) == "X156_W1_neg"):which(colnames(peaklist_neg_W) == "Santa_Clara_W3_neg")
peaklist_pos_W_cut <- peaklist_pos_W[, inds_pos]
peaklist_neg_W_cut <- peaklist_neg_W[, inds_neg]
inds_pos <- which(colnames(peaklist_pos_R) == "X156_R1_pos"):which(colnames(peaklist_pos_R) == "Santa_Clara_R3_pos")
inds_neg <- which(colnames(peaklist_neg_R) == "X156_R1_neg"):which(colnames(peaklist_neg_R) == "Santa_Clara_R3_neg")
peaklist_pos_R_cut <- peaklist_pos_R[, inds_pos]
peaklist_neg_R_cut <- peaklist_neg_R[, inds_neg]

## for G, W, R take only those metabolites that have sd > 0
peaklist_pos_G_cut <- peaklist_pos_G_cut[apply(peaklist_pos_G_cut, 1, sd) > 0, ]
peaklist_neg_G_cut <- peaklist_neg_G_cut[apply(peaklist_neg_G_cut, 1, sd) > 0, ]
peaklist_pos_W_cut <- peaklist_pos_W_cut[apply(peaklist_pos_W_cut, 1, sd) > 0, ]
peaklist_neg_W_cut <- peaklist_neg_W_cut[apply(peaklist_neg_W_cut, 1, sd) > 0, ]
peaklist_pos_R_cut <- peaklist_pos_R_cut[apply(peaklist_pos_R_cut, 1, sd) > 0, ]
peaklist_neg_R_cut <- peaklist_neg_R_cut[apply(peaklist_neg_R_cut, 1, sd) > 0, ]

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

## G, W, R
stat_adj_pos_G <- statistical(as.matrix(peaklist_pos_G_cut), model = models, 
    correlation_adjust = "BH") 
stat_adj_neg_G <- statistical(as.matrix(peaklist_neg_G_cut), model = models, 
    correlation_adjust = "BH")
stat_adj_pos_W <- statistical(as.matrix(peaklist_pos_W_cut), model = models, 
    correlation_adjust = "BH") 
stat_adj_neg_W <- statistical(as.matrix(peaklist_neg_W_cut), model = models, 
    correlation_adjust = "BH")
stat_adj_pos_R <- statistical(as.matrix(peaklist_pos_R_cut), model = models, 
    correlation_adjust = "BH") 
stat_adj_neg_R <- statistical(as.matrix(peaklist_neg_R_cut), model = models, 
    correlation_adjust = "BH")

save(stat_adj_pos_G, stat_adj_pos_W, stat_adj_pos_R, file = "MetNet_strawberry_stat_adj_pos_GWR.RData")
save(stat_adj_neg_G, stat_adj_neg_W, stat_adj_neg_R, file = "MetNet_strawberry_stat_adj_neg_GWR.RData")

## how to set thresholds? use flavonoids, ellagitannins/ellagic acids and 
## hydroxycinnamic acids as benchmark
## (all negative ionization mode)
## quercetin dihexose m/z 625.1596 (rt: 372.06, 376.836, 381.684, 423.684)
tmp <- c("625.175/371.73_26443", "625.144/372.55_20414", "625.175/377.54_16692",
    "625.123/382.33_14256", "625.144/423.85_25837")
## kaempferol hexose m/z 447.0930 (rt: 302.04, 330.78, 462.9)
tmp <- c(tmp, "447.092/302.96_21354", "447.097/331.72_15038", 
    "447.091/463.22_10806", "447.053/463.83_10811")
## kaempferol glucuronide m/z 461.072 (rt: 463.09, 471.3)
tmp <- c(tmp, "461.076/463.44_10806", "461.076/471.67_24", "461.03/471.65_4290")
## quercitin glucuronide m/z 477.0690 (rt: 426.3960)
tmp <- c(tmp, "477.018/426.75_4128", "477.067/426.74_4054")
## isorhamnetin glucuronide m/z 491.0840 (rt: 476.4600)
tmp <- c(tmp, "491.083/476.87_14473")
## quercetin hexose m/z 463.0892 (rt: 428.0700)
tmp <- c(tmp, "463.089/428.34_27")
quantile(unique(stat_adj_neg[[1]][tmp, tmp]), 0.25, na.rm = TRUE) ## rf: 2.53e-06
quantile(unique(stat_adj_neg[[2]][tmp, tmp]), 0.25, na.rm = TRUE) ## clr: 0
quantile(unique(stat_adj_neg[[3]][tmp, tmp]), 0.25, na.rm = TRUE) ## aracne: 0
quantile(unique(stat_adj_neg[[4]][tmp, tmp]), 0.25, na.rm = TRUE) ## pearson: 0.154
quantile(unique(stat_adj_neg[[5]][tmp, tmp]), 0.25, na.rm = TRUE) ## pearson_part: 0.037
quantile(unique(stat_adj_neg[[6]][tmp, tmp]), 0.25, na.rm = TRUE) ## spearman: 0.140
quantile(unique(stat_adj_neg[[7]][tmp, tmp]), 0.25, na.rm = TRUE) ## spearman_part: 0.030

## bis(HHDP) glucose m/z 784.0591 (rt: 342.65) --> nf
## galloyl-HHDP-glucose m/z 633.0740 (rt: 305.64)
tmp <- "633.071/305.64_19"
## di-galloyl HHDP-glucose m/z 785.0830 (rt: 356.64)
tmp <- c(tmp, "785.084/357.04_121")
## ellagic acid deoxyhexose m/z 447.0570 (rt: 399.31, 405.96)
tmp <- c(tmp, "447.059/399.51_10667", "447.059/406.75_11197")
## ellagic acid m/z 300.999 (rt: 420.48)
tmp <- c(tmp, "300.947/420.48_9056", "300.975/420.27_61")
## galloyl-bis(HHDP)-glucose m/z 935.093 (rt: 369.77, 377.52, 404.90, 412.02)
tmp <- c(tmp, "935.09/369.81_26", "935.07/374.66_20696", "935.07/404.63_12",
    "935.09/412.7_18453")
quantile(unique(stat_adj_neg[[1]][tmp, tmp]), 0.25, na.rm = TRUE) ## rf: 2.34e-06
quantile(unique(stat_adj_neg[[2]][tmp, tmp]), 0.25, na.rm = TRUE) ## clr: 0
quantile(unique(stat_adj_neg[[3]][tmp, tmp]), 0.25, na.rm = TRUE) ## aracne: 0
quantile(unique(stat_adj_neg[[4]][tmp, tmp]), 0.25, na.rm = TRUE) ## pearson: 0.161
quantile(unique(stat_adj_neg[[5]][tmp, tmp]), 0.25, na.rm = TRUE) ## pearson_part: 0.041
quantile(unique(stat_adj_neg[[6]][tmp, tmp]), 0.25, na.rm = TRUE) ## spearman: 0.212
quantile(unique(stat_adj_neg[[7]][tmp, tmp]), 0.25, na.rm = TRUE) ## spearman_part: 0.033

## caffeic acid hexose m/z 341.1094 (rt: 297.82)
tmp <- c("341.123/346.84_18251", "341.089/324.39_206", "341.089/297.82_22926") 
## coumaric acid hexose m/z 325.0890 (rt: 311.46, 324.6)
tmp <- c(tmp, "325.094/311.58_3788", "325.035/311.66_3788", 
    "325.054/326.35_12890", "325.035/326.57_12871", "325.094/326.23_12871")
## ferulic acid hexose m/z 355.1040 (rt: 334.2, 420.0, 444.12)
tmp <- c(tmp, "355.037/444.35_2942", "355.069/324.69_206", 
    "355.073/444.41_2851", "355.105/334.47_81", "355.105/420.41_9050",
    "355.105/444.38_2851")
## sinapic acid hexose derivative m/z 385.1511 (rt: 598.344, 676.308)
tmp <- c(tmp, "385.151/598.24_24946", "385.152/675.97_26552")
quantile(unique(stat_adj_neg[[1]][tmp, tmp]), 0.25, na.rm = TRUE) ## rf: 1.59e-06
quantile(unique(stat_adj_neg[[2]][tmp, tmp]), 0.25, na.rm = TRUE) ## clr: 0.443
quantile(unique(stat_adj_neg[[3]][tmp, tmp]), 0.25, na.rm = TRUE) ## aracne: 0
quantile(unique(stat_adj_neg[[4]][tmp, tmp]), 0.25, na.rm = TRUE) ## pearson: 0.355
quantile(unique(stat_adj_neg[[5]][tmp, tmp]), 0.25, na.rm = TRUE) ## pearson_part: 0.044
quantile(unique(stat_adj_neg[[6]][tmp, tmp]), 0.25, na.rm = TRUE) ## spearman: 0.365
quantile(unique(stat_adj_neg[[7]][tmp, tmp]), 0.25, na.rm = TRUE) ## spearman_part: 0.031

## set the final values (min of the three sets)
##rf: 1.5e-06
##clr: 0.4 --> x% percentile from the complete set --> 34%
##aracne: quantile(stat_adj_neg[[3]], 0.34, na.rm = TRUE) --> 0 --> 0.01
##pearson: 0.15
##pearson_part: 0.035
##spearman: 0.14
##spearman_part: 0.03

## optional: set links within a pc_group to NA for pearson*/spearman*,
## clr/aracne/bayes/randomForest to reduce unmeaningful links
## apply function
stat_adj_pos_NA <- setLinkWithinPCgroupTo(stat_adj_pos, set_to = NaN)
stat_adj_neg_NA <- setLinkWithinPCgroupTo(stat_adj_neg, set_to = NaN)
stat_adj_pos_G_NA <- setLinkWithinPCgroupTo(stat_adj_pos_G, set_to = NaN)
stat_adj_neg_G_NA <- setLinkWithinPCgroupTo(stat_adj_neg_G, set_to = NaN)
stat_adj_pos_W_NA <- setLinkWithinPCgroupTo(stat_adj_pos_W, set_to = NaN)
stat_adj_neg_W_NA <- setLinkWithinPCgroupTo(stat_adj_neg_W, set_to = NaN)
stat_adj_pos_R_NA <- setLinkWithinPCgroupTo(stat_adj_pos_R, set_to = NaN)
stat_adj_neg_R_NA <- setLinkWithinPCgroupTo(stat_adj_neg_R, set_to = NaN)

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
table(stat_adj_neg_NA[["pearson"]][upper.tri(stat_adj_neg_NA[["pearson"]])] > 0.15)
table(stat_adj_pos_NA[["pearson"]][upper.tri(stat_adj_pos_NA[["pearson"]])] > 0.15)
table(stat_adj_neg_NA[["pearson_partial"]][upper.tri(stat_adj_neg_NA[["pearson_partial"]])] > 0.035)
table(stat_adj_pos_NA[["pearson_partial"]][upper.tri(stat_adj_pos_NA[["pearson_partial"]])] > 0.035)
table(stat_adj_neg_NA[["spearman"]][upper.tri(stat_adj_neg_NA[["spearman"]])] > 0.14)
table(stat_adj_pos_NA[["spearman"]][upper.tri(stat_adj_pos_NA[["spearman"]])] > 0.14)
table(stat_adj_neg_NA[["spearman_partial"]][upper.tri(stat_adj_neg_NA[["spearman_partial"]])] > 0.03)
table(stat_adj_pos_NA[["spearman_partial"]][upper.tri(stat_adj_pos_NA[["spearman_partial"]])] > 0.03)
table(stat_adj_neg_NA[["clr"]][upper.tri(stat_adj_neg_NA[["clr"]])] > 0.4)
table(stat_adj_pos_NA[["clr"]][upper.tri(stat_adj_pos_NA[["clr"]])] > 0.4)
table(stat_adj_neg_NA[["aracne"]][upper.tri(stat_adj_neg_NA[["aracne"]])] > 0.01)
table(stat_adj_pos_NA[["aracne"]][upper.tri(stat_adj_pos_NA[["aracne"]])] > 0.01)
table(stat_adj_neg_NA[["randomForest"]][upper.tri(stat_adj_neg_NA[["randomForest"]])] > 1.5e-06)
table(stat_adj_pos_NA[["randomForest"]][upper.tri(stat_adj_pos_NA[["randomForest"]])] > 1.5e-06)

args <- list("pearson" = 0.15, "pearson_partial" = 0.035, "spearman" = 0.14,
    "spearman_partial" = 0.03, "clr" = 0.4, "aracne" = 0.01, 
    "randomForest" = 1.5e-06, threshold = 1)
stat_adj_pos_thr <- threshold(statistical = stat_adj_pos_NA, type = "threshold",
    args = args)
stat_adj_neg_thr <- threshold(statistical = stat_adj_neg_NA, type = "threshold",
    args = args)
stat_adj_pos_G_thr <- threshold(statistical = stat_adj_pos_G_NA, type = "threshold",
    args = args)
stat_adj_neg_G_thr <- threshold(statistical = stat_adj_neg_G_NA, type = "threshold",
    args = args)
stat_adj_pos_W_thr <- threshold(statistical = stat_adj_pos_W_NA, type = "threshold",
    args = args)
stat_adj_neg_W_thr <- threshold(statistical = stat_adj_neg_W_NA, type = "threshold",
    args = args)
stat_adj_pos_R_thr <- threshold(statistical = stat_adj_pos_R_NA, type = "threshold",
    args = args)
stat_adj_neg_R_thr <- threshold(statistical = stat_adj_neg_R_NA, type = "threshold",
    args = args)


## type = "top1" 
args_top <- list(n = 500000)
stat_adj_pos_top1 <- threshold(statistical = stat_adj_pos_NA, type = "top1",
    args = args_top)
stat_adj_neg_top1 <- threshold(statistical = stat_adj_neg_NA, type = "top1",
    args = args_top)
stat_adj_pos_G_top1 <- threshold(statistical = stat_adj_pos_G_NA, type = "top1",
    args = args_top)
stat_adj_neg_G_top1 <- threshold(statistical = stat_adj_neg_G_NA, type = "top1",
    args = args_top)
stat_adj_pos_W_top1 <- threshold(statistical = stat_adj_pos_W_NA, type = "top1",
    args = args_top)
stat_adj_neg_W_top1 <- threshold(statistical = stat_adj_neg_W_NA, type = "top1",
    args = args_top)
stat_adj_pos_R_top1 <- threshold(statistical = stat_adj_pos_R_NA, type = "top1",
    args = args_top)
stat_adj_neg_R_top1 <- threshold(statistical = stat_adj_neg_R_NA, type = "top1",
    args = args_top)
 
## type = "top2"
stat_adj_pos_top2 <- threshold(statistical = stat_adj_pos_NA, type = "top2",
    args = args_top)
stat_adj_neg_top2 <- threshold(statistical = stat_adj_neg_NA, type = "top2",
    args = args_top)
stat_adj_pos_G_top2 <- threshold(statistical = stat_adj_pos_G_NA, type = "top2",
    args = args_top)
stat_adj_neg_G_top2 <- threshold(statistical = stat_adj_neg_G_NA, type = "top2",
    args = args_top)
stat_adj_pos_W_top2 <- threshold(statistical = stat_adj_pos_W_NA, type = "top2",
    args = args_top)
stat_adj_neg_W_top2 <- threshold(statistical = stat_adj_neg_W_NA, type = "top2",
    args = args_top)
stat_adj_pos_R_top2 <- threshold(statistical = stat_adj_pos_R_NA, type = "top2",
    args = args_top)
stat_adj_neg_R_top2 <- threshold(statistical = stat_adj_neg_R_NA, type = "top2",
    args = args_top)

## type = "mean"
stat_adj_pos_mean <- threshold(statistical = stat_adj_pos_NA, type = "mean",
    args = args_top)
stat_adj_neg_mean <- threshold(statistical = stat_adj_neg_NA, type = "mean",
    args = args_top)
stat_adj_pos_G_mean <- threshold(statistical = stat_adj_pos_G_NA, type = "mean",
    args = args_top)
stat_adj_neg_G_mean <- threshold(statistical = stat_adj_neg_G_NA, type = "mean",
    args = args_top)
stat_adj_pos_W_mean <- threshold(statistical = stat_adj_pos_W_NA, type = "mean",
    args = args_top)
stat_adj_neg_W_mean <- threshold(statistical = stat_adj_neg_W_NA, type = "mean",
    args = args_top)
stat_adj_pos_R_mean <- threshold(statistical = stat_adj_pos_R_NA, type = "mean",
    args = args_top)
stat_adj_neg_R_mean <- threshold(statistical = stat_adj_neg_R_NA, type = "mean",
    args = args_top)

save(stat_adj_pos_thr, stat_adj_pos_top1, stat_adj_pos_top2, stat_adj_pos_mean, 
    file = "MetNet_strawberry_stat_adj_thr_pos.RData")
save(stat_adj_neg_thr, stat_adj_neg_top1, stat_adj_neg_top2, stat_adj_neg_mean, 
    file = "MetNet_strawberry_stat_adj_thr_neg.RData")
save(stat_adj_pos_G_thr, stat_adj_pos_G_top1, stat_adj_pos_G_top2, stat_adj_pos_G_mean, 
    stat_adj_pos_W_thr, stat_adj_pos_W_top1, stat_adj_pos_W_top2, stat_adj_pos_W_mean, 
    stat_adj_pos_R_thr, stat_adj_pos_R_top1, stat_adj_pos_R_top2, stat_adj_pos_R_mean, 
    file = "MetNet_strawberry_stat_adj_thr_pos_GWR.RData")
save(stat_adj_neg_G_thr, stat_adj_neg_G_top1, stat_adj_neg_G_top2, stat_adj_neg_G_mean, 
    stat_adj_neg_W_thr, stat_adj_neg_W_top1, stat_adj_neg_W_top2, stat_adj_neg_W_mean,
    stat_adj_neg_R_thr, stat_adj_neg_R_top1, stat_adj_neg_R_top2, stat_adj_neg_R_mean, 
    file = "MetNet_strawberry_stat_adj_thr_neg_GWR.RData")


## set links between features belonging to different pcgroups in struct_adj to
## 0 for links that contain the "adduct_" pattern (type = "inter")
## set links between features belonging to the same pcgroup in struct_adj to 0
## for links that contain the "transf_" pattern (type = "intra")
struct_adj_neg_mod <- setLinkPCgroupTo0(struct_adj_neg, type = "inter", 
    pattern = "adduct_")
struct_adj_pos_mod <- setLinkPCgroupTo0(struct_adj_pos, type = "inter",
    pattern = "adduct_")

## set links within pcgroups that are isotope peaks to 0, remove all rows and cols 
## that are isotopes

##struct_adj_neg_mod <- setLinkPCgroupTo0(struct_adj_neg_mod, type = "intra", 
##    pattern = "transf_")
##struct_adj_pos_mod <- setLinkPCgroupTo0(struct_adj_pos_mod, type = "intra",
##    pattern = "transf_")

## set all intra pcgroup links in stat_adj to 1 in order to make sure that all 
## links from struct_adj for intra pcgroups links are taken 
stat_adj_pos_thr_mod <- setLinkWithinPCgroupTo(stat_adj_pos_thr, set_to = 1)
stat_adj_neg_thr_mod <- setLinkWithinPCgroupTo(stat_adj_neg_thr, set_to = 1)
stat_adj_pos_thr_G_mod <- setLinkWithinPCgroupTo(stat_adj_pos_G_thr, set_to = 1)
stat_adj_neg_thr_G_mod <- setLinkWithinPCgroupTo(stat_adj_neg_G_thr, set_to = 1)
stat_adj_pos_thr_W_mod <- setLinkWithinPCgroupTo(stat_adj_pos_W_thr, set_to = 1)
stat_adj_neg_thr_W_mod <- setLinkWithinPCgroupTo(stat_adj_neg_W_thr, set_to = 1)
stat_adj_pos_thr_R_mod <- setLinkWithinPCgroupTo(stat_adj_pos_R_thr, set_to = 1)
stat_adj_neg_thr_R_mod <- setLinkWithinPCgroupTo(stat_adj_neg_R_thr, set_to = 1)

## truncate the structural adjacency matrices for GWR
struct_adj_pos_mod_G <- struct_adj_pos_mod
struct_adj_pos_mod_G[[1]] <- struct_adj_pos_mod[[1]][rownames(stat_adj_pos_thr_G_mod), rownames(stat_adj_pos_thr_G_mod)]
struct_adj_pos_mod_G[[2]] <- struct_adj_pos_mod[[2]][rownames(stat_adj_pos_thr_G_mod), rownames(stat_adj_pos_thr_G_mod)]
struct_adj_neg_mod_G <- struct_adj_neg_mod
struct_adj_neg_mod_G[[1]] <- struct_adj_neg_mod[[1]][rownames(stat_adj_neg_thr_G_mod), rownames(stat_adj_neg_thr_G_mod)]
struct_adj_neg_mod_G[[2]] <- struct_adj_neg_mod[[2]][rownames(stat_adj_neg_thr_G_mod), rownames(stat_adj_neg_thr_G_mod)]
struct_adj_pos_mod_W <- struct_adj_pos_mod
struct_adj_pos_mod_W[[1]] <- struct_adj_pos_mod[[1]][rownames(stat_adj_pos_thr_W_mod), rownames(stat_adj_pos_thr_W_mod)]
struct_adj_pos_mod_W[[2]] <- struct_adj_pos_mod[[2]][rownames(stat_adj_pos_thr_W_mod), rownames(stat_adj_pos_thr_W_mod)]
struct_adj_neg_mod_W <- struct_adj_neg_mod
struct_adj_neg_mod_W[[1]] <- struct_adj_neg_mod[[1]][rownames(stat_adj_neg_thr_W_mod), rownames(stat_adj_neg_thr_W_mod)]
struct_adj_neg_mod_W[[2]] <- struct_adj_neg_mod[[2]][rownames(stat_adj_neg_thr_W_mod), rownames(stat_adj_neg_thr_W_mod)]
struct_adj_pos_mod_R <- struct_adj_pos_mod
struct_adj_pos_mod_R[[1]] <- struct_adj_pos_mod[[1]][rownames(stat_adj_pos_thr_R_mod), rownames(stat_adj_pos_thr_R_mod)]
struct_adj_pos_mod_R[[2]] <- struct_adj_pos_mod[[2]][rownames(stat_adj_pos_thr_R_mod), rownames(stat_adj_pos_thr_R_mod)]
struct_adj_neg_mod_R <- struct_adj_neg_mod
struct_adj_neg_mod_R[[1]] <- struct_adj_neg_mod[[1]][rownames(stat_adj_neg_thr_R_mod), rownames(stat_adj_neg_thr_R_mod)]
struct_adj_neg_mod_R[[2]] <- struct_adj_neg_mod[[2]][rownames(stat_adj_neg_thr_R_mod), rownames(stat_adj_neg_thr_R_mod)]

## use function combine to combine the structural and statistical information
cons_adj_pos <- combine(structural = struct_adj_pos_mod, 
    statistical = stat_adj_pos_thr_mod)
cons_adj_neg <- combine(structural = struct_adj_neg_mod, 
    statistical = stat_adj_neg_thr_mod)
cons_adj_pos_G <- combine(structural = struct_adj_pos_mod_G, 
    statistical = stat_adj_pos_thr_G_mod)
cons_adj_neg_G <- combine(structural = struct_adj_neg_mod_G, 
    statistical = stat_adj_neg_thr_G_mod)
cons_adj_pos_W <- combine(structural = struct_adj_pos_mod_W, 
    statistical = stat_adj_pos_thr_W_mod)
cons_adj_neg_W <- combine(structural = struct_adj_neg_mod_W, 
    statistical = stat_adj_neg_thr_W_mod)
cons_adj_pos_R <- combine(structural = struct_adj_pos_mod_R, 
    statistical = stat_adj_pos_thr_R_mod)
cons_adj_neg_R <- combine(structural = struct_adj_neg_mod_R, 
    statistical = stat_adj_neg_thr_R_mod)

save(cons_adj_pos, file = "MetNet_strawberry_cons_adj_pos.RData")
save(cons_adj_neg, file = "MetNet_strawberry_cons_adj_neg.RData")
save(cons_adj_pos_G, cons_adj_pos_W, cons_adj_pos_R, file = "MetNet_strawberry_cons_adj_pos_GWR.RData")
save(cons_adj_neg_G, cons_adj_neg_W, cons_adj_neg_R, file = "MetNet_strawberry_cons_adj_neg_GWR.RData")

## only retain vertices that link to at least two other mass features
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

pdf("net_coms_comp_pos.pdf")
plot(g_pos, vertex.label.cex = 0.1, vertex.size = 0.1, edge.arrow.width = 0.1,
    edge.arrow.size = 0.1)
dev.off()

pdf("net_coms_comp_neg.pdf")
plot(g_neg, vertex.label.cex = 0.1, vertex.size = 0.1, edge.arrow.width = 0.1,
     edge.arrow.size = 0.1)
dev.off()

g_pos <- graph_from_adjacency_matrix(cons_adj_pos_G[[1]], mode = "directed")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg_G[[1]], mode = "directed")
comp_g_pos <- components(g_pos)
comp_g_neg <- components(g_neg)
inds_cut_pos <- comp_g_pos$membership %in% which(comp_g_pos$csize > 2)
inds_cut_neg <- comp_g_neg$membership %in% which(comp_g_neg$csize > 2)
cons_adj_pos_cut_G <- cons_adj_pos_G
cons_adj_neg_cut_G <- cons_adj_neg_G
cons_adj_pos_cut_G[[1]] <- cons_adj_pos_G[[1]][inds_cut_pos, inds_cut_pos]
cons_adj_pos_cut_G[[2]] <- cons_adj_pos_G[[2]][inds_cut_pos, inds_cut_pos]
cons_adj_neg_cut_G[[1]] <- cons_adj_neg_G[[1]][inds_cut_neg, inds_cut_neg]
cons_adj_neg_cut_G[[2]] <- cons_adj_neg_G[[2]][inds_cut_neg, inds_cut_neg]

g_pos <- graph_from_adjacency_matrix(cons_adj_pos_W[[1]], mode = "directed")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg_W[[1]], mode = "directed")
comp_g_pos <- components(g_pos)
comp_g_neg <- components(g_neg)
inds_cut_pos <- comp_g_pos$membership %in% which(comp_g_pos$csize > 2)
inds_cut_neg <- comp_g_neg$membership %in% which(comp_g_neg$csize > 2)
cons_adj_pos_cut_W <- cons_adj_pos_W
cons_adj_neg_cut_W <- cons_adj_neg_W
cons_adj_pos_cut_W[[1]] <- cons_adj_pos_W[[1]][inds_cut_pos, inds_cut_pos]
cons_adj_pos_cut_W[[2]] <- cons_adj_pos_W[[2]][inds_cut_pos, inds_cut_pos]
cons_adj_neg_cut_W[[1]] <- cons_adj_neg_W[[1]][inds_cut_neg, inds_cut_neg]
cons_adj_neg_cut_W[[2]] <- cons_adj_neg_W[[2]][inds_cut_neg, inds_cut_neg]

g_pos <- graph_from_adjacency_matrix(cons_adj_pos_R[[1]], mode = "directed")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg_R[[1]], mode = "directed")
comp_g_pos <- components(g_pos)
comp_g_neg <- components(g_neg)
inds_cut_pos <- comp_g_pos$membership %in% which(comp_g_pos$csize > 2)
inds_cut_neg <- comp_g_neg$membership %in% which(comp_g_neg$csize > 2)
cons_adj_pos_cut_R <- cons_adj_pos_R
cons_adj_neg_cut_R <- cons_adj_neg_R
cons_adj_pos_cut_R[[1]] <- cons_adj_pos_R[[1]][inds_cut_pos, inds_cut_pos]
cons_adj_pos_cut_R[[2]] <- cons_adj_pos_R[[2]][inds_cut_pos, inds_cut_pos]
cons_adj_neg_cut_R[[1]] <- cons_adj_neg_R[[1]][inds_cut_neg, inds_cut_neg]
cons_adj_neg_cut_R[[2]] <- cons_adj_neg_R[[2]][inds_cut_neg, inds_cut_neg]


## remove edges between M1+H and M2+Na, M3+H and M4+K, in general when there is 
## a link between an adduct and a molecular ion
## node with outgoing transf link should not have outgoing adduct link when
## the adduct-linking vertex has a ingoing adduct link
## node with ingoing transf link should not have ingoing adduct link when the 
## adduct-linking vertex has an outgoing adduct link
cons_adj_pos_cut_rem <- removeFalseLinksAdducts(cons_adj_pos_cut)
cons_adj_neg_cut_rem <- removeFalseLinksAdducts(cons_adj_neg_cut)
cons_adj_pos_cut_rem_G <- removeFalseLinksAdducts(cons_adj_pos_cut_G)
cons_adj_neg_cut_rem_G <- removeFalseLinksAdducts(cons_adj_neg_cut_G)
cons_adj_pos_cut_rem_W <- removeFalseLinksAdducts(cons_adj_pos_cut_W)
cons_adj_neg_cut_rem_W <- removeFalseLinksAdducts(cons_adj_neg_cut_W)
cons_adj_pos_cut_rem_R <- removeFalseLinksAdducts(cons_adj_pos_cut_R)
cons_adj_neg_cut_rem_R <- removeFalseLinksAdducts(cons_adj_neg_cut_R)

## Remove in- and outgoing links based on their relation to the molecular ion
## The function removeFalseLinksCircular deletes ingoing and outgoing links of 
## isotopes/adducts to M if there is no circular relation from the molecular 
## ion M to M+isotope/adduct to M+transf to M+iso/adduct+transformation
cons_adj_pos_cut_rem <- removeFalseLinksCircular(cons_adj_pos_cut_rem)
cons_adj_neg_cut_rem <- removeFalseLinksCircular(cons_adj_neg_cut_rem)
cons_adj_pos_cut_rem_G <- removeFalseLinksCircular(cons_adj_pos_cut_rem_G)
cons_adj_neg_cut_rem_G <- removeFalseLinksCircular(cons_adj_neg_cut_rem_G)
cons_adj_pos_cut_rem_W <- removeFalseLinksCircular(cons_adj_pos_cut_rem_W)
cons_adj_neg_cut_rem_W <- removeFalseLinksCircular(cons_adj_neg_cut_rem_W)
cons_adj_pos_cut_rem_R <- removeFalseLinksCircular(cons_adj_pos_cut_rem_R)
cons_adj_neg_cut_rem_R <- removeFalseLinksCircular(cons_adj_neg_cut_rem_R)

## write to xml file
g_pos <- graph_from_adjacency_matrix(cons_adj_pos_cut[[1]], mode = "directed")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg_cut[[1]], mode = "directed")
write_graph(g_pos, file = "cons_adj_pos_cut_graphml.xml", format = "graphml")
write_graph(g_neg, file = "cons_adj_neg_cut_graphml.xml", format = "graphml")
g_pos <- graph_from_adjacency_matrix(cons_adj_pos_cut_rem[[1]], mode = "directed")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg_cut_rem[[1]], mode = "directed")
write_graph(g_pos, file = "cons_adj_pos_cut_rem_graphml.xml", format = "graphml")
write_graph(g_neg, file = "cons_adj_neg_cut_rem_graphml.xml", format = "graphml")

g_pos <- graph_from_adjacency_matrix(cons_adj_pos_cut_G[[1]], mode = "directed")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg_cut_G[[1]], mode = "directed")
write_graph(g_pos, file = "cons_adj_pos_cut_G_graphml.xml", format = "graphml")
write_graph(g_neg, file = "cons_adj_neg_cut_G_graphml.xml", format = "graphml")
g_pos <- graph_from_adjacency_matrix(cons_adj_pos_cut_rem_G[[1]], mode = "directed")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg_cut_rem_G[[1]], mode = "directed")
write_graph(g_pos, file = "cons_adj_pos_cut_rem_G_graphml.xml", format = "graphml")
write_graph(g_neg, file = "cons_adj_neg_cut_rem_G_graphml.xml", format = "graphml")

g_pos <- graph_from_adjacency_matrix(cons_adj_pos_cut_W[[1]], mode = "directed")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg_cut_W[[1]], mode = "directed")
write_graph(g_pos, file = "cons_adj_pos_cut_W_graphml.xml", format = "graphml")
write_graph(g_neg, file = "cons_adj_neg_cut_W_graphml.xml", format = "graphml")
g_pos <- graph_from_adjacency_matrix(cons_adj_pos_cut_rem_W[[1]], mode = "directed")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg_cut_rem_W[[1]], mode = "directed")
write_graph(g_pos, file = "cons_adj_pos_cut_rem_W_graphml.xml", format = "graphml")
write_graph(g_neg, file = "cons_adj_neg_cut_rem_W_graphml.xml", format = "graphml")

g_pos <- graph_from_adjacency_matrix(cons_adj_pos_cut_R[[1]], mode = "directed")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg_cut_R[[1]], mode = "directed")
write_graph(g_pos, file = "cons_adj_pos_cut_R_graphml.xml", format = "graphml")
write_graph(g_neg, file = "cons_adj_neg_cut_R_graphml.xml", format = "graphml")
g_pos <- graph_from_adjacency_matrix(cons_adj_pos_cut_rem_R[[1]], mode = "directed")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg_cut_rem_R[[1]], mode = "directed")
write_graph(g_pos, file = "cons_adj_pos_cut_rem_R_graphml.xml", format = "graphml")
write_graph(g_neg, file = "cons_adj_neg_cut_rem_R_graphml.xml", format = "graphml")

## go through components and only take those components that have > 1 
## different pcgroup
filter_pc <- function(cons) {
    inds_cut <- rep(FALSE, nrow(cons[[1]]))
    g <- graph_from_adjacency_matrix(cons[[1]], mode = "directed")
    comp <- components(g)
    
    for (i in 1:comp$no) {
        ms_i <- comp$membership[comp$membership == i]
        pc_i <- unlist(lapply(strsplit(names(ms_i), split = "_"), "[", 2))
        pc_i <- unique(pc_i)
        if (length(pc_i) > 1) {
            inds_cut[comp$membership == i] <- TRUE    
        }
    }
    cons_cut <- cons
    cons_cut[[1]] <- cons_cut[[1]][inds_cut, inds_cut]
    cons_cut[[2]] <- cons_cut[[2]][inds_cut, inds_cut]
    return(cons_cut)
    
}

cons_adj_pos_cut_rem_cut <- filter_pc(cons_adj_pos_cut_rem)
cons_adj_neg_cut_rem_cut <- filter_pc(cons_adj_neg_cut_rem)
cons_adj_pos_cut_rem_cut_G <- filter_pc(cons_adj_pos_cut_rem_G)
cons_adj_neg_cut_rem_cut_G <- filter_pc(cons_adj_neg_cut_rem_G)
cons_adj_pos_cut_rem_cut_W <- filter_pc(cons_adj_pos_cut_rem_W)
cons_adj_neg_cut_rem_cut_W <- filter_pc(cons_adj_neg_cut_rem_W)
cons_adj_pos_cut_rem_cut_R <- filter_pc(cons_adj_pos_cut_rem_R)
cons_adj_neg_cut_rem_cut_R <- filter_pc(cons_adj_neg_cut_rem_R)

g_pos <- graph_from_adjacency_matrix(cons_adj_pos_cut_rem_cut[[1]], mode = "directed")
g_neg <- graph_from_adjacency_matrix(cons_adj_neg_cut_rem_cut[[1]], mode = "directed")
g_pos_G <- graph_from_adjacency_matrix(cons_adj_pos_cut_rem_cut_G[[1]], mode = "directed")
g_neg_G <- graph_from_adjacency_matrix(cons_adj_neg_cut_rem_cut_G[[1]], mode = "directed")
g_pos_W <- graph_from_adjacency_matrix(cons_adj_pos_cut_rem_cut_W[[1]], mode = "directed")
g_neg_W <- graph_from_adjacency_matrix(cons_adj_neg_cut_rem_cut_W[[1]], mode = "directed")
g_pos_R <- graph_from_adjacency_matrix(cons_adj_pos_cut_rem_cut_R[[1]], mode = "directed")
g_neg_R <- graph_from_adjacency_matrix(cons_adj_neg_cut_rem_cut_R[[1]], mode = "directed")

write_graph(g_pos, file = "cons_adj_pos_cut_rem_cut_graphml.xml", format = "graphml")
write_graph(g_neg, file = "cons_adj_neg_cut_rem_cut_graphml.xml", format = "graphml")
write_graph(g_pos_G, file = "cons_adj_pos_cut_rem_cut_G_graphml.xml", format = "graphml")
write_graph(g_neg_G, file = "cons_adj_neg_cut_rem_cut_G_graphml.xml", format = "graphml")
write_graph(g_pos_W, file = "cons_adj_pos_cut_rem_cut_W_graphml.xml", format = "graphml")
write_graph(g_neg_W, file = "cons_adj_neg_cut_rem_cut_W_graphml.xml", format = "graphml")
write_graph(g_pos_R, file = "cons_adj_pos_cut_rem_cut_R_graphml.xml", format = "graphml")
write_graph(g_neg_R, file = "cons_adj_neg_cut_rem_cut_R_graphml.xml", format = "graphml")

save(cons_adj_pos_cut_rem_cut, cons_adj_neg_cut_rem_cut, 
    cons_adj_pos_cut_rem_cut_G, cons_adj_neg_cut_rem_cut_G,
    cons_adj_pos_cut_rem_cut_W, cons_adj_neg_cut_rem_cut_W,
    cons_adj_pos_cut_rem_cut_R, cons_adj_neg_cut_rem_cut_R, 
    file = "cons_adj_cut_rem_cut.RData")

## write peaklists
write.table(peaklist_pos[rownames(cons_adj_pos_cut_rem_cut[[1]]), ], 
    file = "cons_adj_pos_cut_rem_cut_peaklist.txt", sep = "\t", dec = ".", 
    quote = FALSE, row.names = TRUE)
write.table(peaklist_neg[rownames(cons_adj_neg_cut_rem_cut[[1]]), ], 
    file = "cons_adj_neg_cut_rem_cut_peaklist.txt", sep = "\t", dec = ".", 
    quote = FALSE, row.names = TRUE)
write.table(peaklist_pos_W[rownames(cons_adj_pos_cut_rem_cut_W[[1]]), ], 
    file = "cons_adj_pos_cut_rem_cut_W_peaklist.txt", sep = "\t", dec = ".", 
    quote = FALSE, row.names = TRUE)
write.table(peaklist_neg_W[rownames(cons_adj_neg_cut_rem_cut_W[[1]]), ], 
    file = "cons_adj_neg_cut_rem_cut_W_peaklist.txt", sep = "\t", dec = ".", 
    quote = FALSE, row.names = TRUE)
write.table(peaklist_pos_G[rownames(cons_adj_pos_cut_rem_cut_G[[1]]), ], 
    file = "cons_adj_pos_cut_rem_cut_G_peaklist.txt", sep = "\t", dec = ".", 
    quote = FALSE, row.names = TRUE)
write.table(peaklist_neg_G[rownames(cons_adj_neg_cut_rem_cut_G[[1]]), ], 
    file = "cons_adj_neg_cut_rem_cut_G_peaklist.txt", sep = "\t", dec = ".", 
    quote = FALSE, row.names = TRUE)
write.table(peaklist_pos_R[rownames(cons_adj_pos_cut_rem_cut_R[[1]]), ], 
    file = "cons_adj_pos_cut_rem_cut_R_peaklist.txt", sep = "\t", dec = ".", 
    quote = FALSE, row.names = TRUE)
write.table(peaklist_neg_R[rownames(cons_adj_neg_cut_rem_cut_R[[1]]), ], 
    file = "cons_adj_neg_cut_rem_cut_R_peaklist.txt", sep = "\t", dec = ".", 
    quote = FALSE, row.names = TRUE)



