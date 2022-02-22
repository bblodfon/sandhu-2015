# take 5-6 min to run on all files
library(dplyr)
library(readr)
library(stringr)
library(AgiMicroRna)

source(file = 'mRNA_funs.R')

# AFE = Agilent Feature Extraction
# I have splitted the raw data txt files to 2 separate directories:
mRNA_dir  = '/media/john/PRIVATE_USB/mRNA_81_samples'
miRNA_dir = '/media/john/PRIVATE_USB/miRNA_83_samples'

dir.exists(mRNA_dir)
dir.exists(miRNA_dir)

mRNA_files  = list.files(path = mRNA_dir, full.names = T)
miRNA_files = list.files(path = miRNA_dir, full.names = T)

# DIY read a file
# miRNA_data = readr::read_tsv(file = miRNA_files[1], skip = 9, col_select = 2:last_col())
# mRNA_data  = readr::read_tsv(file = mRNA_files[1], skip = 9, col_select = 1:last_col())

mRNA_filenames  = stringr::str_match(base::basename(mRNA_files), "_(.*).txt")[,2]
miRNA_filenames = stringr::str_match(base::basename(miRNA_files), "_(.*).txt")[,2]

# take samples for testing
#mRNA_files  = mRNA_files[1:3]
#miRNA_files = miRNA_files[1:3]
#mRNA_filenames = mRNA_filenames[1:3]
#miRNA_filenames = miRNA_filenames[1:3]

########
# mRNA #
########

# take a sample
mRNA_targets = cbind.data.frame(FileName = mRNA_files,
  Treatment = rep('A', length(mRNA_files)), # random fill ins
  GErep = 1:length(mRNA_files), # random
  Subject = 1:length(mRNA_files)) # random
rownames(mRNA_targets) = mRNA_filenames

data_mRNA = AgiMicroRna::readmRnaAFE(mRNA_targets, verbose = T)

# Normalize and background correct (only mean signal)
dd_mRNA = meanSignalMrna(data_mRNA, half = T, verbose = T)
ddNORM_mRNA = meanSignal_Normalization(dd_mRNA, NORMmethod = 'quantile', targets = mRNA_targets, verbose = TRUE)
dim(ddNORM_mRNA$meanS)

# filter data for quality control
data_mRNA_filtered = AgiMicroRna::filterMicroRna(
  ddNORM = ddNORM_mRNA,
  dd = data_mRNA,
  control = TRUE, # remove controls
  IsGeneDetected = FALSE, # no 'IsGeneDetected' feature column!
  wellaboveNEG = TRUE, # =TRUE, removes more genes (more strict)!
  limIsGeneDetected = 50, # mentioned in the paper?
  limNEG = 95, # from paper?!
  makePLOT = FALSE,
  targets = mRNA_targets,
  verbose = TRUE,
  writeout = FALSE)
dim(data_mRNA_filtered$meanS)

# save miRNA expression matrix (Quantile normalized and filtered)
mRNA_expr_mat = data_mRNA_filtered$meanS
colnames(mRNA_expr_mat) = colnames(data_mRNA$meanS) # sample names
rownames(mRNA_expr_mat) = data_mRNA_filtered$genes$GeneName # gene names
saveRDS(object = mRNA_expr_mat, file = 'matrices/mRNA_expr_mat.rds')

# RMA (Robust Multiarray Average) seems like a nice alternative for data pre-processing!
# returned values are log2-transformed already
data_mRNA_rma = AgiMicroRna::rmaMicroRna(data_mRNA, normalize = TRUE, background = TRUE)
all(data_mRNA_rma$TGS == data_mRNA_rma$meanS)
dim(data_mRNA_rma)

# filter data for quality control
data_mRNA_rma_filtered = AgiMicroRna::filterMicroRna(ddNORM = data_mRNA_rma,
  dd = data_mRNA,
  control = TRUE, # remove controls
  IsGeneDetected = FALSE, # no 'IsGeneDetected' feature column!
  wellaboveNEG = TRUE, # =TRUE, removes more genes (more strict)!
  limIsGeneDetected = 50, # mentioned in the paper?
  limNEG = 95, # from paper?!
  makePLOT = FALSE,
  targets = mRNA_targets,
  verbose = TRUE,
  writeout = FALSE)
dim(data_mRNA_rma_filtered$meanS)

# save mRNA expression matrix (RMA normalized and filtered)
mRNA_rma_expr_mat = data_mRNA_rma_filtered$meanS
colnames(mRNA_rma_expr_mat) = colnames(data_mRNA$meanS) # sample names
rownames(mRNA_rma_expr_mat) = data_mRNA_rma_filtered$genes$GeneName # gene names
saveRDS(object = mRNA_rma_expr_mat, file = 'matrices/mRNA_rma_expr_mat.rds')

#########
# miRNA #
#########
miRNA_targets = cbind.data.frame(FileName = miRNA_files,
  Treatment = rep('A', length(miRNA_files)), # random fill ins
  GErep = 1:length(miRNA_files), # random
  Subject = 1:length(miRNA_files)) # random
rownames(miRNA_targets) = miRNA_filenames
data_miRNA = AgiMicroRna::readMicroRnaAFE(miRNA_targets, verbose = T)

# see some stats
if (FALSE) {
  boxplotMicroRna(log2(data_miRNA$meanS), maintitle = 'Log2 Mean Signal (boxplot)')
  plotDensityMicroRna(log2(data_miRNA$meanS), maintitle = 'Log2 Mean Signal (density)')
  # RLE: If the majority of the spots are expected not to be differentially expressed,
  # the boxplots should be centered on zero and all of them with approximately the
  # same dispersion.
  RleMicroRna(log2(data_miRNA$meanS), maintitle = 'log2 Mean Signal Log2 Relative Expression')
  hierclusMicroRna(log2(data_miRNA$meanS),targets.micro$GErep, methdis = "euclidean",
    methclu = "complete", sel = TRUE, size = 100) # 100 highest variance genes for HC

  # takes time to show the MA plot
  ddaux   = data_miRNA
  ddaux$G = log2(data_miRNA$meanS)
  mvaMicroRna(ddaux, maintitle = 'Log2 Mean Signal (MA plot)', verbose = T)
  rm(ddaux)

  # A lower CV median indicates a better array reproducibility
  cvArray(data_miRNA, foreground = "MeanSignal", targets = miRNA_targets, verbose = T)
}

# Normalize and background correct (not done in the paper for miRNA data)
if (FALSE) {
  # AFFECTS ONLY `data_miRNA$TPS`!
  dd_miRNA     = AgiMicroRna::tgsMicroRna(data_miRNA, half = T, makePLOT = F, verbose = T) # not log2 values yet
  dim(dd_miRNA)
  ddNORM_miRNA = AgiMicroRna::tgsNormalization(dd_miRNA, NORMmethod = 'quantile', makePLOTpre = FALSE,
    makePLOTpost = FALSE, targets = miRNA_targets, verbose = TRUE) # log2 values

  RleMicroRna(ddNORM_miRNA$TGS, "RLE TGS NORM", "blue")

  # Note: the other matrices stay the same - they reduce in size because duplicate genes are removed
  all(ddNORM_miRNA$TPS == dd_miRNA$TPS)
  all(ddNORM_miRNA$meanS == dd_miRNA$meanS)
  all(ddNORM_miRNA$procS == dd_miRNA$procS)
}

# RMA (Robust Multiarray Average) pre-processing
# Background corrects and normalizes `meanS` and gives it back as `TGS`!!!
data_miRNA_rma = AgiMicroRna::rmaMicroRna(data_miRNA, normalize = TRUE, background = TRUE)
dim(data_miRNA_rma)
all(data_miRNA_rma$TGS == data_miRNA_rma$meanS) # same!

# filter data for quality control
data_miRNA_rma_filtered = AgiMicroRna::filterMicroRna(ddNORM = data_miRNA_rma,
  dd = data_miRNA,
  control = TRUE, # remove controls
  IsGeneDetected = TRUE, # necessary
  wellaboveNEG = FALSE, # =TRUE, removes more genes (more strict)!
  limIsGeneDetected = 50, # mentioned in the paper?
  limNEG = 25, # maybe!
  makePLOT = FALSE,
  targets = miRNA_targets,
  verbose = TRUE,
  writeout = FALSE)

# save miRNA expression matrix
miRNA_expr_mat = data_miRNA_rma_filtered$meanS
colnames(miRNA_expr_mat) = colnames(data_miRNA$meanS) # sample names
rownames(miRNA_expr_mat) = data_miRNA_rma_filtered$genes$GeneName # gene names
saveRDS(object = miRNA_expr_mat, file = 'matrices/miRNA_expr_mat.rds')

# Write session info
writeLines(capture.output(xfun::session_info()), 'session_info.txt')
