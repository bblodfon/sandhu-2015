# mRNA Quantile Normalized
mRNA_expr_mat = readRDS(file = 'matrices/mRNA_expr_mat.rds')
dim(mRNA_expr_mat) # (34612, 81) #

# get std of all mRNAs
mRNA_stds = apply(mRNA_expr_mat, 1, sd)
# filter to most variable ones
var_mRNAs = mRNA_stds[mRNA_stds > 0.8] # hard-coded
length(var_mRNAs) # 6780

mRNA_expr_mat_red = mRNA_expr_mat[names(var_mRNAs),]
dim(mRNA_expr_mat_red) # (6780, 81) # more close to the publication !!!

# mRNA RMA
mRNA_rma_expr_mat = readRDS(file = 'matrices/mRNA_rma_expr_mat.rds')
dim(mRNA_rma_expr_mat) # (7297, 81) #

# get std of all mRNAs
mRNA_stds = apply(mRNA_rma_expr_mat, 1, sd)
# filter to most variable ones
var_mRNAs = mRNA_stds[mRNA_stds > 0.8] # hard-coded (the higher, the more strict)
length(var_mRNAs) # 3089

mRNA_rma_expr_mat_red = mRNA_rma_expr_mat[names(var_mRNAs),]
dim(mRNA_rma_expr_mat_red) # (3089, 83) # HALF!!!!!!!

# save both reduced matrices
saveRDS(object = mRNA_expr_mat_red, file = 'matrices/mRNA_expr_mat_red.rds')
saveRDS(object = mRNA_rma_expr_mat_red, file = 'matrices/mRNA_rma_expr_mat_red.rds')

# miRNA
miRNA_expr_mat = readRDS(file = 'matrices/miRNA_expr_mat.rds')
dim(miRNA_expr_mat) # (874, 83) #

# get std of all miRNAs
miRNA_stds = apply(miRNA_expr_mat, 1, sd)
# filter to most variable ones
var_miRNAs = miRNA_stds[miRNA_stds > 0.5] # hard-coded
length(var_miRNAs) # 446

miRNA_expr_mat_red = miRNA_expr_mat[names(var_miRNAs),]
dim(miRNA_expr_mat_red) # (446, 83) #

# save reduced matrix
saveRDS(object = miRNA_expr_mat_red, file = 'matrices/miRNA_expr_mat_red.rds')

# What if we centered the data first? (NO NEED TO)
miRNA_expr_mat_centered = scale(miRNA_expr_mat, center = T, scale = F)
dim(miRNA_expr_mat_centered)

# get std of all miRNAs
miRNA_stds = apply(miRNA_expr_mat_centered, 1, sd)
# filter to most variable ones
var_miRNAs = miRNA_stds[miRNA_stds > 0.5]
length(var_miRNAs) # 450, a bit more, not so much different #

# What if we scaled the data first? (NO NEED TO)
miRNA_expr_mat_scaled = scale(miRNA_expr_mat, center = T, scale = T)
dim(miRNA_expr_mat_scaled)

# get std of all miRNAs
miRNA_stds = apply(miRNA_expr_mat_scaled, 1, sd)
# filter to most variable ones
var_miRNAs = miRNA_stds[miRNA_stds > 0.5]
length(var_miRNAs) # 29, A LOT LESS!!! #


