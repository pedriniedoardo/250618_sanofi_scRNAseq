# AIM ---------------------------------------------------------------------
# put together the dataset from shaflick CSF, shafflick blood and full brain dataset WM+CX to compare CHIT and TSPO expression

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ggrepel)
library(cowplot)
library(DESeq2)
library(RNAseqQC)
library(limma)
library(ashr)
library(magick)
library(ComplexHeatmap)
library(pals)
library(scales)
library(patchwork)
library(SeuratWrappers)

# produce the table for the different runs --------------------------------
# read in the shafflick Blood dataset
scobj_BLOOD <- readRDS("../out/object/013_schafflick_blood_full.rds")
DimPlot(scobj_BLOOD,raster = T,group.by = "cell_id.l1",label = T)

# read in the CSF dataset
scobj_CSF <- readRDS("../out/object/013_schafflick_csf_full.rds")
DimPlot(scobj_CSF,raster = T,group.by = "cell_id.l1",label = T)

# read in the full brain dataset
scobj_BRAIN <- readRDS("/beegfs/scratch/ric.cosr/pedrini.edoardo/project_edoardo/220501_scRNAseq_MSbrain_Absinta/out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")
DimPlot(scobj_BRAIN,raster = T,group.by = "expertAnno.l1",label = T)

# simple aggregation of the data by cell type
cts_blood <- AggregateExpression(object = scobj_BLOOD,
                                 # group.by = c("cell_id.l1","sample_id"),
                                 group.by = c("cell_id.l1"),
                                 assays = 'RNA',
                                 slot = "counts",
                                 return.seurat = FALSE)

counts_blood <- cts_blood$RNA

# fix the colnames
colnames(counts_blood) <- paste0("blood_",str_replace_all(colnames(counts_blood),pattern = "\\s","."))
dim(counts_blood)

# simple aggregation of the data by cell type
cts_csf <- AggregateExpression(object = scobj_CSF,
                               # group.by = c("cell_id.l1","sample_id"),
                               group.by = c("cell_id.l1"),
                               assays = 'RNA',
                               slot = "counts",
                               return.seurat = FALSE)

counts_csf <- cts_csf$RNA

# fix the colnames
colnames(counts_csf) <- paste0("csf_",str_replace_all(colnames(counts_csf),pattern = "\\s","."))
dim(counts_csf)

# simple aggregation of the data by cell type
cts_brain <- AggregateExpression(object = scobj_BRAIN,
                               # group.by = c("cell_id.l1","sample_id"),
                               group.by = c("expertAnno.l1"),
                               assays = 'RNA',
                               slot = "counts",
                               return.seurat = FALSE)

counts_brain <- cts_brain$RNA

# fix the colnames
colnames(counts_brain) <- paste0("brain_",str_replace_all(colnames(counts_brain),pattern = "\\s","."))
dim(counts_brain)

# wrangling ---------------------------------------------------------------
# join the three tables make an inner join to so drop the non common genes across the two datasets
counts_full <- purrr::reduce(list(counts_blood %>%
                                    data.frame() %>%
                                    rownames_to_column("gene"),
                                  counts_csf %>%
                                    data.frame() %>%
                                    rownames_to_column("gene"),
                                  counts_brain %>%
                                    data.frame() %>%
                                    rownames_to_column("gene")),
                             inner_join,by = "gene") %>%
  column_to_rownames("gene")

# compare the dimensions
dim(counts_full)
dim(counts_blood)
dim(counts_csf)
dim(counts_brain)

# 2. generate the full metadata from the aggregated sc dataset
meta_full <- data.frame(sample = colnames(counts_full) %>% unname()) %>%
  separate(sample,into = c("dataset","cell_id"),sep="_",remove = F)

# save the metadata and the full matrix
saveRDS(counts_full,"../out/object/015_count_full.rds")
saveRDS(meta_full,"../out/object/015_meta_full.rds")

# buld the DESeq2 object --------------------------------------------------
# build the model, I am not really using it.
batch <- meta_full$dataset
design_full <- model.matrix(~ batch)
colnames(design_full)[1] <- c("intercept")

# save the disign
saveRDS(design_full,"../out/object/015_design_full.rds")

# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_full,
                              colData = meta_full,
                              design = design_full)

# filter
# filter out lowly expressed genes
keep_features <- edgeR::filterByExpr(counts(dds), group = meta_full$batch)
dds_filter <- dds[keep_features,]

# check the dimenstion before and after filtering
dim(dds)
dim(dds_filter)

# scale the data using vst. blind = T
vst_unfiltered <- vst(dds, blind = T)
vst_filter <- vst(dds_filter, blind = T)

# save the objecs
saveRDS(dds_filter,"../out/object/015_dds_filter.rds")
saveRDS(dds,"../out/object/015_dds.rds")

saveRDS(vst_filter,"../out/object/015_vst_filter.rds")
saveRDS(vst_unfiltered,"../out/object/015_vst_unfiltered.rds")

# explore the expression of GOIs ------------------------------------------
# define the GOIs
GOI <- c("TSPO", "CHIT1")

# perform the normalization
dds_full <- DESeq(dds)

# test if the genes are present
# notice If I use the filtered dataset CHIT1 get lost
counts(dds_full,normalized=T)%>%
  data.frame()%>%
  mutate(gene = rownames(.)) %>%
  dplyr::select(gene) %>%
  filter(gene %in% GOI)

# extract the tatal counts per sample
MR <- counts(dds_full,normalized=F)%>%
  data.frame()%>%
  gather(key = sample,value = exp)%>%
  # mutate(sample = str_replace(sample,pattern = "\\.",replacement = "-")) %>%
  group_by(sample)%>%
  summarise(MR = sum(exp)/10^6)

# make the table a long format to plot the expression per gene
count_norm <- counts(dds_full,normalized=T) %>%
  data.frame() %>%
  mutate(gene = rownames(.)) %>%
  dplyr::filter(gene%in%GOI) %>%
  gather(key = sample,value = count_norm,-gene)

count_raw <- counts(dds_full,normalized=F) %>%
  data.frame() %>%
  mutate(gene = rownames(.)) %>%
  dplyr::filter(gene%in%GOI) %>%
  gather(key = sample,value = count_raw,-gene)

df_counts_long <- left_join(count_norm,count_raw,by = c("gene","sample")) %>%
  # add the LUT
  # mutate(sample = str_replace(sample,pattern = "\\.",replacement = "-")) %>%
  left_join(meta_full,by = "sample") %>%
  # add the milion reads per sample
  left_join(MR,by = "sample") %>%
  # fix the sample name
  # left_join(lut,by = c("ens_id"="ensembl"))%>%
  mutate(count_norm_adj = count_norm + 0.5) %>%
  mutate(CPM = count_raw/MR) %>%
  # scale expression gene-wise
  group_by(gene) %>%
  mutate(scaled_exp = scale(count_norm_adj)[,1]) %>%
  ungroup()

# save the table for quicker exploration
write_tsv(df_counts_long,"../out/table/015_df_counts_long.tsv")

# confirm the scaling direction
df_counts_long %>%
  group_by(gene) %>%
  summarise(avg_gene = mean(scaled_exp),sd = sd(scaled_exp))

# plot the data as a dotplot
df_counts_long %>%
  ggplot(aes(x=sample,y = count_norm_adj)) +
  geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~gene,scales = "free") +
  scale_y_continuous(trans = "log10") + 
  theme_bw() +
  theme(strip.background =element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1))
ggsave("../out/plot/015_scatter_expression_normCount.pdf",width = 7,height = 3)

# plot CPM
df_counts_long %>%
  ggplot(aes(x=sample,y = CPM)) +
  geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~gene,scales = "free") +
  scale_y_continuous(trans = "log10") + 
  theme_bw() +
  theme(strip.background =element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1))
ggsave("../out/plot/015_scatter_expression_CPM.pdf",width = 7,height = 3)

# plot the data as an heatmap
df_counts_long %>%
  # force the order
  # mutate(group_id = factor(group_id,levels = c("AST","LYM","VAS","IMM","EXC NEU","INH NEU","OLIGO","OPC","EPENDYMA"))) %>% 
  # mutate(cat = factor(cat,levels = c("ASTRO","B_CELLS","T_CELLS","ENDOTHELIAL","PERICYTE","IMMUNE","NEURONS","OLIGO","OPC"))) %>% 
  ggplot(aes(y = gene,x = sample)) +
  # geom_tile(aes(fill = avg_exp.scaled))+
  geom_tile(aes(fill = scaled_exp), width = 0.95, height = 0.95, size = 0.2, color = "white")+
  # scale_radius(range = c(0, 8)) +
  # facet_grid(~dataset,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle = 90),
        strip.text.x = element_text(angle = 0)) +
  scale_fill_gradientn(colours = viridis::turbo(20),
                       # limits = c(-1,3),
                       oob = scales::squish,
                       name = 'norm exp z-scaled\n within genes') +
  theme(plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"))
ggsave("../out/plot/015_heatmap_expression.pdf",width = 8,height = 3)

# plot with complexheatmap
mat <- assay(vst_unfiltered)[rownames(vst_unfiltered) %in% GOI, ]

# z scale the expresison
# t(scale(t(mat)))
# mat %>%
#   data.frame() %>%
#   rownames_to_column("gene") %>%
#   pivot_longer(names_to = "sample",values_to = "exp",-gene) %>%
#   group_by(gene) %>%
#   mutate(scaled_exp = scale(exp)[,1]) %>%
#   ungroup() %>%
#   select(-exp) %>%
#   pivot_wider(names_from = "sample",values_from = "scaled_exp")

# z scale the table of expression per gene
mat2 <- (mat - rowMeans(mat))/rowSds(mat)

# change the rownames to match the symbol
# reorder the genes
# match(GOI,rownames(mat2))

# build the annotation object 
sample_ordered <- meta_full %>%
  dplyr::slice(match(colnames(mat2),meta_full$sample)) %>%
  pull(dataset)

# build the column annotation
column_ha <- HeatmapAnnotation(treat = sample_ordered, 
                               col = list(treat = c("blood" = "red", "csf" = "yellow","brain"="black")))

# build the heatmap
ht2 <- Heatmap(mat2,
               name = "exp",
               # column_title = "AMD RPE",
               top_annotation = column_ha,
               col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
# cluster_rows = F,
# right_annotation = row_ha,
# row_split = rep(c(1,2,3,4),c(2,3,4,7)))

# draw(ht2,heatmap_legend_side = "left")
pdf("../out/plot/015_heatmap_expression_complexHeatmap.pdf",width = 8,height = 3)
draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 10, 40), "mm"))
dev.off()

# try remove batch effect\ ------------------------------------------------
# create an equivalent object to be filled with the removed batch data
vst_unfilter_batch <- vst_unfiltered

mat_test <- assay(vst_unfiltered)
mat_testBatch <- limma::removeBatchEffect(mat_test, vst_filter$dataset)
assay(vst_unfilter_batch) <- mat_testBatch

# save the object after batch correction
saveRDS(vst_unfilter_batch,"../out/object/015_vst_unfiltered_batch.rds")

# attempt PCA to see the effect fo the batch correction
plot_vst_batch <- plotPCA(vst_unfilter_batch,
                          intgroup = c("dataset","cell_id")) +
  theme_bw()

# extract the pca values
df_plot_PCA_batch <- plot_vst_batch$data %>%
  # harmonize the name of the cells
  mutate(cellid_harmonized = case_when(str_detect(cell_id,pattern = "CD8|CD4")~"T",
                                       str_detect(cell_id,pattern = "NEU")~"NEU",
                                       T ~ cell_id))

# fix the colors, set the costum colors palette for the cell ids
myColors <- alphabet(length(unique(df_plot_PCA_batch$cellid_harmonized)))
names(myColors) <- unique(df_plot_PCA_batch$cellid_harmonized)
show_col(myColors)

custom_scale_color <- scale_color_manual(values = myColors)

# all sample after batch correction
p1_batch <- df_plot_PCA_batch %>%
  # make CSF.MS a single category
  # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
  ggplot(aes(x=PC1,y=PC2,col=cellid_harmonized)) +
  geom_point(size =4,alpha=0.6) +
  # ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vst_batch$labels[1]) + xlab(plot_vst_batch$labels[2]) +
  custom_scale_color +
  ggtitle("All")

# make one coloring by batch
p2_batch <- df_plot_PCA_batch %>%
  # make CSF.MS a single category
  # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
  ggplot(aes(x=PC1,y=PC2,col=dataset)) +
  geom_point(size =4,alpha=0.6) +
  # ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vst_batch$labels[1]) + xlab(plot_vst_batch$labels[2]) +
  ggtitle("Orgin")

# plot the panel
p1_batch + p2_batch + plot_annotation("batch correction")
ggsave("../out/plot/015_PCA_batch.pdf",width = 14,height = 6)

# z scale the table of expression per gene
# plot with complexheatmap
mat3 <- assay(vst_unfilter_batch)[rownames(vst_unfilter_batch) %in% GOI, ]
mat3_batch <- (mat3 - rowMeans(mat3))/rowSds(mat3)

# # build the annotation object 
# sample_ordered <- meta_full %>%
#   dplyr::slice(match(colnames(mat2),meta_full$sample)) %>%
#   pull(dataset)

# # build the column annotation
# column_ha <- HeatmapAnnotation(treat = sample_ordered, 
#                                col = list(treat = c("blood" = "red", "csf" = "yellow","brain"="black")))

# build the heatmap
ht3 <- Heatmap(mat3_batch,
               name = "exp",
               # column_title = "AMD RPE",
               top_annotation = column_ha,
               col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
# cluster_rows = F,
# right_annotation = row_ha,
# row_split = rep(c(1,2,3,4),c(2,3,4,7)))

# draw(ht2,heatmap_legend_side = "left")
pdf("../out/plot/015_heatmap_expression_complexHeatmap_batch.pdf",width = 8,height = 3)
draw(ht3,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 10, 40), "mm"))
dev.off()

# plot the VST as scatter
vst_value <- mat3 %>%
  data.frame() %>%
  mutate(gene = rownames(.)) %>%
  # dplyr::filter(gene%in%GOI) %>%
  pivot_longer(names_to = "sample",values_to = "vst",-gene) %>%
  left_join(meta_full,by = "sample") %>%
  # harmonize the name of the cells
  mutate(cellid_harmonized = case_when(str_detect(cell_id,pattern = "CD8|CD4")~"T",
                                       str_detect(cell_id,pattern = "NEU")~"NEU",
                                       T ~ cell_id)) %>%
  # scale expression gene-wise
  group_by(gene) %>%
  mutate(scaled_vst = scale(vst)[,1]) %>%
  ungroup()

# save the table for quicker exploration
write_tsv(vst_value,"../out/table/015_df_vst_long.tsv")

# plot CPM
vst_value %>%
  ggplot(aes(x=sample,y = vst)) +
  geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~gene,scales = "free") +
  # scale_y_continuous(trans = "log10") + 
  theme_bw() +
  theme(strip.background =element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1))
ggsave("../out/plot/015_scatter_expression_batch.pdf",width = 7,height = 3)

# plot the data as an heatmap
vst_value %>%
  # force the order
  # mutate(group_id = factor(group_id,levels = c("AST","LYM","VAS","IMM","EXC NEU","INH NEU","OLIGO","OPC","EPENDYMA"))) %>% 
  # mutate(cat = factor(cat,levels = c("ASTRO","B_CELLS","T_CELLS","ENDOTHELIAL","PERICYTE","IMMUNE","NEURONS","OLIGO","OPC"))) %>% 
  ggplot(aes(y = gene,x = sample)) +
  # geom_tile(aes(fill = avg_exp.scaled))+
  geom_tile(aes(fill = scaled_vst), width = 0.95, height = 0.95, size = 0.2, color = "white")+
  # scale_radius(range = c(0, 8)) +
  # facet_grid(~dataset,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle = 90),
        strip.text.x = element_text(angle = 0)) +
  scale_fill_gradientn(colours = viridis::turbo(20),
                       limits = c(0,2.5),
                       oob = scales::squish,
                       name = 'norm exp z-scaled\n within genes') +
  theme(plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"))
ggsave("../out/plot/015_heatmap_expression_batch.pdf",width = 8,height = 3)

