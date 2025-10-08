# AIM ---------------------------------------------------------------------
# explore the expression of CHIT1 in the shafflick CSF dataset

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(lemon)
# library(finalfit)
library(cowplot)
library(patchwork)
library(Nebulosa)
# library(viridis)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../data/data/schafflick_csf_harmony_by_patho.rds")
DimPlot(data.combined,label = T,raster = T,group.by = "seurat_clusters")
ggsave("../out/plot/013_UMAP_shafflick_csf_cluster.pdf",width = 5,height = 4)

data.combined@meta.data

# load the label transfer annotation done over the classical PBMC dataset
meta_data.combined <- read_tsv("../data/data/schafflick_csf_metaAzimuth.tsv")

# wrangling ---------------------------------------------------------------
# pull the sample id from the barcode information and add it as metadata
data.combined$sample_id <- data.combined@meta.data %>%
  rownames_to_column("barcode") %>%
  mutate(barcode = str_remove(barcode,pattern = "csf_ms_|csf_ctrl_")) %>%
  separate(barcode,into = c("sample_id","barcode_id"),sep = "_") %>%
  pull(sample_id)

# Identify the most likely assignment for each seurat cluster. this comes from Azimuth on the sc PBMC dataset.
prop_table_seurat_clusters <- meta_data.combined %>%
  group_by(seurat_clusters,predicted.celltype.l1) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.celltype.l1,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.celltype.l1") %>% 
  as.matrix()

pdf("../out/plot/013_heatmap_shafflick_csf_celltype.l1.pdf",height = 3,width = 4)
Heatmap(prop_table_seurat_clusters,
        name = "prop", 
        column_title = "celltype.l1",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

# based on this classification, add the assigned cell classification to the main object.
data.combined$cell_id.l1 <- data.combined@meta.data %>%
  mutate(cell_id.l1 = case_when(seurat_clusters %in% c(10,5,9)~"MONO",
                                seurat_clusters %in% c(8,11)~"B",
                                seurat_clusters %in% c(0,2)~"CD4 T",
                                seurat_clusters %in% c(3)~"CD4 T clu3",
                                seurat_clusters %in% c(6,12)~"DC",
                                seurat_clusters %in% c(8,12)~"other",
                                seurat_clusters %in% c(7)~"NK",
                                seurat_clusters %in% c(1)~"CD8 T clu1",
                                seurat_clusters %in% c(4)~"CD8 T clu4")) %>%
  pull(cell_id.l1)

# confirm the annotation
DimPlot(data.combined,group.by = "cell_id.l1",raster=T,label=T)
ggsave("../out/plot/013_UMAP_shafflick_csf_celltype.l1.pdf",width = 6,height = 5)

# confirm the annotaiton using a dotplot and a panel of known markers
shortlist_features_list <- list(
  "CD4 T" = c("CD4", "IL7R", "CCR7", "CD3E", "TRAC"),
  "CD8 T" = c("CD8A", "CD8B", "GZMB", "GZMH", "PRF1"),
  "B" = c("MS4A1", "CD19", "CD79A", "CD79B", "BANK1"),
  "Mono" = c("CD14", "LYZ", "S100A8", "S100A9", "FCGR3A"),
  "NK" = c("GNLY", "NKG7", "KLRD1", "FCGR3A", "NCAM1"),
  "DC" = c("ITGAX", "LILRA4", "FCER1A", "CLEC9A", "HLA-DPB1")
)

# define the ident
Idents(data.combined) <- "cell_id.l1"

# make the plot
test_long <- DotPlot(data.combined,
                     features = unique(unlist(shortlist_features_list)),
                     dot.scale = 8,cluster.idents = T)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

# build the reference to allow for multiple time the same gene
ref_tab <- test_long$data

# option 1 costume plot
test <- lapply(shortlist_features_list, function(x){
  ref_tab %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "feature.groups")

test %>% 
  filter(id %in% c("B","CD4 T","CD4 T clu3","CD8 T","CD8 T clu1","CD8 T clu4","DC","MONO","NK","other")) %>%
  mutate(id = factor(id, levels = c("B","CD4 T","CD4 T clu3","CD8 T","CD8 T clu1","CD8 T clu4","DC","MONO","NK","other"))) %>% 
  # mutate(id = fct_rev(id)) %>% 
  ggplot(aes(x=features.plot,y=id))+
  geom_point(aes(col=avg.exp.scaled,size=pct.exp),alpha=0.7)+
  facet_wrap(~feature.groups,scales = "free_x",nrow = 1)+
  theme_bw()+
  theme(strip.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_color_gradient(low = "gray",high = "blue")+
  # scale_color_viridis_c(option = "turbo")+
  scale_size_continuous(range = c(0.1, 10))

ggsave("../out/plot/013_dotplot_shafflick_csf_celltype.l1.pdf",width = 10,height = 4)

# save the object with full annotation
saveRDS(data.combined,"../out/object/013_schafflick_csf_full.rds")

# explore gene expression -------------------------------------------------

# str_subset(rownames(data.combined),pattern = "HIF")
# define the gene of interest GOI
# GOI <- c("Irf7","Ddx58")
GOI <- c("CHIT1","TSPO")

# add a broader cell grouping
# data.combined@meta.data$cell_type2 <- data.combined@meta.data |> 
#   mutate(cell_type2 = case_when(# according to label transfer cluster 12 and 6 are Bipolar cells
#     cell_type %in% c("BP_0","BP_10","BP_11","BP_13","BP_16","BP_17","BP_5","BP_8","12","6")~"BP",
#     T~cell_type)) |> 
#   pull(cell_type2)

# generate the table for the plots ----------------------------------------
# get the metadata from the other object
meta <- data.combined@meta.data %>%
  rownames_to_column(var = "barcodes")

# extrac the expression value
df_exp <- FetchData(data.combined, vars = GOI,slot = "data") |> 
  rownames_to_column("barcodes") |> 
  pivot_longer(names_to = "gene",values_to = "exp",-barcodes) |> 
  # try to min/max normalize the count varaible per gene in order to rescale the difference in terms of expression
  group_by(gene) %>%
  # threshold of positiveness is based on the distriubtion of the expression of the signal in tihs case
  mutate(exp_min_max = ((exp - min(exp))/(max(exp)-min(exp))),
         exp_cat = case_when(exp > 0~"pos",
                             T~"neg")) %>%
  ungroup() %>%
  mutate(exp_fix = exp + rnorm(nrow(.))/100000)

# get the coordinates
UMAP1_df <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column(var = "barcodes")

# generate the dataset for mapping the data in the umamp
dim(UMAP1_df)
dim(df_exp)
dim(meta)

df_tot <- purrr::reduce(list(meta,UMAP1_df,df_exp),left_join, by="barcodes")
df_tot_avg <- df_tot %>% group_by(cell_id.l1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

dim(df_tot)

head(data.combined@meta.data)

# plot the average expression per sample use the variable cell tyep per donor as grouping
# data.combined$group <- paste0(data.combined$orig.ident,".",data.combined$cell_type2)
data.combined$group <- paste0(data.combined$pathology,"-",data.combined$cell_id.l1,"-",data.combined$sample_id)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
Idents(data.combined) <- "group"
DefaultAssay(data.combined) <- "RNA"

average_GOI <- AverageExpression(data.combined,features = GOI,group.by = c("group"))

# plot general UMAP -------------------------------------------------------
# build the plot using both info
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=cell_id.l1),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = cell_id.l1),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()
# facet_wrap(~infection)
# ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident.pdf",width = 7,height = 5)

# no lab
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=cell_id.l1),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = cell_type2),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()
# facet_wrap(~infection)
# ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident_noLab.pdf",width = 7,height = 5)

# expression distribution -------------------------------------------------
# crop the 0 expressing cells
df_exp %>%
  ggplot(aes(x=exp))+geom_histogram()+facet_grid(~gene)+theme_bw()+scale_x_log10()+geom_vline(xintercept = 1,col="red",linetype="dotted")

# keep the 0 expressing cells
df_exp %>%
  ggplot(aes(x=exp))+geom_histogram()+facet_wrap(~gene)+theme_bw()+
  # scale_x_log10()+
  geom_vline(xintercept = 2.5,col="red",linetype="dotted")

# library(scales)
# show_col(c("#4662D7FF","#FABA39FF","#7A0403FF"))

# plotting expression -----------------------------------------------------
# by counts
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(gene~pathology) +
  theme_cowplot() +
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_DPP3_count.pdf",width = 13,height = 12)

# do the same using Seurat
FeaturePlot(data.combined,features = GOI,split.by = "pathology",raster = T,order = T,ncol = 2)
# ggsave("../../out/image/06_UMAPSeurat_annotationConfident_DPP3_count.pdf",width = 25,height = 3)

df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(~gene) +
  theme_cowplot() +
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_DPP3_count2.pdf",width = 6,height = 5)

# do the same using Seurat
FeaturePlot(data.combined,features = GOI,raster = T,order = T)
# ggsave("../../out/image/06_UMAPSeurat_annotationConfident_DPP3_count.pdf",width = 6,height = 5)

# try nebulosa option
plot_density(data.combined, GOI,reduction = "umap")

df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp)) + geom_point(alpha = 0.5,size = 0.05) +
  facet_wrap(gene~pathology) +
  theme_cowplot() +
  # scale_color_gradient(low = "gray",high = "blue") +
  scale_color_viridis_c(option = "turbo") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_DPP3_count_alt.pdf",width = 13,height = 12)

# by min max normalized counts
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_min_max)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(gene~pathology) +
  theme_cowplot() +
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_DPP3_minmax.pdf",width = 13,height = 12)

# plot the category. being 0 or non zero per cell
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_cat)) + geom_point(alpha = 0.5,size = 0.05) +
  # facet_wrap(gene~NMDA_time,nrow = 2) +
  # facet_rep_wrap(gene~treat,repeat.tick.labels = "all",nrow=3)+
  facet_wrap(gene~pathology)+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_cowplot() +
  scale_color_manual(values = c("gray","blue")) +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_DPP3_proppos.pdf",width = 13,height = 12)

df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_cat)) + geom_point(alpha = 0.5,size = 0.05) +
  # facet_wrap(gene~NMDA_time,nrow = 2) +
  # facet_rep_wrap(gene~treat,repeat.tick.labels = "all",nrow=3)+
  facet_wrap(~gene)+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_cowplot() +
  scale_color_manual(values = c("gray","blue")) +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_DPP3_proppos2.pdf",width = 6,height = 5)

# violin plot for GOI expression use macro categories
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  # this is the processing shown in the violinplot function
  # mutate(exp_fix = exp + rnorm(nrow(.))/100000) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=pathology,y=exp_fix)) + 
  geom_violin(scale = "width")+
  geom_point(position=position_jitter(width = 0.2),alpha=0.01) +
  facet_wrap(~cell_id.l1) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1,angle = 90)) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# ggsave("../../out/image/00_violin_annotationConfident_DPP3.pdf",width = 15,height = 10)

# try to depict the average expression there is roughly one sample per condition
df_avg <- average_GOI$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI) |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  # filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(pathology = str_extract(group,pattern = c("CSF.CTRL|CSF.MS"))) |> 
  mutate(donor = str_extract(group,pattern = c("P\\w+|MS\\w+"))) |> 
  mutate(cell_id.l1 = str_extract(group,pattern = c("B|CD4.T|CD8.T|DC|MONO|NK|other")))

# plot the average expresison by cell annotation
df_avg |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=cell_id.l1,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~gene,scales = "free")+
  scale_y_continuous(trans = "log1p")

# do the same as above but split by condition
df_avg |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=pathology,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~gene,scales = "free") +
  scale_y_continuous(trans = "log1p")

# rank the cell_id.l1 category by the average expression per sample (regardelss of the gene)
rank_cell_id <- df_avg %>%
  group_by(cell_id.l1) %>%
  summarise(avg = mean(avg_exp)) %>%
  arrange(desc(avg)) %>%
  pull(cell_id.l1)

# generate the average estimates per cell annotation per gene
summary_exp_cell_id <- df_avg %>%
  group_by(cell_id.l1, gene) %>%
  summarise(avg = mean(avg_exp)) %>%
  mutate(cell_id.l1 = factor(cell_id.l1,levels = rank_cell_id))

# plot splitting by treat full
df_avg |>
  mutate(cell_id.l1 = factor(cell_id.l1,levels = rank_cell_id)) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=pathology,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1,height = 0),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_grid(gene~cell_id.l1,scales = "free") +
  geom_hline(data = summary_exp_cell_id,aes(yintercept = avg),col="red",linetype = "dashed")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/06_dotplot_annotationConfident_DPP3_expressionAvg_treatFull.pdf",width = 9,height = 9)

# plot splitting by treat
df_avg |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=pathology,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(gene~cell_id.l1,scales = "free")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/00_dotplot_annotationConfident_DPP3_expressionAvg_treat.pdf",width = 9,height = 9)

# -------------------------------------------------------------------------
# martina was also interested in explorign the relative correlation of expression for the two genes.
# try to check a general correlation across all the averaged samples
df_avg_wide <- df_avg %>%
  pivot_wider(names_from = gene,values_from = avg_exp)

df_avg_wide %>%
  ggplot(aes(x=CHIT1,y=TSPO)) +
  geom_smooth(method = "lm") +
  geom_point() +
  scale_x_continuous(trans = "log1p") +
  scale_y_continuous(trans = "log1p") +
  theme_bw() +
  facet_wrap(~cell_id.l1,
             scales="free") +
  theme(strip.background = element_blank())

df_cor_stats <- df_avg %>%
  split(f = .$cell_id.l1) %>%
  lapply(function(x){
    test <- x %>%
      pivot_wider(names_from = gene,values_from = avg_exp)
    
    cor.test(test$CHIT1,test$TSPO) %>%
      broom::tidy()
  }) %>%
  bind_rows(.id = "cell_id.l1")
