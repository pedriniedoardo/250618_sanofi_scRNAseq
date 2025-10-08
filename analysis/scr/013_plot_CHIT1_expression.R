# AIM ---------------------------------------------------------------------
# explore the expression of CHIT1 in the WM + CX dataset

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(lemon)
library(finalfit)
library(cowplot)
library(patchwork)
library(Nebulosa)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../../out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")
DimPlot(data.combined,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(data.combined,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(data.combined,label = T,raster = T,group.by = "expertAnno.l1",split.by = "origin")
DimPlot(data.combined,label = T,raster = T,group.by = "origin")
DimPlot(data.combined,label = T,raster = T,group.by = "pathology_class")

# str_subset(rownames(data.combined),pattern = "HIF")
# define the gene of interest GOI
# GOI <- c("Irf7","Ddx58")
GOI <- c("CHIT1","TSPO")

table(data.combined@meta.data$expertAnno.l1)

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
df_tot_avg <- df_tot %>% group_by(expertAnno.l1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

dim(df_tot)


head(data.combined@meta.data)

# plot the average expression per sample use the variable cell tyep per donor as grouping
# data.combined$group <- paste0(data.combined$orig.ident,".",data.combined$cell_type2)
data.combined$group <- paste0(data.combined$pathology_class,"-",data.combined$expertAnno.l1,"-",data.combined$orig.ident)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
Idents(data.combined) <- "group"
DefaultAssay(data.combined) <- "RNA"

average_GOI <- AverageExpression(data.combined,features = GOI,group.by = c("group"))

# plot general UMAP -------------------------------------------------------
# build the plot using both info
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=expertAnno.l1),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = expertAnno.l1),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()
# facet_wrap(~infection)
# ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident.pdf",width = 7,height = 5)

# no lab
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=expertAnno.l1),size=0.3) +
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
  facet_wrap(gene~pathology_class) +
  theme_cowplot() +
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_DPP3_count.pdf",width = 13,height = 12)

# do the same using Seurat
FeaturePlot(data.combined,features = GOI,split.by = "pathology_class",raster = T,order = T,ncol = 3)
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
  facet_wrap(gene~pathology_class) +
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
  facet_wrap(gene~pathology_class) +
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
  facet_wrap(gene~pathology_class)+
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
  ggplot(aes(x=pathology_class,y=exp_fix)) + 
  geom_violin(scale = "width")+
  geom_point(position=position_jitter(width = 0.2),alpha=0.01) +
  facet_wrap(~expertAnno.l1) +
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
  mutate(pathology_class = str_extract(group,pattern = c("CX_Ctrl|CX_Demye|CX_Mye|WM_CA|WM_CI|WM_Core|WM_Ctrl|WM_NAWM"))) |> 
  mutate(donor = str_extract(group,pattern = c("s\\d+"))) |> 
  mutate(expertAnno.l1 = str_extract(group,pattern = c("AST|EPENDYMA|EXC|NEU|IMM|INH|NEU|LYM|OLIGO|OPC|VAS")))

# plot the average expresison by cell annotation
df_avg |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=expertAnno.l1,y=avg_exp))+
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
  ggplot(aes(x=pathology_class,y=avg_exp))+
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

# rank the expertAnno.l1 category by the average expression per sample (regardelss of the gene)
rank_cell_id <- df_avg %>%
  group_by(expertAnno.l1) %>%
  summarise(avg = mean(avg_exp)) %>%
  arrange(desc(avg)) %>%
  pull(expertAnno.l1)

# generate the average estimates per cell annotation per gene
summary_exp_cell_id <- df_avg %>%
  group_by(expertAnno.l1, gene) %>%
  summarise(avg = mean(avg_exp)) %>%
  mutate(expertAnno.l1 = factor(expertAnno.l1,levels = rank_cell_id))

# plot splitting by treat full
df_avg |>
  mutate(expertAnno.l1 = factor(expertAnno.l1,levels = rank_cell_id)) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=pathology_class,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_grid(gene~expertAnno.l1,scales = "free") +
  geom_hline(data = summary_exp_cell_id,aes(yintercept = avg),col="red",linetype = "dashed")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/06_dotplot_annotationConfident_DPP3_expressionAvg_treatFull.pdf",width = 9,height = 9)

# plot splitting by treat
df_avg |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=pathology_class,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(gene~expertAnno.l1,scales = "free")
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
  facet_wrap(~expertAnno.l1,
             scales="free") +
  theme(strip.background = element_blank())

df_cor_stats <- df_avg %>%
  split(f = .$expertAnno.l1) %>%
  lapply(function(x){
    test <- x %>%
      pivot_wider(names_from = gene,values_from = avg_exp)
    
    cor.test(test$CHIT1,test$TSPO) %>%
      broom::tidy()
  }) %>%
  bind_rows(.id = "expertAnno.l1")
