# AIM ---------------------------------------------------------------------
# check the expression of CD40 and CD40LG in the different dataset available


# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# read in the data --------------------------------------------------------
# brain dataset
brain <- readRDS("/beegfs/scratch/ric.cosr/pedrini.edoardo/project_edoardo/220501_scRNAseq_MSbrain_Absinta/out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")

DimPlot(brain,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(brain,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(brain,label = T,raster = T,group.by = "expertAnno.l1",split.by = "origin")
DimPlot(brain,label = T,raster = T,group.by = "origin")
DimPlot(brain,label = T,raster = T,group.by = "pathology_class")

# read in the shafflick CSF dataset and its new annotation (from label transfer)
shafflick_CSF <- readRDS("../data/data/schafflick_csf_harmony_by_patho.rds")
DimPlot(shafflick_CSF,reduction = "umap",label = T)

# pull the meta with the cell type annotation
meta_shafflick_CSF <- read_tsv("../data/out/table/schafflick_csf_metaAzimuth.tsv") %>%
  mutate(robust_score_celltype.l1 = case_when(predicted.celltype.l1.score>0.70&mapping.score>0.70~predicted.celltype.l1,
                                              T~"uncertain"),
         robust_score_celltype.l2 = case_when(predicted.celltype.l2.score>0.70&mapping.score>0.70~predicted.celltype.l2,
                                              T~"uncertain"),
         robust_score_celltype.l3 = case_when(predicted.celltype.l3.score>0.70&mapping.score>0.70~predicted.celltype.l3,
                                              T~"uncertain"))

# confirm the dimensions are congruent
dim(shafflick_CSF@meta.data)
dim(meta_shafflick_CSF) 

# update the object
shafflick_CSF@meta.data <- 
  meta_shafflick_CSF %>%
  separate(barcodes,into = c("id1","id2","donor","barcode_id"),remove=F,sep = "_") %>%
  column_to_rownames("barcodes")

# check the plot with the new annotation
DimPlot(shafflick_CSF,reduction = "umap",label = T,group.by = "robust_score_celltype.l1")
DimPlot(shafflick_CSF,reduction = "umap",label = T,group.by = "seurat_clusters")


# read in the shafflick CSF dataset and its new annotation (from label transfer)
shafflick_BLOOD <- readRDS("../data/data/schafflick_blood_harmony_by_origident.rds")
DimPlot(shafflick_BLOOD,reduction = "umap",label = T)

# pull the meta with the cell type annotation
meta_shafflick_BLOOD <- read_tsv("../data/out/table/schafflick_blood_metaAzimuth.tsv") %>%
  mutate(robust_score_celltype.l1 = case_when(predicted.celltype.l1.score>0.70&mapping.score>0.70~predicted.celltype.l1,
                                              T~"uncertain"),
         robust_score_celltype.l2 = case_when(predicted.celltype.l2.score>0.70&mapping.score>0.70~predicted.celltype.l2,
                                              T~"uncertain"),
         robust_score_celltype.l3 = case_when(predicted.celltype.l3.score>0.70&mapping.score>0.70~predicted.celltype.l3,
                                              T~"uncertain"))

# confirm the dimensions are congruent
dim(shafflick_BLOOD@meta.data)
dim(meta_shafflick_BLOOD) 

# update the object
shafflick_BLOOD@meta.data <- 
  meta_shafflick_BLOOD %>%
  separate(barcodes,into = c("id1","id2","donor","barcode_id"),remove=F,sep = "_") %>%
  column_to_rownames("barcodes")

# check the plot with the new annotation
DimPlot(shafflick_BLOOD,reduction = "umap",label = T,group.by = "robust_score_celltype.l1")
DimPlot(shafflick_BLOOD,reduction = "umap",label = T,group.by = "seurat_clusters")

# define the genes of interest
GOI <- c("CD40","CD40LG")

# wraingling --------------------------------------------------------------
# brain
brain$group <- paste0(brain$pathology_class,"-",brain$expertAnno.l1,"-",brain$orig.ident)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
Idents(brain) <- "group"
DefaultAssay(brain) <- "RNA"

average_GOI_brain <- AverageExpression(brain,features = GOI,group.by = c("group"),layer = "data")


# try to depict the average expression there is roughly one sample per condition
df_avg_brain <- average_GOI_brain$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI) |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  # filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(pathology_class = str_extract(group,pattern = c("CX.Ctrl|CX.Demye|CX.Mye|WM.CA|WM.CI|WM.Core|WM.Ctrl|WM.NAWM"))) |> 
  mutate(donor = str_extract(group,pattern = c("s\\d+"))) |> 
  mutate(expertAnno.l1 = str_extract(group,pattern = c("AST|EPENDYMA|EXC|NEU|IMM|INH|NEU|LYM|OLIGO|OPC|VAS")))

# csf 
# shafflick_CSF@meta.data$donor %>% table()
shafflick_CSF$group <- paste0(shafflick_CSF$pathology,"-",shafflick_CSF$robust_score_celltype.l1,"-",shafflick_CSF$donor)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
Idents(shafflick_CSF) <- "group"
DefaultAssay(shafflick_CSF) <- "RNA"

average_GOI_CSF <- AverageExpression(shafflick_CSF,features = GOI,group.by = c("group"),layer = "data")

# try to depict the average expression there is roughly one sample per condition
df_avg_CSF <- average_GOI_CSF$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI) |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  # filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(pathology_class = str_extract(group,pattern = c("CSF.CTRL|CSF.MS"))) |> 
  mutate(donor = str_extract(group,pattern = c("MS19270|MS49131|MS58637|MS60249|MS71658|MS74594|PST45044|PST95809|PTC32190|PTC41540|PTC85037"))) |>
  mutate(expertAnno.l1 = str_extract(group,pattern = c("B|CD4.T|CD8.T|DC|Mono|NK|other|other.T|uncertain")))

df_avg_CSF %>%
  filter(is.na(expertAnno.l1)|is.na(donor)|is.na(pathology_class))

# blood
# shafflick_BLOOD@meta.data$donor %>% table()
shafflick_BLOOD$group <- paste0(shafflick_BLOOD$pathology,"-",shafflick_BLOOD$robust_score_celltype.l1,"-",shafflick_BLOOD$donor)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
Idents(shafflick_BLOOD) <- "group"
DefaultAssay(shafflick_BLOOD) <- "RNA"

average_GOI_BLOOD <- AverageExpression(shafflick_BLOOD,features = GOI,group.by = c("group"),layer = "data")

# try to depict the average expression there is roughly one sample per condition
df_avg_BLOOD <- average_GOI_BLOOD$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI) |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  # filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(pathology_class = str_extract(group,pattern = c("blood.ctrl|blood.MS"))) |> 
  mutate(donor = str_extract(group,pattern = c("MS19270|MS49131|MS60249|MS71658|MS74594|PST83775|PST95809|PTC32190|PTC41540|PTC85037"))) |>
  mutate(expertAnno.l1 = str_extract(group,pattern = c("B|CD4.T|CD8.T|DC|Mono|NK|other|other.T|uncertain")))

df_avg_BLOOD %>%
  filter(is.na(expertAnno.l1)|is.na(donor)|is.na(pathology_class))

# plotting ----------------------------------------------------------------
# brain
# plot the average expresison by cell annotation
df_avg_brain |>
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

# add splitting by condition
df_avg_brain |>
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
  facet_grid(expertAnno.l1~gene,scales = "free")+
  scale_y_continuous(trans = "log1p")

# CSF
# plot the average expresison by cell annotation
df_avg_CSF |>
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

# add splitting by condition
df_avg_CSF |>
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
  facet_grid(expertAnno.l1~gene,scales = "free")+
  scale_y_continuous(trans = "log1p")

# BLOOD
# plot the average expresison by cell annotation
df_avg_BLOOD |>
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

# add splitting by condition
df_avg_BLOOD |>
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
  facet_grid(expertAnno.l1~gene,scales = "free")+
  scale_y_continuous(trans = "log1p")
