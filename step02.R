#################################
######-Ambient RNA removal-######
#################################
library(SoupX)
library(Seurat)
library(cowplot)

tod <- Read10X("./raw_matrix",gene.column=1)
toc <- Read10X("./filter_matrix",gene.column=1)
tod <- tod[rownames(toc),]
all <- toc
all <- CreateSeuratObject(all)
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(all)
all <- ScaleData(all, features = all.genes)
all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = F)
all <- FindNeighbors(all, dims = 1:20)
all <- FindClusters(all, resolution = 0.5)
all <- RunUMAP(all, dims = 1:20)
matx <- all@meta.data
sc = SoupChannel(tod, toc)
sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
#sc = setContaminationFraction(sc, 0.2)
sc = autoEstCont(sc)
out = adjustCounts(sc)
DropletUtils:::write10xCounts("./soupX_matrix", out,version="3")

#################################
########-Doublet removal-########
#################################
library(DoubletFinder)
library(Seurat)
library(cowplot)

input=Read10X(data.dir = "./soupX_matrix",gene.column = 1)
input <- CreateSeuratObject(counts = input, project = "input", min.cells = 3, min.features = 200)
input[["percent.mt"]]<-PercentageFeatureSet(input,pattern = "^MT-")
input <- subset(input, subset =  nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)
input <- NormalizeData(input, normalization.method = "LogNormalize", scale.factor = 10000)
input <- FindVariableFeatures(input, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(input)
input <- ScaleData(input, features = all.genes)
input <- RunPCA(input, features = VariableFeatures(object = input))
DimPlot(object = input, group.by="orig.ident",reduction = "pca")
input <- JackStraw(input, num.replicate = 100)
input <- ScoreJackStraw(input, dims = 1:20)
n_PC <- max(which(input@reductions[["pca"]]@jackstraw@overall.p.values[,2]<0.01))
input <- RunUMAP(input, reduction = "pca", dims = 1:n_PC)
input <- RunTSNE(input, reduction = "pca", dims = 1:n_PC)
input_double <- paramSweep_v3(input, PCs = 1:20, sct = FALSE)
input_double <- summarizeSweep(input_double, GT = FALSE)
input_double <- find.pK(input_double)
pK_value <- as.numeric(as.character(input_double$pK[input_double$BCmetric == max(input_double$BCmetric)]))
annotations <- input@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.075*length(input@meta.data$orig.ident))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
input <- doubletFinder_v3(input, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
df.classifications <- paste("DF.classifications", pN_value, pK_value, nExp_poi, sep="_")
input <- subset(input, subset = df.classifications == "Singlet")

##############################################################
########-Cell clustering and cell type identification-########
##############################################################
library(Seurat)
library(cowplot)

input <- FindVariableFeatures(input, selection.method = "vst", nfeatures = 2000)
input <- ScaleData(input, features = VariableFeatures(object = input))
input <- RunPCA(input, verbose = FALSE)
n_PC <- max(which(input@reductions[["pca"]]@jackstraw@overall.p.values[,2]<0.01))
input <- RunUMAP(input,reduction = "pca", dims = 1:n_PC)
input <- RunTSNE(input,reduction = "pca",dims = 1:n_PC)
input <- FindNeighbors(input, dims = 1:n_PC)
input <- FindClusters(input, resolution = 1.0)
DimPlot(input, reduction = "umap",pt.size = 1,label = T)
DimPlot(input, reduction = "tsne",pt.size = 1,label = T)
input_makers<- FindAllMarkers(input, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FeaturePlot(input, features = c("XXXX"), min.cutoff = "q9",pt.size = 1,reduction = "umap",raster=FALSE, order = T)

##################################
########-Merge all tissue-########
##################################
library(Seurat)
library(cowplot)

Sample<-readRDS("./Sample.rds")
Sample@meta.data$tissue <- "Sample"
tissue_all<-merge(x=Sample_1,y=Sample_list,project = "Donkey")
tissue_all <- NormalizeData(object = tissue_all, normalization.method = "LogNormalize")
tissue_all<- FindVariableFeatures(object = tissue_all)
tissue_all <- ScaleData(tissue_all)
tissue_all <- RunPCA(tissue_all)
tissue_all <- RunUMAP(tissue_all, reduction = "pca", dims = 1:20)
tissue_all <- RunTSNE(tissue_all, reduction = "pca", dims = 1:20)
tissue_all <- FindNeighbors(tissue_all, reduction = "pca", dims = 1:20)
tissue_all <- FindClusters(tissue_all, resolution = 0.5)
DimPlot(tissue_all, reduction = "tsne", group.by = "tissue", label = T,pt.size = 0.1,raster=FALSE)

#########################
########-hdWGCNA-########
#########################
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(cowplot)
library(MetBrewer)

Sample<-readRDS("./Sample.rds")
WGCNA_Sample <- SetupForWGCNA(
   Sample,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = "WGCNA")

WGCNA_Sample <- MetacellsByGroups(
    seurat_obj = WGCNA_Sample,
    group.by = c("system"), 
    reduction = "pca", 
    k = 25, 
    max_shared = 10, 
    ident.group = "system")

WGCNA_Sample <- NormalizeMetacells(WGCNA_Sample)
WGCNA_Sample <- SetDatExpr(
  WGCNA_Sample,
  group_name = c("Circulatory", "Digestive", "Endocrine", 
  "Immune", "Integumentary", "Muscular", "Nervous", 
  "Female Reproductive", "Male Reproductive", "Respiratory", "Urinary"),
  group.by='system', 
  assay = 'RNA', 
  slot = 'data')

WGCNA_Sample <- TestSoftPowers(
    WGCNA_Sample,
    networkType = "signed")
	
plot_list <- PlotSoftPowers(WGCNA_Sample)
wrap_plots(plot_list, ncol=2)
WGCNA_Sample <- ConstructNetwork(
    WGCNA_Sample, soft_power = 7,
    setDatExpr = FALSE,
    tom_name = 'INH')

PlotDendrogram(WGCNA_Sample)
TOM <- GetTOM(WGCNA_Sample)
WGCNA_Sample <- ScaleData(WGCNA_Sample, features = VariableFeatures(WGCNA_Sample))
WGCNA_Sample <- ModuleEigengenes(WGCNA_Sample, group.by.vars = "system")
WGCNA_Sample <- ModuleConnectivity(WGCNA_Sample,
  group.by = 'system', group_name = c("Circulatory", "Digestive", "Endocrine", 
  "Immune", "Integumentary", "Muscular", "Nervous", 
  "Female Reproductive", "Male Reproductive", "Respiratory", "Urinary"))

modules <- GetModules(WGCNA_Sample)
hub_df <- GetHubGenes(seurat_obj = WGCNA_Sample, n_hubs = 10)
WGCNA_Sample <- ModuleExprScore(WGCNA_Sample,n_genes = 25,method = "Seurat")
plot_list <- ModuleFeaturePlot(WGCNA_Sample,features = 'hMEs',order = TRUE)
wrap_plots(plot_list, ncol = 6)
ModuleNetworkPlot(seurat_obj = WGCNA_Sample, mods = "M1")
WGCNA_Sample <- RunModuleUMAP(WGCNA_Sample, 
    n_hubs = 10, 
    n_neighbors = 15, 
    min_dist = 0.1)
	
ModuleUMAPPlot(WGCNA_Sample,
    edge.alpha = 0.25,
    sample_edges = T,
    edge_prop = 0.1,
    label_hubs = 2,
    keep_grey_edges = F)

Idents(WGCNA_Sample) <- WGCNA_Sample$system
markers <- Seurat::FindAllMarkers(
  WGCNA_Sample,
  only.pos = TRUE,
  logfc.threshold=1)

overlap_df <- OverlapModulesDEGs(
  WGCNA_Sample,
  deg_df = markers,
  fc_cutoff = 1)
  
OverlapDotPlot(overlap_df)
  
##########################
########-pySCENIC-########
##########################
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)

Sample<-readRDS("./Sample.rds")
write.csv(t(as.matrix(Sample@assays$RNA@counts)), file = "./SCENIC_Sample.csv")

#import loompy as lp;
#import numpy as np;
#import scanpy as sc;
#x=sc.read_csv("SCENIC_Sample.csv");
#row_attrs = {"Gene": np.array(x.var_names),};
#col_attrs = {"CellID": np.array(x.obs_names)};
#lp.create("SCENIC_Sample.loom",x.X.transpose(),row_attrs,col_attrs);

pyscenic grn \
--num_workers 40 \
--output Sample.tsv \
--method grnboost2 \
SCENIC_Sample.loom  \
./tfs.txt

pyscenic ctx \
Sample.tsv \
./genes_vs_motifs.rankings.feather \
--annotations_fname ./motifs.tbl \
--expression_mtx_fname SCENIC_Sample.loom  \
--mode "dask_multiprocessing" \
--output Sample.csv \
--num_workers 40  \
--mask_dropouts

pyscenic aucell \
SCENIC_Sample.loom \
Sample.csv \
--output SCENIC.loom \
--num_workers 40 

loom <- open_loom('SCENIC.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulonAUC <- get_regulons_AUC(loom,column.attr.name = 'RegulonsAUC')
cellInfo <-Sample@meta.data
cellTypes <-  as.data.frame(subset(cellInfo,select = 'celltype'))
selectedResolution <- "celltype"
cellAnnotation = cellTypes[colnames(regulonAUC), selectedResolution]
cellAnnotation = na.omit(cellAnnotation)
rss <- calcRSS(AUC = getAUC(regulonAUC), cellAnnotation = cellAnnotation)
rss=na.omit(rss)

###############################
########-Calculate CSI-########
###############################
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

loom <- open_loom('SCENIC.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulon <- regulons_incidMat@assays@data$AUC
regulon <- t(regulon)
module_TF <- CSI_matrix_cal(regulon = regulon,
                            CSI_threshold = 0.9,
                            module_k = 8,
                            legend_parm = "number",
                            rect_color="red")
draw(module_TF)

######################################################
########-Species homeotic gene transformation-########
######################################################
library(Seurat)
library(biomaRt)

count_raw <- as.matrix(Species@assays$RNA@counts)
usGenes  <- rownames(count_raw)
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="https://dec2021.archive.ensembl.org/")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
Species = useMart("ensembl", dataset = "Species_gene_ensembl")
data = getLDS(attributes = c("Species_symbol"), filters = "Species_symbol", values = usGenes , mart = Species, attributesL = c("hgnc_symbol","ensembl_gene_id"), martL = human, uniqueRows=T)
data <- data[!duplicated(data[,1]),]
data <- data[!duplicated(data[,2]),]

count_raw <- as.data.frame(count_raw)  
count_raw <- subset(count_raw, rownames(count_raw) %in% data$Gene.name)
count_raw$Gene.name <- rownames(count_raw)
count_final <- merge(data, count_raw, by="Gene.name", all=T) 
rownames(count_final) <- count_final$Gene.name.1
count_final[,c("Gene.name","Gene.name.1")] <- NULL 
length(grep(pattern = "-1",colnames(count_final)))
data <- na.omit(count_final)
nonzero <- data > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- data[keep_genes, ]
dim(filtered_counts) 
cellInfo <- Species@meta.data
Species <- CreateSeuratObject(counts = filtered_counts,
                              meta.data =cellInfo,
                              project = "Species",
                              min.features = 200,
                              min.cells = 3)

##########################################################
########-Cross-species scRNA-seq data integration-########
##########################################################
library(Seurat)

human<-readRDS("./human.rds")
donkey<-readRDS("./donkey.rds")
mouse<-readRDS("./mouse.rds")
pig<-readRDS("./Pig.rds")
Merge <- merge(human,y=c(donkey,mouse,pig), project = "Merge")
species.list <- SplitObject(Merge, split.by = "species")
species.list <- lapply(X = species.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = species.list)
species.list <- lapply(X = species.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
species.anchors <- FindIntegrationAnchors(object.list = species.list, anchor.features = features, reduction = "rpca")
species.combined <- IntegrateData(anchorset = species.anchors)
DefaultAssay(species.combined) <- "integrated"
species.combined <- ScaleData(species.combined, verbose = FALSE)
species.combined <- RunPCA(species.combined, npcs = 30, verbose = FALSE)
species.combined <- RunUMAP(species.combined, reduction = "pca", dims = 1:30)
species.combined <- RunTSNE(species.combined, reduction = "pca", dims = 1:30,check_duplicates=F)
species.combined <- FindNeighbors(species.combined, reduction = "pca", dims = 1:30)
species.combined <- FindClusters(species.combined, resolution = 0.5)
DimPlot(species.combined, reduction = "tsne", group.by = "species", label = TRUE, repel = TRUE,raster=FALSE)
DimPlot(species.combined, reduction = "tsne", group.by = "tissue", label = TRUE, repel = TRUE,raster=FALSE)

##############################
########-MetaNeighbor-########
##############################
library(MetaNeighbor)
library(SingleCellExperiment)

Sample<-readRDS("./Sample.rds")
Sample.sce<-as.SingleCellExperiment(Sample)
aurocs <- MetaNeighborUS(var_genes = gene_list,
                         dat = All.sce,
                         study_id = Sample.sce$specie,
                         cell_type = Sample.sce$celltype,
                         fast_version=TRUE)
plotHeatmap(aurocs,cex=0.3)

####################################################
########-Cross-species correlation analysis-########
####################################################
library(ggcorrplot)
library(reshape2)

corr<-read.csv("./celltype_AverageExpression.csv")
cor_data<-cor(corr,method="spearman")
head(cor_data)
corr_data<-melt(cor_dat)

#################################################
########-Single-Cell-Metabolic-Landscape-########
#################################################
library(stringr)
library(reshape2)
library(scales)
library(scater)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(rjson)

Sample<-readRDS("./Sample.rds")
Sample$newcelltype <- paste0(Sample$tissue,"_", Sample$celltype)
all_cell_types <- as.vector(Sample$newcelltype)
cell_types <- unique(all_cell_types)
metabolism_activaty <- Pathway_act_Score(Sample,
                                         pathways=pathways,
                                         assay = "RNA",
                                         filterGene=T,
                                         Mean_cut = 0.001,
                                         percent_cut =0.1,
                                         all_cell_types = all_cell_types,
                                         cell_types = cell_types)

data <- metabolism_activaty[[1]]
all_NA <- rowAlls(is.na(data))
data <- data[!all_NA,]

KEGG_pathway <- rjson::fromJSON(file = "./br08901.json",simplify=F)
KEGG_pathway <- KEGG_pathway[[2]]
KEGG_pathway <- KEGG_pathway[[1]]
KEGG_pathway <- KEGG_pathway[[2]]

cate1 <- c()
for (i in 1:length(KEGG_pathway)) {
  length_names <- length(KEGG_pathway[[i]][[2]])
  cate1 <- append(cate1, rep(KEGG_pathway[[i]][[1]], length_names))
   }
 
cate2 <- c()
for (i in 1:length(KEGG_pathway)) {
  cate2 <- append(cate2, KEGG_pathway[[i]][[2]])
  }

cate3 <- c()
for (i in 1:length(cate2)) {
   cate3 <- append(cate3, as.character(cate2[[i]]))
  }

vec = cate3
vec <- str_replace(vec, "[0-9]+", "")
vec <- str_replace(vec, "  ", "")
KEGG_pathways_Metabolism_anno <- data.frame(cat1  = cate1,
                                        cat2  = cate3,
                                        vec  = vec)
 
rownames(KEGG_pathways_Metabolism_anno) <- KEGG_pathways_Metabolism_anno$vec
Anno_pathway <- KEGG_pathways_Metabolism_anno[rownames(data),]
data$pathway_anno <- Anno_pathway$directory2
data$pathway <- rownames(data)
melt_data = melt(data,
                id.vars = c("pathway","pathway_anno"),
                measure.vars = 1:(length(data)-2),
                variable.name = c('celltype'),
                value.name = 'score')

#######################
########-scFEA-########
#######################
###The single-cell data from donkeys should initially undergo conversion to their human homologous genes.
library(Seurat)
Sample<-readRDS("./Sample.rds")
Sample <- as.matrix(Sample@assays$RNA@counts)
write.csv(Sample, file = "./scFEA_Sample.csv")

python src/scFEA.py --data_dir data --input_dir input/ \
                    --test_file scFEA_Sample.csv \
                    --moduleGene_file module_gene_m168.csv \
                    --stoichiometry_matrix cmMat_c70_m168.csv \
                    --res_dir output \
                    --sc_imputation True

#####################
########-PCA-########
#####################
library(ggplot2)
library(tidyr)
library(dplyr)

data=read.csv("./metabolome.csv")
data <- t(data)  
group=c(rep('species1',5),
        rep('species2',5))

pca_data <- prcomp(data)
screeplot(pca_data, type = "lines")
summary(pca_data)
rownames(pca_data$x)
x_per <- round(summary(pca_data)$importance[2, 1]*100, 1)
y_per <- round(summary(pca_data)$importance[2, 2]*100, 1)
metabolome_pca <- data.frame(name=rownames(pca_data$x), pca_data$x) %>%
  left_join(group, by = "name")

ggplot(data=metabolome_pca, 
             aes(x=PC1,y=PC2,fill= species,colour=species)) + 
  geom_point(aes(color=species),size=3)+
  geom_hline(aes(yintercept=0), linetype=3) + 
  geom_vline(aes(xintercept=0), linetype=3) +  
  stat_ellipse(aes(fill = species), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 
  theme_test()

#########################
########-OPLS-DA-########
#########################
library(ropls)

data=read.csv("./metabolome.csv")
data = t(data)
group=c(rep('Species1',5),
        rep('Species2',5))
oplsda <- opls(data, group, predI = 1, orthoI = NA)
vip <- getVipVn(oplsda)
vip<-data.frame(vip)
data_vip <- merge(data, vip, by="row.names", all=T) 
data_vip$species1_mean<-rowMeans(data_vip[2:6])
data_vip$species2_mean<-rowMeans(data_vip[7:11])
data_vip$Foldchange<-data_vip$species1_mean/data_vip$species2_mean
data_vip$LOG2FC<-log2(data_vip$Foldchange)

##################################
########-Metabolome_WGCNA-########
##################################
library(WGCNA)

data=read.csv("./metabolome.csv")
data = t(data)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
SoftThreshold = pickSoftThreshold(data, powerVector = powers, verbose = 5)
par(mfrow = c(1,2))
plot(SoftThreshold$fitIndices[,1], -sign(SoftThreshold$fitIndices[,3])*SoftThreshold$fitIndices[,2],type="n")
text(SoftThreshold$fitIndices[,1], -sign(SoftThreshold$fitIndices[,3])*SoftThreshold$fitIndices[,2],labels=powers,cex=0.9,col="red")
plot(SoftThreshold$fitIndices[,1], SoftThreshold$fitIndices[,5],type="n")
text(SoftThreshold$fitIndices[,1], SoftThreshold$fitIndices[,5],labels=powers,cex=0.9,col="red")
net = blockwiseModules(data, power = SoftThreshold$powerEstimate,
                       maxBlockSize = 1000,TOMType = "unsigned", 
                       minModuleSize = 5,reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       verbose = 3)
					   
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
					
design=model.matrix(~0+ data$type)
colnames(design)=levels(data$type)
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(data, moduleColors)$eigengenes
MEs = orderMEs(MEs0)				
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(data))
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
colnames(design)=levels(as.factor(data$type))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,0.7))

nSamples <- dim(data)[1]
geneModuleMembership <- cor(data, MEs, use = "p")
MMPvalue <- corPvalueStudent(geneModuleMembership, nSamples)
geneSignificanceCor <- cor(data, design, use = "p")
geneSignificanceP <- corPvalueStudent(geneSignificanceCor, nSamples)
column <- paste0("ME", module)
moduleGenes <- names(net$colors)[which(moduleColors == module)]
MM <- abs(geneModuleMembership[moduleGenes, column])
GS <- abs(geneSignificanceCor[moduleGenes, 2])
verboseScatterplot(
  MM, GS,
  abline = TRUE,
  pch = 21,
  cex.main = 1.2,
  cex.lab = 1.2,
  cex.axis = 1.2,
  col = "black",
  bg = module)

##########################
########-CellChat-########
##########################
###The single-cell data from donkey skin should initially undergo conversion to their human homologous genes.
library(CellChat)
library(Seurat)
library(dplyr)
library(ggalluvial)
library(svglite)

Sample<-readRDS("./Sample.rds")
data.input<-GetAssayData(Sample,assay = "RNA",slot = "data")
meta<-Sample@meta.data
cellchat <- createCellChat(object = data.input,meta=meta,group.by = "celltype")
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human
colnames(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
counts<-cellchat@net$count
weight<-cellchat@net$weight
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F)
				 
netVisual_circle(cellchat@net$counts,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F)
				 
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

###########################
########-nichenetr-########
###########################
###The single-cell data from donkey skin should initially undergo conversion to their human homologous genes.
library(nichenetr)
library(Seurat)

Sample<-readRDS("./Sample.rds")
Idents(Sample) <- "celltype"
ligand_target_matrix = readRDS('./ligand_target_matrix_nsga2r_final.rds')
lr_network = readRDS('./lr_network_human_21122021.rds')
weighted_networks = readRDS('./weighted_networks_nsga2r_final.rds')
receiver = "Sebocyte"
expressed_genes_receiver = get_expressed_genes(receiver, Sample, pct = 0.10,assay_oi = 'RNA')
background_expressed_genes = rownames(ligand_target_matrix)
sender_celltypes = c("Endothelial cell","Fibroblast", "Keratinocyte", "Macrophage", "Mast cell", "Melanocyte", 'Pericyte','Smooth muscle cell','Sweat gland epithelial cell')
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, Sample, 0.10,assay_oi = 'RNA')
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

seurat_obj_receiver= subset(Sample, idents = receiver)
geneset_oi = expressed_genes_receiver
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
best_upstream_ligands = ligand_activities %>% top_n(30, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 1000) %>% bind_rows() %>% drop_na()
write.table(active_ligand_target_links_df,"./active_ligand_target_links_df.txt", quote = FALSE)

########################
########-GENIE3-########
########################
library(GENIE3)
library(Seurat)

Sample<-readRDS("./Sample.rds")
Idents(Sample) <- "celltype"
Sebocyte<-subset(Sample,idents=c("Sebocyte"))
Sebocyte <- as.matrix(Sebocyte@assays$RNA@counts)
Sample_nichenetr <-read.table("./Sample.txt")
Sample_nichenetr <- subset(Sebocyte, rownames(Sebocyte) %in% Sample_nichenetr)
LinkList <- getLinkList(Sample_nichenetr, threshold=0.01)
write.table(LinkList,"./LinkList.txt", quote = FALSE)
