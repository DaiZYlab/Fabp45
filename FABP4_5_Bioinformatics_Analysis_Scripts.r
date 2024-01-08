### Analysis of FABP4/5 Paper

## BulkRNAseq analysis
#If there are any modules loaded, remove them.
module purge
module load samtools/1.6
module load pigz/2.4
module load boost/1.56.0
module load gcc/4.8.3

# Running FastQ_screen to look for contamination.
module load bowtie2/2.2.6
module load perl/5.16

# Check all reference genomes currently installed
perl ~/bin/fastq_screen_v0.11.4/fastq_screen --threads 12 --aligner bowtie2 --conf ~/bin/fastq_screen_v0.11.4/fastq_screen.allRefs.conf --outdir ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/fastqc ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/fastq/WT_Rep1_R1.fastq.gz

# Setting up trimming for sample WT_Rep1 for paired end reads.
# Trim PE poor quality sequence with TRAILING:30 MINLEN:20 (see Trimmomatic documentation)
java -jar ~/bin/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 12 -phred33 ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/fastq/WT_Rep1_R1.fastq.gz ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/fastq/WT_Rep1_R2.fastq.gz /projects/b1042/PanosLab/Yalu/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/fastq/WT_Rep1_R1.fastq.trimmed.gz /projects/b1042/PanosLab/Yalu/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/fastq/WT_Rep1_R1U.fastq.trimmed.gz /projects/b1042/PanosLab/Yalu/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/fastq/WT_Rep1_R2.fastq.trimmed.gz ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/fastq/WT_Rep1_R2U.fastq.trimmed.gz TRAILING:30 MINLEN:20

./bin/FastQC/fastqc -t 6 -o ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/fastqc ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/fastq/WT_Rep1*.fastq.gz ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/fastq/WT_Rep1*.fastq.trimmed.gz
date

# Align reads with STAR
module load STAR/2.5.2
# Adding Readgroups from rgstring ID:WT_Rep1 LB:WT_Rep1 PU:nextseq DT:2023-12-05 SM:WT_Rep1 CN:ASH PL:illumina
STAR --runMode alignReads --genomeDir ~/anno/STAR_indexes/mm10/ --runThreadN 12 --readFilesIn ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/fastq/WT_Rep1_R1.fastq.trimmed.gz ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/fastq/WT_Rep1_R2.fastq.trimmed.gz --readFilesCommand zcat -c --outFileNamePrefix ./Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/STAR_aln/WT_Rep1 --outSAMtype BAM Unsorted --chimSegmentMin 20 --quantMode TranscriptomeSAM --outReadsUnmapped Fastq --outMultimapperOrder Random --outSAMattrRGline ID:WT_Rep1 LB:WT_Rep1 PU:nextseq DT:2023-12-05 SM:WT_Rep1 CN:ASH PL:illumina --outSAMmapqUnique 60 --outFilterMultimapNmax 20 --outFilterMismatchNmax 2

# Sort the output of STAR (outputting sorted BAMs from STAR took too much memory)
samtools sort -@ 12 -o ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/STAR_aln/WT_Rep1Aligned.sortedByCoord.out.bam ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/STAR_aln/WT_Rep1Aligned.out.bam
date
ln -s ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/STAR_aln/WT_Rep1Aligned.sortedByCoord.out.bam ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/bam/WT_Rep1.bam
date
cd ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/bam
samtools index -@ 12 ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/bam/WT_Rep1.bam
module unload mpi
module load python/anaconda3

# Run htseq-count for WT_Rep1, stranded = 1
htseq-count -f bam -q -m intersection-nonempty -s reverse -t exon -i gene_id ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/bam/WT_Rep1.bam ~/anno/Ens/mm10.Ens_78/mm10.Ens_78.cuff.gtf > ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/bam/WT_Rep1.htseq.counts

# Wait for htseq jobs to finish, then start more.
wait
module unload python/anaconda3
module load gcc/4.8.3
module load R/3.2.2

# Make Analysis directory for all analysis files.
mkdir ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/analysis

# Make HTseq counts table.
perl ~/tools/bin/makeHTseqCountsTable.pl ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/bam// ~/anno/Ens/mm10.Ens_78/mm10.Ens_78.cuff.gtf ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/bam//

# Run EdgeR, using comparisons file, without MDS plot (which sometimes crashes), htseq.
Rscript ~/tools/bin/runEdgeRrnaSeq.2.R --assembly=mm10 --countFile= ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/bam//htseq.all.counts.txt --comparisonFile= ~/Data/FABP4_Dai/bulk_RNAseq/pipeline/comparison.csv --numCores=1 --outputDirectory= ~/Data/FABP4_Dai/bulk_RNAseq/pipeline//bulk_RNAseq/analysis --runMDS=0


## PCA plot
library(ggplot2)
exprsFile <- "~/Data/FABP45_bulkRNAseq/analysis/htseq.normCounts.txt"
data.matrix <- as.matrix(read.table(exprsFile, header = TRUE, sep="\t", row.names = 1, as.is = TRUE))

# PCA plot of gene expression
pca <- prcomp(t(data.matrix), scale = TRUE)
plot(pca$x[,1], pca$x[,2])
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
barplot(pca.var.per, main = "Scree Plot", xlab = "Principal Component", ylab = "Percent Variation")

pca.data <- data.frame(Sample=rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])
pca.data

pca.data$color <- factor(pca.data$color, level = c("red", "blue", "green", "purple"))
colr <- as.character(unique(pca.data$color))
colr

ggplot(pca.data, aes(X,Y))+geom_point(aes(color=color))+ggtitle("PCA plot")+ylim(-80,60)+xlim(-120,110)+xlab("PC1")+ylab("PC2")+theme_bw()+theme(panel.border=element_rect(size=1), 
panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(colour="black"), legend.title=element_blank(), legend.box.background=element_rect(size=0))
+scale_color_manual(labels=c(paste("CKO"),paste("TKO"),paste("WT")), values=colr)


## Complex Heatmap
# Load packages
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyverse)

# Load and modify data
exprsFile <- "~/Data/FABP45_bulkRNAseq/Heatmaps/Heatmap_PH_Glycolysis_related_genes.txt"
exprs <- read.table(exprsFile, header = TRUE, sep = "\t", row.names = 1)

data <- exprs[,2:10]
head(data)

# Scale data by row
mat_scaled = t(scale(t(data)))

# Generate row annotation
Hallmark <- exprs[,1]
activity=list(Hallmark=c("PH" = "#d07b70", "Glycolysis" = "#5fb6c1"))
ha_d = HeatmapAnnotation(Hallmark=Hallmark, which = "row", width = unit(4, "mm"), col = activity)

# Draw Heatmaps
Heatmap(mat_scaled, name = "Z-score", col=brewer.pal(11,"RdBu"), left_annotation = ha_d, show_row_names = TRUE, row_order = 28:1, row_names_gp = gpar(fontsize = 10), column_order= 1:9)

## Heatmap with pheatmap
# Load package
library(pheatmap)

## Load data
# Load expression data
exprsFile <- "~/Data/FABP45_bulkRNAseq/Heatmap_PKO_vs_WT_adj_0.5_avg_reads_30.txt"
exprs <- as.matrix(read.table(exprsFile, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE))
exp <- log2(exprs+1)

# Load metadata
pData <- "~/Data/FABP45_bulkRNAseq/heatmap_metadata.txt"
annotation <- read.table(pData, header = TRUE, sep = "\t", row.names = 1)

## define color
hmcol <- colorRampPalette(c("blue", "white", "red"))(256)
ann_colors = list(sampletype = c(WT="blue", PKO="red", TKO="green"))
pheatmap(exp, color = hmcol, annotation = annotation, cluster_cols = F, cluster_rows = T, show_rownames=F, annotation_colors= ann_colors, fontsize = 10, scale="row",fontsize_row = 10, height=20)


## Dotplot enrichment analysis
# Load library and data
library(ggplot2)

exprsFile <- "~/Data/FABP45_scRNAseq/Analysis/Pathway_Enrichment/Bubble_Chart_Hallmark_Up_Genes_EC_CKO_vs_WT.txt"
data <- read.table(exprsFile, sep ="\t", header = TRUE, stringsAsFactors = FALSE)

S1<- ggplot(data, aes(x= NES, y=reorder(Pathway,-FDR), size=DEGs, color=FDR)) + geom_point(alpha = 0.8) 
S1 = S1+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab", limit = c(9e-19, 3.20e-05))+scale_size(range = c(2, 8))
S1+theme_bw(base_size = 24) +
    theme(
        legend.position = 'right',
        legend.background = element_rect(),
        plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
        plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
        plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
        
        axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
        axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
        axis.title = element_text(size = 12, face = 'bold'),
        axis.title.x = element_text(size = 12, face = 'bold'),
        axis.title.y = element_text(size = 12, face = 'bold'),
        axis.line = element_line(colour = 'black'),
        
        #Legend
        legend.key = element_blank(), # removes the border
        legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
        legend.text = element_text(size = 14, face = "bold"), # Text size
        title = element_text(size = 14, face = "bold")) +xlab('Enrichment score') + ylab('')

## Enrichment Bar Charts with ggplot2
# load libraries
library(ggplot2)

# load data
exprsFile <- "~/Data/FABP45_scRNAseq/Analysis/Pathway_Enrichment/Bar_Chart_GO_Dn_Genes_EC_CKO_vs_WT.txt"
df <- read.table(exprsFile, header=TRUE, as.is=TRUE, sep = "\t")

# with unique color
ggplot(df, aes(x=reorder(Term_Description,adjust_p), y=adjust_p))+geom_bar(stat="identity", fill="blue")+geom_col(position="dodge", width=0.4, size=0.7)+coord_flip()+xlab("")+ylab("-log10(adj_p)")+theme_bw()


## Volcano Plot
# Load package and file
library(EnhancedVolcano)
exprsFile <- "~/Data/FABP45_bulkRNAseq/Volcano_Plot_DEGs_TKO_vs_WT.txt"
res1 <- read.table(exprsFile, row.names = 1, header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Setup cutoff
FC <- 0.585
p <- 0.05

# Define colors
keyvals <- rep('grey75', nrow(res1))
names(keyvals) <- rep('NS', nrow(res1))

keyvals[which(abs(res1$log2FoldChange) < FC & res1$pvalue < p)] <- 'grey75'
names(keyvals)[which(abs(res1$log2FoldChange)  < FC & res1$pvalue < p)] <- 'NS'

# Down-regulated
keyvals[which(res1$log2FoldChange < -FC & res1$pvalue < p)] <- 'blue2'
names(keyvals)[which(res1$log2FoldChange  < -FC & res1$pvalue < p)] <- 'Significantly down-regulated'

#Up-regulated
keyvals[which(res1$log2FoldChange > FC & res1$pvalue < p)] <- 'red2'
names(keyvals)[which(res1$log2FoldChange > FC & res1$pvalue < p)] <- 'Significantly up-regulated'

# Check key values
unique(keyvals)
unique(names(keyvals))

# Draw plot
EnhancedVolcano(res1,
+                 lab = NA,
+                 x = 'log2FoldChange',
+                 y = 'pvalue',
+                 #selectLab = rownames(res2)[which(names(keyvals) %in% c('NS','log2FoldChange','-Log10Q','low','high'))],
+                 xlim = c(-12,12),
+                 ylim = c(-0.5,7120),
+                 xlab = bquote(~Log[2]~ 'fold change'),
+                 ylab = bquote(~-Log[10] ~ italic(P)),
+                 title = 'CKO_vs_WT',
+                 pCutoff = 0.05,
+                 FCcutoff = 0.585,
+                 pointSize = 2.5,
+                 labSize = 4.5,
+                 colCustom = keyvals,
+                 colAlpha = 0.75,
+                 legendPosition = 'right',
+                 legendLabSize = 15,
+                 legendIconSize = 5.0,
+                 drawConnectors = FALSE,
+                 widthConnectors = 0.5,
+                 colConnectors = 'grey75',
+                 gridlines.major = TRUE,
+                 gridlines.minor = FALSE,
+                 border = 'partial',
+                 borderWidth = 1.5,
+                 borderColour = 'black')


### scRNAseq analysis
## Load libraries
library(Seurat)
library(SeuratData)
library(patchwork)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(scDblFinder)
library(SingleCellExperiment)

## scRNAseq data were processed followed the Seurat analysis pipeline.

## Violin Plot
DefaultAssay(lung.combined) <- "RNA"
lung.combined <- NormalizeData(lung.combined)

VlnPlot(lung.combined, features = "Fabp4")
VlnPlot(lung.combined, features = "Fabp5")

## Glycolysis Score by violin plot
glycolysis.genes <- c("Eno1", "Ldha", "Pgk1", "Pgm2", "Ppf1a4", "Slc16a3", "Tpi1", "Pgam1", "Pkm", "Artn", "Adora2b", "Pfkp")
lung.combined <- AddModuleScore(lung.combined, features = list(glycolysis.genes), name="Glycolysis_Score")
VlnPlot(lung.combined, features = "Glycolysis_Score", group.by = "orig.ident")
VlnPlot(lung.combined, features = "Glycolysis_Score", split.by = "orig.ident")
DoHeatmap(ec.combined, features = glycolysis.genes)

## EC subclusters
Idents(lung.combined) <- "celltype"
lung.ec <- subset(lung.combined, idents = "EC")
DefaultAssay(lung.ec) <- "RNA"
??? ec.list <- split(lung.ec, f = lung.ec$orig.ident)

ec.list <- lapply(X = ec.list, FUN = function(x) {
          x <- NormalizeData(x)
          x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
      })

features <- SelectIntegrationFeatures(object.list = ec.list)
ec.anchors <- FindIntegrationAnchors(object.list = ec.list, anchor.features = features)
ec.combined <- IntegrateData(anchorset = ec.anchors)

## Dimention reduction and clustering
DefaultAssay(ec.combined) <- "integrated"
ec.combined <- ScaleData(ec.combined, verbose = FALSE)
ec.combined <- RunPCA(ec.combined, npcs = 30, verbose = FALSE)
DimPlot(ec.combined, reduction = "pca")
ElbowPlot(ec.combined, ndims = 30)
ec.combined <- RunUMAP(ec.combined, reduction = "pca", dims = 1:15)
ec.combined <- FindNeighbors(ec.combined, reduction = "pca", dims = 1:15)
ec.combined <- FindClusters(ec.combined, resolution = 0.1)
DimPlot(ec.combined, reduction = "umap", split.by = "orig.ident")

## tacked Bar Plot for cluster proportions in scRNAseq
library(ggplot2)
library(RColorBrewer)

propFile <- "~/Data/FABP45_ecRNAseq/EC_Subclusters_Proportion.txt"
prop <- read.table(propFile, header = TRUE, sep = "\t")
prop$Clusters <- as.character(prop$Clusters)

# Draw Stacked Plot
ggplot(prop, aes(x=orig.ident, y = Freq, fill = Clusters)) 
+ theme_bw(base_size = 15) + geom_col(position = "fill", width = 0.5) 
+ xlab("Sample") + ylab("Proportion") 
+ scale_fill_manual(values = brewer.pal(12, "Paired")) 
+ theme(legend.title = element_blank())

## Single cell proportion test
library("scProportionTest")
prop_test <- sc_utils(ec.combined)
prop_test <- permutation_test(
	prop_test, cluster_identity = "ec.sub",
	sample_1 = "WT", sample_2 = "CKO",
	sample_identity = "orig.ident"
)
permutation_plot(prop_test)

## Dotplot 
markers.to.plot <- c("Cdh5", "Pecam1", "Emcn", "Aqp5", "Ager", "Rtkn2", "Sftpb", "Lamp3", "Notch3",
    "Acta2", "Myh11", "Cnn1", "Fn1", "Col1a1", "Gm867", "Cdhr3", "Scgb3a2", "Olr1", "F7", "Atp6v0d2", 
    "C1qa", "C1qb", "Ms4a7", "Ifitm6", "Gpr141", "Xcr1", "Sept3", "S100a9", "S100a8", "Mmp9", "Mmp9", 
    "Ighm", "Cd79a", "Cd3g", "Cd3e", "Trac", "Klrb1c", "Klre1", "Ncr1", "Ccl21a", "Prox1", "Reln", "Top2a", "Mki67")

DotPlot(lung.combined, features = markers.to.plot, dot.scale = 8) + RotatedAxis()

## EC markers
gCapillary.enriched  <- c("Aplnr","Plvap", "Gpihbp1","Socs3", "Cd93", "Sema3c")
lung.combined <- AddModuleScore(lung.combined, features = list(gCapillary.enriched), name="gCapillary_enriched")
FeaturePlot(ec.combined, features = "gCapillary_enriched", label = TRUE, repel = TRUE) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

aCapillary.enriched  <- c("Chst1", "Tbx2", "Car4", "Igfbp7", "Fibin", "Ednrb", "Tbx2", "Cdkn2b", "Ptprb")
lung.combined <- AddModuleScore(lung.combined, features = list(aCapillary.enriched), name="aCapillary_enriched")
FeaturePlot(ec.combined, features = "aCapillary_enriched", label = TRUE, repel = TRUE) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

Artery.enriched  <- c("Mgp", "Cxcl12", "Tm4sf1", "Plat", "Gja4", "Plac8", "Mecom", "Atp13a3", "Gja5")
lung.combined <- AddModuleScore(lung.combined, features = list(Artery.enriched), name="Artery_enriched")
FeaturePlot(ec.combined, features = "Artery_enriched", label = TRUE, repel = TRUE) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

Vein.enriched  <- c("Prss23", "Slc6a2", "Vwf", "Cd200", "Bgn", "Ackr3", "Lyve1", "Tmem252", "Amigo2", "Vcam1")
lung.combined <- AddModuleScore(lung.combined, features = list(Vein.enriched), name="Vein_enriched")
FeaturePlot(ec.combined, features = "Vein_enriched", label = TRUE, repel = TRUE) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))




