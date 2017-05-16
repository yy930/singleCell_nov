library(monocle)
library(reshape2)
celltype <- read.delim("~/Documents/sc_count/celltype_all.txt")
#celltype_T <- celltype[celltype$cell.type== "T",]
#genecount <- read.delim("~/Documents/sc_count/agg_count.txt", row.names=1,check.names=FALSE)
genetpm <- read.delim("~/Documents/sc_count/agg_tpm.txt", row.names=1,check.names=FALSE)
genetpm<-data.matrix(genetpm)
#genetpm_T <-genetpm[fData(QCmet)$use_man,(celltype$cell.type== "T")&(pData(QCmet)$use_man)]
gene_biotype <- read.delim("~/Documents/sc_count/gene_biotype.txt",row.names=1,check.names=FALSE)
gene_biotype <- gene_biotype[rownames(genetpm),]
colnames(gene_biotype) <- c("gene_short_name","biotype")
#filter_genes_T <- apply(genetpm_T, 1, 
#                      function(x) length(x[x > 1]) >= 2)
#genetpm_T<-genetpm_T[filter_genes_T,]
#gene_biotype_T<-gene_biotype[filter_genes_T,]
###############
ercc_gene <- rownames(gene_biotype)[grep("^ERCC", rownames(gene_biotype), perl=TRUE, value=FALSE)]
endo_biotype<-gene_biotype[-which(rownames(gene_biotype) %in% ercc_gene),]
ind_ercc <- grep("^ERCC", rownames(gene_biotype), perl=TRUE, value=FALSE)
#endo_biotype<-gene_biotype[-ind_ercc,]
ercctpm<-genetpm[which(rownames(genetpm) %in% ercc_gene),]
#endotpm_T<-genetpm[-ind_ercc,celltype$cell.type== "T"]
endotpm<-genetpm[-which(rownames(genetpm) %in% ercc_gene),]
library("readr")
ERCC_annotation<- as.data.frame(read_delim("~/Documents/NSC/concentrationERCC.txt",
                             "\t", escape_double = FALSE, trim_ws = TRUE))
colnames(ERCC_annotation)[4]<-"conc_attomoles_ul_Mix1"
colnames(ERCC_annotation)[5]<-"conc_attomoles_ul_Mix2"
ERCC_annotation<-ERCC_annotation[match(rownames(ercctpm),ERCC_annotation$`ERCC ID`),]
rownames(ERCC_annotation)<-ERCC_annotation$`ERCC ID`

abs_ercc<-8*ERCC_annotation$conc_attomoles_ul_Mix1/800000000
mean_ercc<-rowMeans(ercctpm)
#correlation
corr<-apply(ercctpm,2,function(x) cor(log(x)[x>0],log(abs_ercc)[x>0]))
#logistic regression model of each sample's detection limit
#the probability of detecting a spike-in at a given input level was modeled by the logistic function:
logercctpm<-log(cbind.data.frame(mean_ercc,abs_ercc))[-which(mean_ercc==0),]
lercctpm <- replace(ercctpm,ercctpm>0,1)
lercctpm<-lercctpm[,-which(colSums(lercctpm)<=8)]
#lercctpm<-cbind(lercctpm,abs_ercc)
library(reshape2)
lercctpm<-melt(lercctpm)
lercctpm<-cbind(lercctpm,rep(log(abs_ercc),478))
colnames(lercctpm)<-c("ercc_id","cell_id","detected_ercc","log_abs_ercc")
model <- glm(detected_ercc ~ log_abs_ercc,family=binomial(link='logit'),data=lercctpm)
library(ggplot2)
ggplot(lercctpm,aes(x=log_abs_ercc,y=detected_ercc))+geom_point()+
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)

#plot for pearson correlation between absolute number of ercc & ercctpm 
scatter_plot <- ggplot(logercctpm, aes(abs_ercc,mean_ercc))
scatter_plot + geom_point() + labs(x = "log input ERCC molecules", y = "log mean TPM") 
#voilin plot for accuracy and Sensitivity 
#correlation (accuracy)
corr<-apply(ercctpm,2,function(x) cor(log(x)[x>0],log(abs_ercc)[x>0]))
detection_limit<-apply(ercctpm,2,function(x) abs_ercc[which.min(x[x>0])])
mat<-cbind(corr,detection_limit)
mat<-melt(mat)
p<-ggplot(mat, aes(x=Var2, y=value))
+ geom_dotplot(binaxis='y', stackdir='center', dotsize=1)?????????
# Calculating Pearson's product-moment correlation
cor.test(logercctpm$abs_ercc, logercctpm$mean_ercc, method = "pearson", conf.level = 0.95)

#pd <- new("AnnotatedDataFrame", data = celltype_T)
pd <- new("AnnotatedDataFrame", data = celltype)
gene_biotype_T<-gene_biotype[rownames(endotpm_T),]
#fd <- new("AnnotatedDataFrame", data = gene_biotype_T)
fd <- new("AnnotatedDataFrame", data = gene_biotype)#[-which(rownames(genetpm) %in% ercc_gene),])
# # First create a CellDataSet from the relative expression levels
HSMM <- newCellDataSet(as.matrix(genetpm),#endotpm_T
                        phenoData = pd,
                        featureData = fd)
rpc_matrix <- relative2abs(HSMM)

# Now, make a new CellDataSet using the RNA counts
HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit=1,
                       expressionFamily=negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

# # Next, use it to estimate RNA counts #HSMM including ERCC
#rpc_matrix <- relative2abs(HSMM, ERCC_controls =ercctpm[,celltype$cell.type== "T"],ERCC_annotation=ERCC_annotation)
rpc_matrix_byercc<-relative2abs(HSMM, t_estimate = estimate_t(genetpm),
             modelFormulaStr = "~1", ERCC_controls = ercctpm, ERCC_annotation = ERCC_annotation,
             volume = 8000, dilution = 800000000, mixture_type = 1,
             detection_threshold = 800, expected_capture_rate = 0.25,
             verbose = FALSE, return_all = FALSE, cores = 1)

rownames(rpc_matrix_byercc)<-rownames(genetpm)
colnames(rpc_matrix_byercc)<-colnames(genetpm)
#write.table(rpc_matrix, file = "~/Documents/sc_count/rpc_matrix.txt",sep="\t", row.names = FALSE, qmethod = "double")
write.table(rpc_matrix, file = "~/Documents/sc_count/agg_rpc_matrix.txt",sep="\t", row.names = FALSE, qmethod = "double")

# # Now, make a new CellDataSet using the RNA counts
HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit=1,
                       expressionFamily=negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
#2.5 Filtering low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 2))

pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) + 2*sd(log10(pData(HSMM)$Total_mRNAs))) 
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) - 2*sd(log10(pData(HSMM)$Total_mRNAs))) 
qplot(Total_mRNAs, data=pData(HSMM), color=batch, geom="density") +
  geom_vline(xintercept=lower_bound) +
  geom_vline(xintercept=upper_bound)

#valid_cells <- QCmet$cell
#drop.lib=!(pData(HSMM)$Total_mRNAs > 10000 & pData(HSMM)$Total_mRNAs < upper_bound)
#id.drop.lib=which(drop.lib==TRUE)
id.drop<-list(id.drop.pct.mapped,id.drop.n.reads,id.drop.n.genes,
              id.drop.pct.mt,id.drop.pct.ercc,id.drop.control,id.drop.t.outlier)
x <- expand.grid(seq(length(id.drop)),seq(length(id.drop)))
dropmat <- mapply(function(m1, m2) length(intersect(id.drop[[m1]],id.drop[[m2]])) 
                  , x[, 1] 
                  , x[, 2] )
dropmat<- matrix(dropmat, nrow = length(id.drop)) 
print(dropmat)

HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
               pData(HSMM)$Total_mRNAs < upper_bound]
#HSMM <- HSMM[,!(drop.pct.mapped|drop.n.reads|drop.n.genes|drop.pct.mt|drop.pct.ercc|drop.control|drop.t.outlier)[celltype$cell.type== "T"]]
HSMM <- detectGenes(HSMM, min_expr = 0.1)
#expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 2))
HSMM <- HSMM[expressed_genes,]
#HSMM <- HSMM[expressed_genes,(celltype$cell.type== "T")&(pData(QCmet)$use_man)]
dim(HSMM)# 23602      134  

#3.1 Classifying cells with CellTypeHierarchy
TH1_id <- row.names(subset(fData(HSMM), gene_short_name %in% c("TBX21","IFNG")))#marker for Th1:"IL2","IL12","TNFA"
TH2_id <- row.names(subset(fData(HSMM), gene_short_name %in% c("GATA3","EPAS1"))) #Th2:"BATF3","IL3","IL4","IL5"
TH17_id <- row.names(subset(fData(HSMM), gene_short_name =="RORC")) #Th17:"IL17", "IL22","IL23"
TFH_id <- row.names(subset(fData(HSMM), gene_short_name =="IL21"))
TREG_id <- row.names(subset(fData(HSMM), gene_short_name %in% c("FOXP3", "IL10")))#Treg:"TGFB",
cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Th1", classify_func=function(x) {colSums(x[TH1_id,] > 1)>=1})
cth <- addCellType(cth, "Th1", classify_func=function(x) {x[TH1_id,] >=1})
cth <- addCellType(cth, "Th2", classify_func=function(x) {colSums(x[TH2_id,] > 1)>=1})
cth <- addCellType(cth, "Th2", classify_func=function(x) {x[TH2_id,] >=1})
cth <- addCellType(cth, "Th17", classify_func=function(x) {x[TH17_id,] > 1})
cth <- addCellType(cth, "Tfh", classify_func=function(x) {x[TFH_id,] > 1})
cth <- addCellType(cth, "Treg", classify_func=function(x) {x[TREG_id,] > 1})

HSMM <- classifyCells(HSMM, cth, 0.1)
table(pData(HSMM)$CellType)

pie <- ggplot(pData(HSMM), aes(x = factor(1), fill = factor(CellType))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

#Unsupervised cell clustering
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.5 & dispersion_empirical >= 2 * dispersion_fit)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)#2484 genes
plot_ordering_genes(HSMM)
#clustering the cells:
HSMM <- clusterCells(HSMM, num_clusters=2)
plot_cell_trajectory(HSMM, 1, 2, color="CellType")
#Semi-supervised cell clustering with known marker genes
marker_diff <- markerDiffTable(HSMM[unsup_clustering_genes,],
                               cth,
                               residualModelFormulaStr="~ batch",
                               cores=1)
candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.01))
marker_spec <- calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth) 
head(selectTopMarkers(marker_spec, 3))
semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 20)$gene_id) 
HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes) 
plot_ordering_genes(HSMM)

HSMM <- clusterCells(HSMM,
                     num_clusters=5,
                     clustering_genes=semisup_clustering_genes,
                     residualModelFormulaStr="~batch")
plot_cell_trajectory(HSMM, 1, 2, color="CellType", markers = c("TBX21"))
plot_cell_trajectory(HSMM, 1, 2, color="CellType")
plotPhenoData(QCmet, aesth = aes_string(x = "total_features",
                                        y = "n_detected_feature_controls_ERCC", colour = "batch"))
