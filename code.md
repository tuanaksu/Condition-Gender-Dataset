```R
directory <- "C:/Users/haksu/Desktop/Adebali Lab/RUNRUNRUN"
sampleFiles <- grep("SRR",list.files(directory),value=TRUE)
sampleCondition <- sub(".*count_","\\1",sampleFiles)
sampleCondition <- sub("_.*","\\1",sampleCondition)
sampleName <- sub("_count.*","\\1",sampleFiles)
sampleDataset <- sub(".*control_","\\1",sampleFiles)
sampleDataset <- sub(".*md_","\\1",sampleDataset)
sampleDataset <- sub(".*suicidal_","\\1",sampleDataset)
sampleDataset <- sub("_.*","\\1",sampleDataset)
sampleGender <- sub(".*_","\\1",sampleFiles)
sampleGender <- sub(".txt","\\1",sampleGender)
```

```R
sampleTable <- data.frame(sampleName = sampleName,fileName = sampleFiles,condition = sampleCondition, sampleDataset = sampleDataset, sampleGender=sampleGender)
sampleTable$condition <- factor(sampleTable$condition)
sampleTable$sampleDataset <- factor(sampleTable$sampleDataset)
sampleTable$sampleGender <- factor(sampleTable$sampleGender)
```

```R
library("DESeq2")
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = directory,design= ~ sampleDataset + sampleGender + condition)
dds <- DESeq(dds)
results <- results(dds)
View(data.frame(results))
```

```R
vsd <- vst(dds, blind = T)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$sampleDataset)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$sampleGender)
plotPCA(vsd, intgroup=c("condition","sampleDataset","sampleGender"))
```

```R
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat) 
heatmap(vsd_cor)
```

```R
contrast <- c("condition","md","control")
results(dds, contrast = contrast, alpha = 0.05)
res_tableOE_unshrunken <- results(dds, contrast=contrast, alpha = 0.05)
res_tableOE <- lfcShrink(dds, contrast=contrast, type = "normal",  res=res_tableOE_unshrunken)
```

```R
padj.cutoff <- 0.05
res_tableOE_df <- data.frame(res_tableOE)
res_tableOE_df$ensembl_gene_id <- rownames(res_tableOE_df)
rownames(res_tableOE_df) <- NULL
sig <- subset(res_tableOE_df, padj <= padj.cutoff)
```

```R
FUN <- function(x){
    x$oe = ifelse(x$padj <= 0.05, TRUE, FALSE)
    return(x)
}
library("dplyr")
library("ggplot2")
res_tableOE_volcano <- res_tableOE_df %>% mutate(oe = padj < 0.05)
ggplot(res_tableOE_volcano) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = oe)) +
    ggtitle("Control vs MD (Gender-Specific)") +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))
```
```R
> library("AnnotationDbi")
> library("org.Hs.eg.db")
> columns(org.Hs.eg.db)
> res_tableOE$symbol <- mapIds(org.Hs.eg.db,
+                          keys=row.names(res_tableOE), 
+                          column="SYMBOL",
+                          keytype="ENSEMBL",
+                          multiVals="first")
> res_tableOE$entrez = mapIds(org.Hs.eg.db,
+                             keys=row.names(res_tableOE), 
+                             column="ENTREZID",
+                             keytype="ENSEMBL",
+                             multiVals="first")
> res_tableOE$name =   mapIds(org.Hs.eg.db,
+                             keys=row.names(res_tableOE), 
+                             column="GENENAME",
+                             keytype="ENSEMBL",
+                             multiVals="first")
> library(pathview)
> library(gage)
> library(gageData)
> data(kegg.sets.hs)
> data(sigmet.idx.hs)
> kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
> foldchanges =res_tableOE$log2FoldChange
> names(foldchanges) = res_tableOE$entrez
> keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
> keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
+     tbl_df() %>% 
+     filter(row_number()<=5) %>% 
+     .$id %>% 
+     as.character()
> keggresids = substr(keggrespathways, start=1, stop=8)
> plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
> tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))
```
