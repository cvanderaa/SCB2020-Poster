
library(tidyverse)
library(xtable)
library(scp)
library(nipals)
library(SingleCellExperiment)
library(scater)

####---- Generate scpdata table ----####

read.csv("../scpdata/inst/extdata/metadata.csv") %>% 
    select(Title, Description, Species, PublicationDate, NumberAssays) %>%
    mutate(Description = paste0(substr(Description, 1, 30), "...")) %>%
    dplyr::rename(Date = PublicationDate,
           `# Assays` = NumberAssays) %>%
    xtable()


####---- Typical scp pipeline ----####

readSCP(quantTable = quantData, metaTable = metaData,
        channelCol = "Channel", batchCol = "Set") %>%
    zeroIsNA(i = 1:4) %>%
    filterFeatures(~ Reverse != "+" & Potential.contaminant != "+") %>%
    subsetByAssay(dims(.)[1, ] > 150) %>%
    computeSCR(i = 1:3, colDataCol = "SampleType",
               carrierPattern = "Carrier",
               samplePattern = "Macrophage|Monocyte") %>%
    filterFeatures(~ !is.na(.meanSCR) & .meanSCR < 0.1) %>%
    aggregateFeaturesOverAssays(i = 1:3, fcol = "peptide",
                                name = paste0("peptides_", names(.)),
                                fun = robustSummary) %>%
    joinAssays(i = 4:6, name = "peptides") %>%
    computeMedianCV(i = "peptides_filter1", proteinCol = "protein",
                    peptideCol = "peptide", batchCol = "Set") %>%
    normalize(i = "peptides", name = "peptides_norm",
              method = "median", na.rm = TRUE) %>%
    logTransform(i = "peptides_norm", name = "peptides_log",
                 base = 2) %>%
    aggregateFeatures(i = "peptides_log", name = "proteins",
                      fcol = "protein", robustSummary) %>%
    impute(i = "proteins_norm", method = "knn", k = 3) ->
    scp

####---- Batch effects ----####

x <- specht2019v2[["peptides"]]
nbatch <- 3
ncell <- 4
totcell <- nbatch * ncell
x <- x[, 1:(totcell)]
x <- x[do.call(order, 
               lapply(seq(1, totcell, by = ncell), 
                      function(i) rowSums(is.na(assay(test)[, i:(i+ncell-1)])))), ]

png(filename = "figs/batch_effects.png", res = 400, height = 2000, width = 2000)
image(t(assays(test)[[1]][seq(3200, 1, -30), ]), axes = FALSE, 
      xlab = "", ylab = "Features (peptides)", 
      mgp = c(1,0,0), 
      col = colorRampPalette(c(rgb(221/256, 230/256, 163/256), 
                               rgb(122/256, 128/256, 89/256)))(1000),
      main = "Peptide expression data",
      useRaster = TRUE)
marg <- 1 / nbatch / 2
mtext(text = paste0("Batch ", 1:3), side = 1, line = 0, 
      at = marg + 2 * marg * ((1:nbatch) - 1))
dev.off()


####---- Imputation ----####

## Get protein data before imputation
specht2019v2 <- aggregateFeatures(specht2019v2,
                                  i = "peptides",  
                                  name = "proteins_NA", 
                                  fcol = "protein", 
                                  fun = colMedians, na.rm = TRUE)
## Compute and plot PCA (NIPALS)
pcaRes <- nipals(assay(specht2019v2[["proteins_NA"]]),
                 ncomp = 5)
reducedDim(specht2019v2[["proteins_NA"]], "nipals") <- pcaRes$loadings
plotReducedDim(getWithColData(specht2019v2, "proteins_NA"), 
               dimred = "nipals", 
               ncomponents = c(3, 5), 
               colour_by = "SampleType") +
    ggtitle("PCA (NIPALS) on missing data") +
    ggsave(filename = "figs/PCA_missing.png", height = 3, width = 4.5,
           dpi = 300)

## Show batch effect 
plotReducedDim(getWithColData(specht2019v2, "proteins_NA"), 
               dimred = "nipals", 
               ncomponents = c(1, 2), 
               colour_by = "lcbatch") +
    ggsave(filename = "figs/PCA_batch_effect.png", 
           height = 3, width = 4.5, dpi = 300)


## Perform imputation
scp_imputeKNN <- function(obj, i, name = "KNNimputedAssay", k = 3){
    
    oldi <- i
    exp <- obj[[i]]
    dat <- assay(exp)
    
    # Create a copy of the data, NA values to be filled in later
    dat.imp<-dat
    
    # Calculate similarity metrics for all column pairs (default is Euclidean distance)
    dist.mat<-as.matrix( dist(t(dat)) )
    #dist.mat<-as.matrix(as.dist( dist.cosine(t(dat)) ))
    
    # Column names of the similarity matrix, same as data matrix
    cnames<-colnames(dist.mat)
    
    # For each column in the data... 
    for(X in cnames){
        
        # Find the distances of all other columns to that column 
        distances<-dist.mat[, X]
        
        # Reorder the distances, smallest to largest (this will reorder the column names as well)
        distances.ordered<-distances[order(distances, decreasing = F)]
        
        # Reorder the data matrix columns, smallest distance to largest from the column of interest
        # Obviously, first column will be the column of interest, column X
        dat.reordered<-dat[ , names(distances.ordered ) ]
        
        # Take the values in the column of interest
        vec<-dat[, X]
        
        # Which entries are missing and need to be imputed...
        na.index<-which( is.na(vec) )
        
        # For each of the missing entries (rows) in column X...
        for(i in na.index){
            
            # Find the most similar columns that have a non-NA value in this row
            closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )
            
            # If there are more than k such columns, take the first k most similar
            if( length(closest.columns)>k ){
                # Replace NA in column X with the mean the same row in k of the most similar columns
                vec[i]<-mean( dat[ i, closest.columns[1:k] ] )
            }
            
            # If there are less that or equal to k columns, take all the columns
            if( length(closest.columns)<=k ){
                # Replace NA in column X with the mean the same row in all of the most similar columns
                vec[i]<-mean( dat[ i, closest.columns ])
            }
        }
        # Populate a the matrix with the new, imputed values
        dat.imp[,X]<-vec
    }
    
    assay(exp) <- dat.imp
    obj <- addAssay(obj, exp, name = name)
    addAssayLinkOneToOne(obj, from = oldi, to = name)
}
specht2019v2 <- scp_imputeKNN(specht2019v2,
                              i = "proteins_NA", 
                              name = "proteins_impd", 
                              k = 3)
## Compute and plot PCA (NIPALS)
pcaRes <- nipals(assay(specht2019v2[["proteins_impd"]]),
                 ncomp = 5)
reducedDim(specht2019v2[["proteins_impd"]], "nipals") <- pcaRes$loadings
plotReducedDim(getWithColData(specht2019v2, "proteins_impd"), 
               dimred = "nipals", 
               ncomponents = c(3, 5),
               colour_by = "SampleType") +
    theme(legend.position = "none") +
    ggtitle("PCA (NIPALS) after imputation") +
    ggsave(filename = "figs/PCA_imputed.png", height = 3, width = 3,
           dpi = 300)

