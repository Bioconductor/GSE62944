### =========================================================================
### GSE62944 metadata 
### -------------------------------------------------------------------------
###

meta <- data.frame(
    Title = c(paste0("RNA-Sequencing and clinical data for 7706 tumor ",
                   "samples from The Cancer Genome Atlas"),
              paste0("RNA-Sequencing and clinical data for 9246 tumor ",
                   "samples from The Cancer Genome Atlas"), 
              paste0("RNA-Sequencing and clinical data for 741 normal ",
                   "samples from The Cancer Genome Atlas")),
    Description = c(paste0("TCGA RNA-seq Rsubread-summarized raw count data ",
                         "for 7706 tumor samples, represented as an ", 
                         "ExpressionSet. R data representation derived from ",
                         "GEO accession GSE62944."), 
                    paste0("TCGA RNA-seq Rsubread-summarized raw count data ",
                         "for 9246 tumor samples, represented as a ",
                         "SummarizedExperiment. R data representation derived from ",
                         "GEO accession GSE62944."), 
                    paste0("TCGA RNA-seq Rsubread-summarized raw count data ",
                         "for 741 normal samples, represented as a ",
                         "SummarizedEXperiment. R data representation derived from ",
                         "GEO accession GSE62944.")),
    BiocVersion = c("3.2", "3.3", "3.3"),
    Genome = rep("hg19", 3), 
    SourceType = rep("tar.gz", 3), 
    SourceUrl = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944",
    SourceVersion = "Jan 28 2015",
    Species = "Homo sapiens",
    TaxonomyId = 9606,
    Coordinate_1_based = TRUE,
    DataProvider = "GEO",
    Maintainer = "Bioconductor Package Maintainer <maintainer@bioconductor.org>",
    RDataClass = c("ExpressionSet", rep("SummarizedExperiment", 2)) ,
    DispatchClass = c("ExpressionSet", rep("SummarizedExperiment", 2)),
    ResourceName = c(paste0("GSE62944_GSM1536837_TCGA_20.Illumina.",
                          "tumor_Rsubread_FeatureCounts.ExpressionSet.Rda"), 
             "GSE62944_GSM1536837_TCGA_24.tumor_Rsubread_FeatureCounts.SE.Rda", 
             "GSE62944_GSM1697009_TCGA_24.normal_Rsubread_FeatureCounts.SE.Rda")
)

write.csv(meta, file="metadata.csv", row.names=FALSE)
