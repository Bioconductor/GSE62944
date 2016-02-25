### =========================================================================
### GSE62944 metadata 
### -------------------------------------------------------------------------
###

meta <- data.frame(
    Title = paste0("RNA-Sequencing and clinical data for 7706 tumor ",
                   "samples from The Cancer Genome Atlas"),
    Description = paste0("TCGA RNA-seq Rsubread-summarized raw count data ",
                         "for 7706 tumor samples, represented as an ", 
                         "ExpressionSet. R data representation derived from ",
                         "GEO accession GSE62944."),
    BiocVersion = "3.2,3.3",
    Genome = "hg19",
    SourceType = "tar.gz",
    SourceUrl = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944",
    SourceVersion = "Jan 28 2015",
    Species = "Homo sapiens",
    TaxonomyId = 9606,
    Coordinate_1_based = TRUE,
    DataProvider = "GEO",
    Maintainer = "Bioconductor Package Maintainer <maintainer@bioconductor.org>",
    RDataClass = "ExpressionSet",
    DispatchClass = "ExpressionSet",
    Tags = "TCGA,RNA-seq,Expression,Count",
    ResourceName = paste0("GSE62944_GSM1536837_TCGA_20.Illumina.",
                          "tumor_Rsubread_FeatureCounts.ExpressionSet.Rda")
)

write.csv(meta, file="metadata.csv", row.names=FALSE)
