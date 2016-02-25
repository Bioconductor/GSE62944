### =========================================================================
### GSE62944 ExpressionSet 
### -------------------------------------------------------------------------
###

##  Download files
acc <- "GSE62944"
rootDir <- getwd() 
outputPath <- file.path(rootDir, tempfile())
getGEOSuppFiles(acc, baseDir=file.path(rootDir, "geo"))

##  Parse
fl <- "GSE62944_TCGA_20_420_Clinical_Variables_7706_Samples.txt.gz"
m <- scan(fl, what=character(), sep="\t", quote="")
m <- matrix(m, 7707)
dimnames(m) <- list(m[,1], m[1,])
df <- as.data.frame(m[-1, -1])

fl <- "GSE62944_TCGA_20_CancerType_Samples.txt.gz"
ct <- read.delim(fl, header=FALSE,
                 colClasses=c("character", "factor"),
                 col.names=c("sample", "type"))
idx <- match(rownames(df), ct$sample)
stopifnot(!anyNA(idx))
df$CancerType <- ct$type[idx]
clinvar <- df

## Create ExpressionSet 
GSE62944ToExpressionSet <- function(ahm) {
message("counts")
fl <- "GSM1536837_TCGA_20.Illumina.tumor_Rsubread_FeatureCounts.txt.gz"
if (!file.exists(fl))
    untar("GSE62944_RAW.tar", fl)
m <- scan(fl, what=character(), sep="\t", quote="")
m <- matrix(m, 7707)
dimnames(m) <- list(m[,1], m[1,])
m <- t(m[-1, -1])
mode(m) <- "integer"
counts <-  m

    adf <- Biobase::AnnotatedDataFrame(clinvar)
    eset <- Biobase::ExpressionSet(counts, adf)

    save(eset, file=outputPath)
    outputFile(ahm)
}

makeExperimentHubMetadata(
    Title = paste0("RNA-Sequencing and clinical data for 7706 tumor ",
                   "samples from The Cancer Genome Atlas"),
    Description = paste0("TCGA RNA-seq Rsubread-summarized raw count data ",
                         "for 7706 tumor samples, represented as an ", 
                         "ExpressionSet. R data representation derived from ",
                         " GEO accession GSE62944."),
    BiocVersion = c("3.2", "3.3"),
    Genome = "hg19",
    SourceType = "tar.gz",
    SourceUrl = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944",
    SourceVersion = "Jan 28 2015",
    Species = "Homo sapiens",
    TaxonomyId = 9606L,
    RDataPath = paste0("geo/GSE62944/GSE62944_GSM1536837_TCGA_20.",
                       "Illumina.tumor_Rsubread_FeatureCounts.",
                       "ExpressionSet.Rda"),
    Coordinate_1_based = TRUE,
    DataProvider = "GEO",
    Maintainer = "Bioconductor Maintainer <maintainer@bioconductor.org>",
    RDataClass = "ExpressionSet",
    Tags = c("TCGA", "RNA-seq", "Expression", "Count")
)

## This function or something similar would go in ExperimentHub
makeExperimentHubMetadata <- function(Title, Description, 
    BiocVersion = biocVersion(), Genome, SourceType, SourceUrl, 
    SourceVersion, Species, TaxonomyID, Coordinate_1_based, 
    DataProvider, RDataClass, RDataDate, Tags)
{
    ## check tags against biocviews
    ## check RDataClass against object 
    Maintainer = ## taken from Maintainer field
    RdataPath = ## ?

    data.frame(Title, Description, BiocVersion, Genome, SourceType,
               SourceUrl, SourceVersion, Species, TaxonomyID,
               RDataPath, Coordinate_1_based, DataProvider, 
               Maintainer, RDataClass, RdataDateAdded = Sys.time(), 
               Location_Prefix = "http://s3.amazonaws.com/experimenthub/",
               Tags, Recipe = NULL, DispatchClass = NULL)
}
