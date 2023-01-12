##### This script downloads all fastqfiles from a serie in NCBI GEO with RNA-Seq data #####
##### and accesion to SRA. The function links NCBI GEO and ENA browser, get the links #####
##### and downloads fastq files. The input should be a Serie from GEO like GSE167923  #####

##### Example #####

# run in the terminal:

# Rscript getGEOFasqt.R --serie GSE167923

library(optparse)

option_list <- list(make_option("--serie", type = "character", default = "GSE167923"))
option_list <- parse_args(OptionParser(option_list=option_list))

downloadFastqFiles <- function(srx){
  require(GEOquery)
  require(glue)
  require(vroom)
  serie_info <- getGEO(serie,GSEMatrix = F)
  serie_info <- serie_info@gsms
  dir.create(serie, showWarnings = F)
  lapply(names(serie_info), function(gsm){
    info <- getGEO(gsm)@header$relation
    info_sra <- info[grep("SRA:",info)]
    srx <- unlist(strsplit(info_sra,"term="))[2]
    baseurl = glue("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={srx}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true&limit=0")
    download.file(baseurl,destfile = "info_table.tsv")
    info_table <- vroom("info_table.tsv")
    fastq_files <- unlist(strsplit(info_table$fastq_ftp,";"))
    basenames <- basename(fastq_files)
    basenames <- gsub(info_table$run_accession, gsm, basenames)
    basenames <- paste0(serie,"/",basenames)
    mapply(download.file, url = fastq_files, destfile = basenames)
  })
}

downloadFastqFiles(option_list$serie)
