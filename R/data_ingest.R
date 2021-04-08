#' Download GDSC Data
#'
#' Function to download the latest version of RNA expression, drug response,
#' and metadata from the Genomics of Drug Sensitivity in Cancer database.
#' link [here](https://www.cancerrxgene.org)
#'
#' @inheritParams clean_gdsc
#'
#'
#' @export

get_gdsc <- function(cache = NULL) {
  if(is.null(cache)) {
    stop("Please specify a filepath specifying where to store data")
  }

  #get expression data
  exp_df <- get_exp(cache = cache)

  #get drug data
  drug_df <- get_drug()

  #get metadata
  meta_df <- get_meta(cache = cache)

  dna_df <- get_dna()

  return(list(exp_df, drug_df, meta_df, dna_df))

}

#' helper function to download expression data specifically
#'
#' @inheritParams clean_gdsc
#'

get_exp <- function(cache) {

  # check if file has already been downloaded at the provided cache
  if(!file.exists(paste0(cache, "/Cell_line_RMA_proc_basalExp.txt"))) {
    download.file("https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/preprocessed/Cell_line_RMA_proc_basalExp.txt.zip",
                  destfile = paste0(cache, "/gdsc_exp.zip"))
    unzip(zipfile = paste0(cache, "/gdsc_exp.zip"), exdir = cache)
    file.remove(paste0(cache, "/gdsc_exp.zip")) #remove the zombie zip file
  } else {
    message("using cached gdsc data")
  }

  exp_df <- readr::read_tsv(paste0(cache, "/Cell_line_RMA_proc_basalExp.txt"))
  return(exp_df)

}

#' helper function to get drug response data
#'
#'
#'
#'


get_drug <- function() {

  #download the data
  url1<- "https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS4A.xlsx"
  httr::GET(url1, httr::write_disk(tf <- tempfile(fileext = ".xlsx")))

  #read into R
  df <- readxl::read_excel(tf, skip = 5)

  #rename the first 2 columns because of an excel generated parsing error
  colnames(df)[1:2] <- c("cosmic_identifier", "sample_name")
  return(df)

}

#' helper function to get the cell line metadata
#'
#' @inheritParams clean_gdsc
#'

get_meta <- function(cache) {

  if(!file.exists(paste0(cache, "/Cell_Lines_Details.xlsx"))) {
    curl::curl_download("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/Cell_Lines_Details.xlsx",
                        destfile = paste0(cache, "/Cell_Lines_Details.xlsx"))
  } else {
    message("using cached metadata")
  }
  df <- readxl::read_excel(paste0(cache, "/Cell_Lines_Details.xlsx"))
  return(df)
}

#' helper function to get DNA seq data
#'
#'

get_dna <- function() {
  #download the data
  url1<- "https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS2C.xlsx"
  httr::GET(url1, httr::write_disk(tf <- tempfile(fileext = ".xlsx")))

  #read into R
  df <- suppressWarnings(readxl::read_excel(tf, skip = 18)) # do it quietly
  return(df)
}

#' clean the GDSC data
#'
#'
#' @param cache filepath specifying where to store store downloaded/ processed data
#'
#' @export

clean_gdsc <- function(cache = NULL) {
  if(is.null(cache)) {
    stop("Please specify a filepath specifying where to store data")
  }

  df_list <- get_gdsc(cache = cache)

  #grab the metadata out of the list of dirty dataframes
  meta_df <- df_list[[3]]
  meta_df <- clean_meta(meta_df)

  #grab and clean the response data
  response_df <- df_list[[2]]
  response_df <- clean_response(response_df)

  #grab and clean expression data
  exp_df <- df_list[[1]]
  exp_df <- clean_expression(exp_df)

  #grab and clean dna variants data
  dna_df <- df_list[[4]]
  dna_df <- clean_dna(dna_df)

  return(list(meta_df, response_df, exp_df, dna_df))

}

#' helper function for processing the gdsc metadata
#'
#' Here we read in the meta data for GDSC. This contains all the cell lines by
#' various identifiers (we use `COSMIC_ID`) and descriptions of each line,
#' including tissue and histology of origin. We also do a quick cleaning step by
#' adjusting a variable name containing a forward slash (`/`) to avoid errors later.
#' We will also split out epithelial vs non epithelial cell lines.
#'
#' @param df unclean data frame
#'

clean_meta <- function(df) {

  df <- janitor::clean_names(df)
  levels(df$cancer_type_matching_tcga_label) <-
    c(levels(df$cancer_type_matching_tcga_label), "COAD&READ")
  df$cancer_type_matching_tcga_label[df$cancer_type_matching_tcga_label=="COAD/READ"] <- "COAD&READ"

  gdsc_interest <- c("head and neck", "oesophagus", "breast", "biliary_tract",
                     "digestive_system_other","large_intestine", "stomach",
                     "lung_NSCLC_adenocarcinoma", "lung_NSCLC_carcinoid",
                     "lung_NSCLC_large cell", "lung_NSCLC_not specified",
                     "lung_NSCLC_squamous_cell_carcinoma", "Lung_other",
                     "pancreas", "skin_other", "thyroid", "Bladder", "cervix",
                     "endometrium", "ovary", "prostate", "testis",
                     "urogenital_system_other", "uterus", "liver", "kidney")

  #Add a variable for if the cell line is epithelial in nature.
  df <-  dplyr::mutate(df,
                       epi_origin = ifelse(
                         gdsc_tissue_descriptor_2 %in% gdsc_interest,
                         TRUE,
                         FALSE
                       ))

  return(df)

}

#' Import and clean drug response data
#'
#' Here we import the drug response data from the second version of the GDSC dataset (GDSC2).
#' Each cell line can still be identified using `COSMIC_ID`.
#' We clean the data by adjusting data types and converting IC50 from natural log to log2.
#'
#' @inheritParams clean_meta
#'
#' @importFrom magrittr %>%

clean_response <- function(df) {

  #all but the first two columns (the two identifiers, we want to be numeric)
  df <- df %>%
    tidyr::pivot_longer(cols = -c(1,2), names_to = "drug", values_to = "IC50") %>%
    dplyr::mutate(IC50 = as.numeric(.data$IC50),
                  IC50 = exp(.data$IC50),
                  IC50 = log2(.data$IC50))

  z_df <- df %>%
    dplyr::group_by(drug) %>%
    dplyr::group_modify(
      ~data.frame(
        cosmic_identifier = .x$cosmic_identifier,
        sample_name = .x$sample_name,
        IC50 = .x$IC50,
        z_score = as.numeric(scale(.x$IC50))
      ), .keep = TRUE)

  return(z_df)
}

#' Function to clean expression data from gdsc
#' We'll read in the GDSC expression data, which comes from microarray experiments.
#' Details on the experimental protocols can be found
#' [here](https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html).
#'
#' @inheritParams clean_meta
#'
#' @importFrom magrittr %>%

clean_expression <- function(df) {
  df <- df %>%
    tidyr::pivot_longer(-c(1,2), names_to = "cell_line",
                        values_to = "expression",
                        names_prefix = "DATA\\.") %>%
    janitor::clean_names() %>%
    dplyr::mutate(log_expression = log2(expression))
 return(df)
}

#' Function to clean dna variant data from gdsc
#' data can be fouund [here](https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html).
#'
#' @inheritParams clean_meta
#'
#' @importFrom magrittr %>%
#'

clean_dna <- function(df) {
  df <- df %>%
    janitor::clean_names() %>%
    dplyr::rename(cell_line = cosmic_id)
  return(df)
}


