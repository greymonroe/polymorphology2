#' strelka2_filter
#'
#' This function reads in results from strelka2 VCF files and applies various filters to identify de novo mutations.
#'
#' @param directory A string specifying the path to the directory containing the VCF files.
#' @param output A string specifying the path to the output directory.
#'
#' @return The function writes two CSV files in the output directory: "strelka2_all_calls.csv" and "strelka2_filtered_calls.csv".
#'
#' @importFrom data.table rbindlist fread fwrite
#' @export
strelka2_filter <- function(directory, output) {
  # load necessary libraries
  require(data.table)

  bp <- c("A","C","G","T") # base pairs

  normals <- list.files(directory, full.names = FALSE)
  tumors <- normals

  results <- rbindlist(lapply(normals, function(norm){
    message("Processing normal sample: ", norm)
    mutations <- rbindlist(lapply(tumors, function(tumor){
      if (norm != tumor){
        message("\tComparing with tumor sample: ", tumor)
        vcf_file<-paste0(directory, norm, "/", tumor, "/results//variants/somatic.snvs.vcf.gz")
        if(file.exists(vcf_file)){
        calls <- read.VCF(paste0(directory, norm, "/", tumor, "/results//variants/somatic.snvs.vcf.gz"), caller = "strelka2")
        } else return(NULL)
        calls[, c("NORMAL_SAMPLE", "TUMOR_SAMPLE") := .(norm, tumor)]
        return(calls)
      }
    }))

    message("Processing mutations for sample: ", norm)
    mutations[, unique := paste(CHROM, POS, ALT, sep = "_")]
    mutations[, Mut_N := .N, by = .(unique)]
    return(mutations)
  }))

  message("Finalizing the results...")

  # Check if the output file already exists
  if (file.exists(paste0(output, "/strelka2_all_calls.csv"))) {
    warning("The file already exists. Overwriting...")
  }

  results[, unique2 := paste(unique, TUMOR_SAMPLE, sep = "_")]
  results[, results_N := .N, by = .(unique2)]

  message("Writing all calls to the output file...")
  fwrite(results, paste0(output, "/strelka2_all_calls.csv"))



  message("Applying filters...")
  PASS <- results[Mut_N == 1 & FILTER == "PASS"]

  PASS[, N := .N, by = .(unique2)]

  PASS <- PASS[N == length(normals) - 1]

  # calculate depths of variants called


  PASS$tumor_ref_depth<-as.numeric(apply(PASS, 1, function(x) unlist(strsplit(unlist(strsplit(x["TUMOR"], split = ":"))[4+which(bp==x["REF"])],split=","))[1]))
  PASS$tumor_alt_depth<-as.numeric(apply(PASS, 1, function(x) unlist(strsplit(unlist(strsplit(x["TUMOR"], split = ":"))[4+which(bp==x["ALT"])],split=","))[1]))
  PASS[, tumor_depth := tumor_ref_depth + tumor_alt_depth]
  PASS[, depth_pct := tumor_alt_depth / tumor_depth]

  PASS <- unique(PASS[, .(CHROM, POS, REF, ALT, TUMOR_SAMPLE, unique, tumor_alt_depth, tumor_depth, depth_pct, N, EVS = mean(EVS)), by=unique2])

  message("Calculating Chi-square tests and distances...")
  PASS[, chi := apply(.SD, 1, function(x) chisq.test(c(as.numeric(x["tumor_alt_depth"]), as.numeric(x["tumor_depth"]) - as.numeric(x["tumor_alt_depth"])))$p.value)]

  PASS[, dist := sapply(1:.N, function(i) {
    x <- .SD[i]
    chr <- x[["CHROM"]]
    POS <- as.numeric(x[["POS"]])
    f <- x[["TUMOR_SAMPLE"]]
    dist <- PASS[CHROM == chr & TUMOR_SAMPLE == f]$POS - POS
    min(abs(dist[dist != 0]))
  })]

  message("Classifying mutation types...")
  PASS[, type := ifelse(depth_pct > 0.5 & chi < 0.05, "homozygous", "somatic")]
  PASS[type == "somatic" & chi > 0.001, type := "heterozygous"]

  PASS <- PASS[!duplicated(unique2)]

  # Check if the output file already exists
  if (file.exists(paste0(output, "/strelka2_filtered_calls.csv"))) {
    warning("The file already exists. Overwriting...")
  }

  message("Writing filtered calls to the output file...")
  fwrite(PASS, paste0(output, "/strelka2_filtered_calls.csv"))

  message("Completed the strelka2_filter function.")
  return(PASS)
}

