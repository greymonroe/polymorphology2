#' Extract genes of interest associated with a specific GO term
#'
#' This function extracts the genes from your `gene_list` that are associated with a given GO term (GO.ID)
#' from the `topGOdata` object generated during a GO enrichment analysis.
#'
#' @param GOdata A `topGOdata` object created during the GO enrichment analysis.
#' @param go_term_of_interest A character string specifying the GO term (GO.ID) of interest.
#' @param gene_list A character vector of gene names that represents your genes of interest.
#'
#' @return A character vector containing the genes from `gene_list` associated with the specified GO term.
#'
#' @examples
#' # Assuming `GOdata` is your topGOdata object and `gene_list` is your list of genes:
#' # genes_of_interest <- extract_genes_for_go_term(GOdata, "GO:0030941", gene_list)
#'
#' @export
extract_genes_for_go_term <- function(GOdata, go_term_of_interest, gene_list) {

  # Step 2: Extract the genes associated with this GO.ID from the topGOdata object
  genes_in_go_term <- genesInTerm(GOdata, go_term_of_interest)

  # Step 3: Filter the genes to include only those that are in your gene_list
  genes_of_interest_in_go <- intersect(genes_in_go_term[[go_term_of_interest]], gene_list)

  # Return the filtered list of genes
  return(genes_of_interest_in_go)
}

# Example usage:
# genes_of_interest <- extract_genes_for_go_term(GOdata, "GO:0030941", gene_list)
# print(genes_of_interest)
