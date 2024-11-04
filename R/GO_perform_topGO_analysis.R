#' Perform GO enrichment analysis using the topGO package
#'
#' This function performs a Gene Ontology (GO) enrichment analysis using the `topGO` package.
#' It tests for GO term enrichment in a list of genes compared to a background set of genes.
#'
#' @param gene_go_table A `data.table` containing two columns:
#'   \itemize{
#'     \item `GENE`: Gene identifiers.
#'     \item `GO`: Corresponding Gene Ontology (GO) terms.
#'   }
#' @param gene_list A character vector of gene names (in the same format as the `GENE` column of `gene_go_table`)representing the genes of interest (e.g., differentially expressed genes).
#' @param ontology A character string specifying the ontology to use. Options are "MF" (Molecular Function), "BP" (Biological Process), or "CC" (Cellular Component). Default is ALL... c("MF","BP","CC").
#' @param nodeSize Minimum number of genes per GO term
#' @param topNodes Number of nodes to report for fishertest (ranked by p-value)
#'
#' @return A data.frame containing the top enriched GO terms with both raw and adjusted p-values(Benjamini-Hochberg corrected), for both the classic and elim algorithms.
#'
#' @examples
#' # Assuming gene_go_table is a data.table of gene-to-GO mappings and gene_list is a vector of genes of interest:
#' # result <- perform_topGO_analysis(gene_go_table, gene_list, ontology = "BP")
#'
#' @import topGO
#' @export
GO_perform_topGO_analysis <- function(gene_go_table, gene_list, ontology = c("MF","BP","CC"), nodeSize=3, topNodes=100) {

  # Create the background gene universe
  gene_universe <- unique(gene_go_table$GENE)  # Background gene set

  # Create the named vector for gene status (1 for interesting genes, 0 for others)
  gene_status <- factor(as.integer(gene_universe %in% gene_list))  # 1 if gene in gene_list, 0 otherwise
  names(gene_status) <- gene_universe

  # Create the gene-to-GO mapping
  gene2GO <- split(gene_go_table$GO, gene_go_table$GENE)

  GO_results<-lapply(ontology, function(Ont){
    # Initialize topGOdata object
    GOdata <- new("topGOdata",
                  ontology = Ont,             # User-defined ontology (BP, MF, or CC)
                  allGenes = gene_status,          # The named gene vector (1 = interesting gene, 0 = background)
                  geneSel = function(p) p == 1,    # Select genes of interest (those with a 1)
                  annot = annFUN.gene2GO,          # Custom annotation function
                  gene2GO = gene2GO,               # Gene-to-GO mapping
                  nodeSize = nodeSize)             # Minimum number of genes per GO term

    # Perform Fisher’s exact test with the "classic" algorithm
    resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

    # Perform the "elim" algorithm (accounts for the GO hierarchy)
    resultElim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

    # Get p-values from the "classic" and "elim" tests
    classic_pvalues <- score(resultClassic)
    elim_pvalues <- score(resultElim)

    # Adjust p-values using the Benjamini-Hochberg (BH) method
    classic_pvalues_adj <- p.adjust(classic_pvalues, method = "BH")
    elim_pvalues_adj <- p.adjust(elim_pvalues, method = "BH")

    # Generate results table with adjusted p-values
    allRes <- GenTable(GOdata,
                       classicFisher = resultClassic,
                       elimFisher = resultElim,
                       topNodes = topNodes)

    # Add adjusted p-values to the results table
    allRes$classicFisher_adj <- classic_pvalues_adj[allRes$GO.ID]
    allRes$elimFisher_adj <- elim_pvalues_adj[allRes$GO.ID]
    allRes$Enrich=allRes$Significant/allRes$Expected
    allRes<-data.table(allRes)
    return(list(GOdata=GOdata, GO_results=allRes))
  })
  names(GO_results)<-ontology
  # Return the results
  return(GO_results)
}
