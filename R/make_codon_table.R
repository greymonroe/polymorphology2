#' Create a Standard Genetic Codon Table
#'
#' This function generates a data table representing the standard genetic code.
#' The table maps codons to their corresponding amino acids and also provides
#' the count of each nucleotide (A, T, C, G) in the codon.
#'
#' @return A data.table containing the following columns:
#' \itemize{
#'   \item \strong{AA}: Amino acid corresponding to the codon.
#'   \item \strong{codon}: The codon sequence.
#'   \item \strong{A}: Count of adenine in the codon.
#'   \item \strong{T}: Count of thymine in the codon.
#'   \item \strong{C}: Count of cytosine in the codon.
#'   \item \strong{G}: Count of guanine in the codon.
#' }
#'
#' @examples
#' codon_tbl <- make_codon_table()
#' head(codon_tbl)
#'
#' @export

make_codon_table<-function(){
# Define the standard genetic code
codons <- c(
  "TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
  "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG",
  "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG",
  "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG",
  "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG",
  "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG",
  "TGT", "TGC", "TGA", "TGG", "CGT", "CGC", "CGA", "CGG",
  "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG"
)

amino_acids <- c(
  "F", "F", "L", "L", "L", "L", "L", "L",
  "I", "I", "I", "M", "V", "V", "V", "V",
  "S", "S", "S", "S", "P", "P", "P", "P",
  "T", "T", "T", "T", "A", "A", "A", "A",
  "Y", "Y", "*", "*", "H", "H", "Q", "Q",
  "N", "N", "K", "K", "D", "D", "E", "E",
  "C", "C", "*", "W", "R", "R", "R", "R",
  "S", "S", "R", "R", "G", "G", "G", "G"
)

# Create the data.table
codon_table <- data.table(
  AA = amino_acids,
  codon = codons,
  A = nchar(gsub("[^A]", "", codons)),
  T = nchar(gsub("[^T]", "", codons)),
  C = nchar(gsub("[^C]", "", codons)),
  G = nchar(gsub("[^G]", "", codons))
)

return(codon_table)
}
