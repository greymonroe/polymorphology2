#' Nonsynonymous to Synonymous Mutation Expectation
#'
#' This function calculates the expected ratio of nonsynonymous to synonymous mutations
#' based on a given mutation dataset and a coding sequence.
#'
#' @param muts A data.table with REF and ALT columns indicating the reference and alternate nucleotides.
#' @param CDS A coding sequence fasta object read by seqinr.
#' @param reps Number of repetitions for the simulation.
#' @param n Number of samples per repetition.
#'
#' @return A list containing the results for each repetition and the overall ratio of nonsynonymous to synonymous mutations.
#' @export
#'
#' @examples
#' muts <- fread("~/Dropbox/Research/rice mutation paper/data/Strelka2_mutations.csv")[trt=="MSH6"]
#' CDS <- read.fasta("~/Documents/TAIR10_cds_20110103_representative_gene_model_updated.txt")
#' result <- nonsyn_syn_exp(muts, CDS, 10, 1000)
#'
nonsyn_syn_exp <- function(muts, CDS, reps, n) {

  nucleotides <- c("A", "T", "C", "G")


  amino_acids <- c("G", "A", "V", "L", "I", "P", "M", "F", "Y", "W", "S", "T", "C", "N", "Q", "K", "R", "H", "D", "E")
  properties <- c(rep("Non-polar\naliphatic", 7), rep("Aro-\nmatic", 3), rep("Polar\nuncharged", 5), rep("Pos", 3), rep("Neg", 2))

  properties_table <- data.frame(AA = amino_acids, Property = properties)

  # Create an empty data frame to store the mutations
  mutations <- data.frame(REF = character(0), ALT = character(0))

  # Loop through each nucleotide and find its possible mutations
  for (ref in nucleotides) {
    alts <- nucleotides[nucleotides != ref]
    temp_df <- data.frame(REF = ref, ALT = alts)
    mutations <- rbind(mutations, temp_df)
  }

  muts<-data.table(prop.table(table(REF=muts$REF, ALT=muts$ALT)))[REF!=ALT]
  muts$PROB=muts$N

  mutations<-merge(muts, mutations)
  unique_refs <- unique(mutations$REF)
  ref_probs <- sapply(unique_refs, function(ref) sum(mutations[REF == ref, PROB]))

  results <- list()

  all_results <- list(simmuts=list(), reps=list(), ratio=numeric())

  for (rep in 1:reps) {
    message(rep)
    results <- list()

    for (i in 1:n) {
      # 1. Randomly select a DNA sequence from CDS
      random_seq_name <- sample(names(CDS), 1)
      random_seq <- toupper(CDS[[random_seq_name]])

      # 2. Randomly select a REF based on its total probability

      selected_ref <- sample(unique_refs, 1, prob = ref_probs)

      # 3. Find a random base pair in random_seq with the selected REF
      positions <- which(random_seq == selected_ref)
      if(length(positions) == 0) next # Skip if the selected REF is not in the sequence
      pos <- sample(positions, 1)

      # 4. Choose a mutation based on the selected base pair and its relative probabilities
      possible_mutations <- mutations[mutations$REF == selected_ref, ]
      mutated_base <- sample(possible_mutations$ALT, 1, prob = possible_mutations$PROB)

      # 5. Introduce the mutation
      mutated_sequence <- random_seq
      mutated_sequence[pos] <- tolower(mutated_base)

      # 6. Translate both sequences
      original_protein <- seqinr::translate(random_seq)
      mutated_protein <- seqinr::translate(mutated_sequence)

      # 7. Compare the two sequences
      non_syn <- paste0(original_protein, collapse="") != paste0(mutated_protein, collapse="")

      # Break random_seq into codons
      codons <- sapply(seq(1, length(random_seq), by=3), function(x) paste0(random_seq[x:(x+2)], collapse=""))

      # Identify which codon the selected_ref belongs to
      codon_index <- ceiling(pos / 3)
      ref_codon <- codons[codon_index]

      # Identify the amino acid of ref_codon
      ref_AA <- seqinr::translate(unlist(strsplit(ref_codon, split="")))

      # Identify which position (1, 2, or 3) the mutation site belongs to
      ref_codon_pos <- pos - (3 * (codon_index - 1))

      # Identify the new codon after mutation
      mut_codon <- ref_codon
      substr(mut_codon, ref_codon_pos, ref_codon_pos) <- mutated_base

      # Identify the amino acid of mut_codon
      mut_AA <- seqinr::translate(unlist(strsplit(mut_codon, split="")))


      REF_AA_Property<-properties_table$Property[match(ref_AA, properties_table$AA)]
      mut_AA_Property<-properties_table$Property[match(mut_AA, properties_table$AA)]


      # Add these as columns to rep_results
      rep_results <- data.table(
        mut = paste(selected_ref, mutated_base, sep=">"),
        seq_name = random_seq_name,
        selected_ref,
        mutated_base,
        ref_codon,
        ref_AA,
        REF_AA_Property,
        ref_codon_pos,
        mut_codon,
        mut_AA,
        mut_AA_Property,
        non_syn,
        Property_mut=REF_AA_Property!=mut_AA_Property
      )

      # Store results
      results[[i]] <- rep_results

      # Update the progress bar
    }

    final_results <- rbindlist(results)

    ns_s <- data.table(table(mut=final_results$mut, ns=final_results$non_syn))
    wide_ns_s <- dcast(ns_s, mut ~ ns, value.var = "N")
    wide_ns_s$Total<-wide_ns_s$`TRUE`+wide_ns_s$`FALSE`
    wide_ns_s$NS_S = wide_ns_s$`TRUE`/wide_ns_s$`FALSE`
    wide_ns_s$rep=rep

    ns_s_rate <- table(final_results$non_syn)[2]/table(final_results$non_syn)[1]

    all_results$simmuts[[rep]] <-  final_results
    all_results$reps[[rep]] <-  wide_ns_s
    all_results$ratio[rep] = ns_s_rate
  }

  all_results
}
