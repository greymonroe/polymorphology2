
#' Nonsynonymous to Synonymous Mutation Expectation
#'
#' This function calculates the expected ratio of nonsynonymous to synonymous mutations
#' based on a given mutation dataset and a coding sequence, accounting for base pair composition of genome. For pre-calculated bp_freq for a genome sequence use nonsyn_syn_exp4()
#'
#' @param muts A data.table with REF and ALT columns indicating the reference and alternate nucleotides.
#' @param CDS A coding sequence fasta object read by seqinr.
#' @param genome A fasta object read in by seqinr read.fasta()
#' @param reps Number of repetitions for the simulation.
#'
#' @return A list containing the results for each repetition and the overall ratio of nonsynonymous to synonymous mutations.
#' @export

nonsyn_syn_exp3 <- function(muts, genome, CDS, reps) {

  bp_freq<-rbindlist(lapply(1:length(genome), function(i){
    bp_freq<-data.table(table(REF=genome[[i]]))[REF!="n"]
    bp_freq$REF<-toupper(bp_freq$REF)
    return(bp_freq)
  }))
  bp_freq<-bp_freq[,.(N=sum(N)), by="REF"]

mutations<-data.table((table(REF=muts$REF, ALT=muts$ALT)))[REF!=ALT]
mutations<-merge(mutations, bp_freq, by="REF")
mutations$PROB=prop.table(mutations$N.x/mutations$N.y)
mutations$PROB2=prop.table(mutations$N.x)
unique_refs <- unique(mutations$REF)
ref_probs <- sapply(unique_refs, function(ref) sum(mutations[REF == ref, PROB]))
ref_probsdt<-data.table(REF=names(ref_probs), PROB=ref_probs)

CDS_dt_all<-rbindlist(lapply(sample(1:length(CDS), reps), function(prot){
  random_seq<-toupper(CDS[[prot]])
  CDS_dt<-data.table(POS=1:length(random_seq), REF=random_seq)
  CDS_dt<-merge(CDS_dt, mutations[,.(REF, ALT, PROB)], by="REF",  allow.cartesian=TRUE)[order(POS),.(MUT=.I, POS, REF, ALT, PROB=prop.table(PROB))]
  mutation<-CDS_dt[sample(CDS_dt$MUT, round(length(random_seq)*.02), prob = CDS_dt$PROB)]

  # 6. Translate both sequences
  original_protein <- seqinr::translate(random_seq)

  # 5. Introduce the mutation

  mutation$effects<-apply(mutation, 1, function(r){
    mutated_sequence <- random_seq
    pos<-as.numeric(r["POS"])
    mutated_base<-r["ALT"]
    mutated_sequence[pos] <- mutated_base
    mutated_protein <- seqinr::translate(mutated_sequence)
    non_syn <- paste0(original_protein, collapse="") != paste0(mutated_protein, collapse="")
    non_syn<-ifelse(non_syn, "Non-Synonymous","Synonymous")
  })
  return(mutation)
}))

}

