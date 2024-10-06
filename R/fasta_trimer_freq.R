
fasta_trimer_freq<-function(fasta, stop = NULL){
  fasta_windows<-rbindlist(lapply(1:length(fasta), function(i){
    data.table(CHROM=names(fasta)[i], START=1, STOP=length(fasta[[i]]), ID=i)
    data.table(CHROM=names(fasta)[i], START=1, STOP=ifelse(is.null(stop), length(fasta[[i]]), stop), ID=i)

  }))[!is.na(CHROM)]
  trimer_freq<-Nmerfrequency(features = fasta_windows, fasta=fasta, Nmer = 3, mode="counts")
  trimer_freq<-rbindlist(lapply(1:ncol(trimer_freq), function(i){
    data.table(TRI=names(trimer_freq)[i], N=sum(trimer_freq[[i]]))
  }))

  trimer_freq<-trimer_freq[!grepl("ID|N|V", TRI)]
}

