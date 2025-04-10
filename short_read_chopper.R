

library(data.table)
library(dplyr)

option_list = list(
  make_option(c("--tsv",), type="character", default=NULL, 
              help="path to tsv with simulated mRNAs", metavar="character"),
  make_option(c("--insertsize"), type="numeric", default="200,300", 
              help="library insert size range", metavar="numeric"),
  
  make_option(c("--read_length"), type="numeric", default="100", 
              help="length of each individual read", metavar="numeric"),
  
  make_option(c("--single_or_paired"), type="character", default="SE", 
              help="is this a SE or a PE library?", metavar="character"),
  
  make_option(c("--library_direction"), type="character", default="rf", 
              help="library strandness (rf,fr or unstranded)", metavar="numeric"),
  
  make_option(c("-t", "--threads"), type="numeric", default="1", 
              help="Number of threads", metavar="numeric"),
  make_option(c("-o", "--dir_out"), type="character", default=".", 
              help="dir out to create temp files", metavar="character")
)

#Parse input
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$tsv)){
  stop("path to tsv with simulated mRNAs is missing", call.=FALSE)
}

setDTthreads(threads = opt$threads)
insertsize <- opt$insertsize
file <- opt$tsv

#File importing and formating
file_importer <- function(file){
  full_reads <- fread(file, sep = '\t', header = T, stringsAsFactors = FALSE)
  full_reads$transcript_id <- as.numeric(1:nrow(full_reads))
  full_reads$sequence_length <- nchar(full_reads$full_molecule_sequence)
  return(full_reads)
}
#Chopper function
get_reads <- function(lengths, eta_val = 200, insertsize){
  deltas = log10(lengths)
  ns_minus_1 = pmax(round(lengths/eta_val/gamma(1/deltas + 1)) - 1, 0)
  xis = lapply(ns_minus_1, function(n) {diff(sort(c(runif(n), 0, 1)))})
  xis_transformed = mapply(function(x, d) {x^(1/d)}, xis, deltas, SIMPLIFY = F)
  delta_is = mapply(function(len, x_t) {round(len*x_t/sum(x_t))}, lengths, xis_transformed, SIMPLIFY = F)

# get all the start and end points of the fragments 
  starts = lapply(delta_is, function(d) {
    if (length(d) > 1) {
      c(sample(min(insertsize[1], d[1]), 1), cumsum(d[1:(length(d)-1)]))
    } else{
      sample(min(insertsize[1], d), 1)
    }
  })
  ends = lapply(delta_is, function(d) {
    if (length(d) > 1) {
      c(cumsum(d[1:(length(d)-1)]), sum(d)-sample(min(insertsize[1], sum(d) - d[length(d)]), 1))
    } else{
      d
    }
  })
  fragments = data.frame(transcript = rep(1:length(deltas), unlist(lapply(delta_is, length))),
                       start = unlist(starts),
                       end = unlist(ends))
  fragments$length = fragments$end - fragments$start

  # Filter fragments by length and return
  fragments = fragments[fragments$length >= insertsize[1] & fragments$length <= insertsize[2],]
  reads_list <- fragments[c('transcript', 'start','end')]
  colnames(reads_list) <- c('transcript_id', 'read_start','read_end')
  reads_list$transcript_id <- as.numeric(reads_list$transcript_id)
  return(reads_list)
}

random_string_gen <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

full_reads <- file_importer(file)
reads_list <- lapply(full_reads$sequence_length, function(x) get_reads(x, insertsize = insertsize))
reads_list <- dplyr::bind_rows(reads_list,.id = 'transcript_id')
full_reads <- merge(full_reads,reads_list,by='transcript_id')

#Add unique read identifier
full_reads$molecule_id <- paste0(full_reads$molecule_id, "_", replicate(nrow(full_reads), paste0(sample(letters, 5, replace = TRUE), collapse = "")))

#Export
temp_dir <- paste(opt$dir_out,'/temp/mRNAs_with_fragments',sep = '')
system(paste("mkdir ", temp_dir,sep = ''))
path_for_file <- paste(temp_dir,'/',sub(".*/([^/_]+)_.*", "\\1", file),'_fragments.tsv',sep = '')
fwrite(full_reads,path_for_file,sep = '\t',col.names = T,row.names = F,quote = F)

#Run the python chopper
path_to_this_script <- rstudioapi::getSourceEditorContext()$path


paste('python ',gsub("/[^/]+$", "", path_to_this_script),' ')




path_for_file
opt$read_length
opt$single_or_paired
opt$library_direction






data.frame(chr=c(2,19,17),
           start=c(71266561,12823951,19998865),
           end=c(71455069,12894952,20329026),
           name=c('ZNF638','MAST1','SPECC1'),
           score=c('.','.','.'),
           strand=c('+','+','+')) %>% write.table('/pi/athma.pai-umw/data/ezequiel/TSS-TES/PAS_deletions/NextSeq_3_17_25/bed_file_ref.bed',sep = '\t',row.names = F,col.names = F,quote = F)






























