library(data.table)
library(dplyr)
library(stringr)
library(reshape2)

select = dplyr::select
summarise = dplyr::summarise

read_mixcr_old = function(cmd) {
  fread(cmd) %>%
    filter(nchar(`AA. Seq. CDR1`) > 0, nchar(`AA. Seq. CDR2`) > 0, nchar(`AA. Seq. CDR3`) >= 5,
           !grepl("[*_]", paste(`AA. Seq. CDR1`, `AA. Seq. CDR2`, `AA. Seq. CDR3`, `AA. Seq. FR1`, `AA. Seq. FR2`, `AA. Seq. FR3`))) %>%
    group_by(v = str_split_fixed(`All V hits`, fixed("*"), 2)[,1],
             j = str_split_fixed(`All J hits`, fixed("*"), 2)[,1],
             isotype = str_split_fixed(`All C hits`, fixed("*"), 2)[,1],
             cdr1aa = `AA. Seq. CDR1`,
             cdr2aa = `AA. Seq. CDR2`,
             cdr3aa = `AA. Seq. CDR3`,
             fr1aa = `AA. Seq. FR1`,
             fr2aa = `AA. Seq. FR2`,
             fr3aa = `AA. Seq. FR3`,
             fr4aa = `AA. Seq. FR4`,
             contignt = `Clonal sequence(s)`,
             cdr3nt = `N. Seq. CDR3`,
             vAlign = `All V alignments`,
             jAlign = `All J alignments`,
             refs = `Ref. points`) %>%
    summarise(count = sum(`Clone count`), freq = sum(`Clone fraction`)) %>%
    transform_refpoints %>%
    transform_alignment %>%
    merge_aa_seqs
}

read_mixcr_new = function(cmd) {
  fread(cmd) %>%
    filter(nchar(aaSeqCDR1) > 0, nchar(aaSeqCDR2) > 0, nchar(aaSeqCDR3) >= 5,
           !grepl("[*_]", paste(aaSeqCDR1, aaSeqCDR2, aaSeqCDR3, aaSeqFR1, aaSeqFR2, aaSeqFR3))) %>%
    group_by(v = str_split_fixed(allVHitsWithScore, fixed("*"), 2)[,1],
             j = str_split_fixed(allJHitsWithScore, fixed("*"), 2)[,1],
             isotype = str_split_fixed(allCHitsWithScore, fixed("*"), 2)[,1],
             cdr1aa = aaSeqCDR1,
             cdr2aa = aaSeqCDR2,
             cdr3aa = aaSeqCDR3,
             fr1aa = aaSeqFR1,
             fr2aa = aaSeqFR2,
             fr3aa = aaSeqFR3,
             fr4aa = aaSeqFR4,
             contignt = clonalSequence,
             cdr3nt = nSeqCDR3,
             vAlign = allVAlignments,
             jAlign = allJAlignments,
             refs = refPoints) %>%
    summarise(count = sum(cloneCount), freq = sum(cloneFraction)) %>%
    transform_refpoints %>%
    transform_alignment %>%
    merge_aa_seqs
}

transform_refpoints = function(data) {
  tmp = as.data.table(str_split_fixed(data$refs, ":", n = 22)[,c(10,12,13,16,17,6,7,8,9,19)])
  
  colnames(tmp) = c("cdr3Start", "vEnd", "dStart", "dEnd", "jStart",
                    "cdr1Start", "fr2Start", "cdr2Start", "fr3Start", "fr4Start")
  
  tmp$cdr3Start = as.integer(tmp$cdr3Start)
  tmp$vEnd = as.integer(tmp$vEnd)
  tmp$dStart = as.integer(tmp$dStart)
  tmp$dEnd = as.integer(tmp$dEnd)
  tmp$jStart = as.integer(tmp$jStart)
  tmp$cdr1Start = as.integer(tmp$cdr1Start)
  tmp$fr2Start = as.integer(tmp$fr2Start)
  tmp$cdr2Start = as.integer(tmp$cdr2Start)
  tmp$fr3Start = as.integer(tmp$fr3Start)
  tmp$fr4Start = as.integer(tmp$fr4Start)
  
  data = cbind(as.data.table(data), tmp)
  data$refs = NULL
  data
}

transform_alignment = function(data){
  data$vMutations <- sapply(data$vAlign, function(x) str_match(x, '\\d+\\|\\d+\\|\\d+\\|\\d+\\|\\d+\\|([\\w\\d>:]+)')[1,2])
  data$jMutations <- sapply(data$jAlign, function(x) str_match(x, '\\d+\\|\\d+\\|\\d+\\|\\d+\\|\\d+\\|([\\w\\d>:]+)')[1,2])
  
  data$vAlignStart <- sapply(data$vAlign, function(x) str_match(x, "(\\d+)|")[1])
  data$jAlignStart <- sapply(data$jAlign, function(x) str_match(x, "(\\d+)|")[1])
  
  data$vAlign <- NULL
  data$jAlign <- NULL
  data
}

merge_aa_seqs = function(data){
  data$aa.seq <- data %>% select(fr1aa, cdr1aa, fr2aa, cdr2aa, fr3aa, cdr3aa, fr4aa) %>%
    apply(1, function(i) paste0(i, collapse=""))
  
  data$cdr1aa <- NULL
  data$cdr2aa <- NULL
  data$fr1aa <- NULL
  data$fr2aa <- NULL
  data$fr3aa <- NULL
  data$fr4aa <- NULL
  data
}

read_mixcr = function(sample, path, mixcr_v, replica, compressed) {
  format = ifelse(compressed, ".txt.gz", ".txt")
  tool = ifelse(compressed, "zcat ", "cat ")
  
  if (mixcr_v == "new") {
    res = read_mixcr_new(paste0(tool, path, "/", sample, format))
  } else {
    res = read_mixcr_old(paste0(tool, path, "/", sample, format))
  }
  res %>% mutate(sample_name = sample, replica = replica)
}


downsampling = function(data){
  df.new <- data.frame()
  
  df.size = data %>% group_by(sample_name) %>%
    summarise(n = sum(count))
  sample.size = min(df.size$n)
  
  for (s in unique(data$sample_name)){
    .df <- filter(data, sample_name == s)
    .df$clone = 1:nrow(.df)
    row.names(.df) <- .df$clone
    .df <- .df[rep(row.names(.df), .df$count), ] %>%
      sample_n(size = sample.size) %>%
      group_by_at(vars(-count, -freq)) %>%
      summarise(count = n()) %>%
      mutate(freq = count/nrow(.df))
    
    df.new <- rbind(df.new, as.data.frame(.df))
  }
  
  df.new
}