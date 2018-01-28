source("read_mixcr.R")

dt.meta = fread("metadata.txt")

dt.clones = list()

for (i in 1:nrow(dt.meta)) {
  dt.clones = with(dt.meta, 
                   c(dt.clones, 
                     list(
                       read_mixcr(sample_name[i], path[i], mixcr_v[i], replica[i], compressed[i])
                     )
                   )
  )
}

dt.clones = rbindlist(dt.clones) %>%
  merge(dt.meta) %>%
  mutate(freq = ifelse(vaccination == "yf", freq/2, freq)) %>%
  downsampling

#final <- data.frame()
#for (i in 1:10){
#  .f <- downsampling(dt.clones) %>%
#    mutate(downsampling = i)
#  final = rbind(final, .f)
#}

#dt.clones <- final

save(dt.clones, file = "mixcr_processed_downsampled.Rda")