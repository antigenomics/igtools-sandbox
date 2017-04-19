# AGING DATASET

load('sp.Rda')
aging.shm <- df %>% filter(type == 'RNA') %>%
  unnest(clonotype.id = clonotype.ids) %>% 
  mutate(rate = rate/total.clonotypes, total.clonotypes = 1)

old <- c("Abdulain", "Ilgen", "Mamaev", "Smirnov", "Vlasov")
young <- c("Antipyat", "Epifancev", "Hadjibekov", "Koshkin", "Kovalchuk")
samples <- data.frame(sample = c(old,young), proj = c(rep('old', length(old)), rep('young', length(young))))

clones <- data.frame()

for (i in 1:nrow(samples)){
  .clones <- read.table(paste0('../yf/data/yf_', samples$proj[i], '_RNA/', samples$sample[i], '.clones.txt'), header=T, sep="\t")
  .clones$proj <- samples$proj[i]
  .clones$sample <- samples$sample[i]
  .clones$clonotype.id <- (1:nrow(.clones))-1
  clones <- rbind(clones, .clones)
}

aging.shm.1 <- merge(aging.shm, dplyr::select(clones, clonotype.id, proj, sample, freq, contignt), 
                     by=c('clonotype.id', 'proj', 'sample'))
aging.shm.1$contignt <- as.character(aging.shm.1$contignt)
aging.shm.1 <- mutate(aging.shm.1, len = nchar(contignt),
                      context = str_sub(contignt, pmax(0,pos.nt-3), pmin(pos.nt+3, len)),
                      cells = 'P')


# RAJI DATASET

melt_mutations <- function(df){
  df %>% mutate(clonotype.id = (1:nrow(df))-1) %>%
    dplyr::select(clonotype.id, freq, starts_with('mutations.nt.'), contignt) %>% 
    mutate_each(funs(as.character(.)), starts_with('mutations.nt.')) %>%
    melt(id.vars = c('clonotype.id', 'freq','contignt')) %>%
    mutate(region = str_sub(variable, 14), muts = strsplit(value, ",")) %>% 
    unnest(muts) %>%
    filter(str_sub(muts, 1,1) == 'S') %>%
    mutate(pos = as.integer(str_extract(muts, '\\d+')),
           from = str_sub(str_extract(muts, '.>'),1,1), 
           to = str_sub(str_extract(muts, '>.'), 2)) %>%
    dplyr::select(-variable, -value, -muts)
}

raji.shm <- read.table('data/raji_R12.txt', header=T, sep="\t") %>%
  melt_mutations()
raji.shm$contignt <- as.character(raji.shm$contignt)
raji.shm <- mutate(raji.shm, len = nchar(contignt),
                      context = str_sub(contignt, pmax(0,pos-3), pmin(pos+3, len)),
                   cells = 'P', proj = 'raji', sample ='raji')

#YF+FLU DATASET

samples <- data.frame(type = c(rep('S', 8), rep('YF', 4)),
                      cells = c(rep('M', 6), 'P', 'P', 'M', 'M', 'P', 'P'),
                      day = c('0_1', '0_2', '14_1', '14_2', '7_1', '7_2', '14_1', '14_2', '7', '9', '7', '9'))

yff.shm <- data.frame()

for (i in 1:nrow(samples)){
  .yff.shm <- read.table(paste0('data/yf+flu/', samples$type[i], samples$cells[i], samples$day[i], '_q23.txt'), header=T, sep="\t")
  .yff.shm <- melt_mutations(.yff.shm)
  .yff.shm$proj <- ifelse(samples$type[i] == 'S', 'flu', 'yf')
  .yff.shm$cells <- samples$cells[i]
  .yff.shm$sample <- samples$day[i]
  .yff.shm$clonotype.id <- (1:nrow(.yff.shm))-1
  yff.shm <- rbind(yff.shm, .yff.shm)
}

yff.shm$contignt <- as.character(yff.shm$contignt)
yff.shm <- mutate(yff.shm, len = nchar(contignt),
                   context = str_sub(contignt, pmax(0,pos-3), pmin(pos+3, len)))

# SAVE

shm <- rbind(dplyr::select(aging.shm.1, proj, sample, cells, clonotype.id, freq, region, 
                           from=from.nt, to=to.nt, pos=pos.nt, context),
             dplyr::select(raji.shm, proj, sample, cells, clonotype.id, freq, region, 
                           from, to, pos, context),
             dplyr::select(yff.shm, proj, sample, cells, clonotype.id, freq, region, 
                           from, to, pos, context))

save(df, file = "shm.Rda", compress = T)

