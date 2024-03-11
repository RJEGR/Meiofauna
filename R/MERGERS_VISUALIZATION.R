
path <- "C:/Users/Israel V/Documents/MEIOFAUNA//raw-seqs/"

setwd(path)

library(tidyverse)

mergers <- read_rds("mergers.rds")

seqtab.nochim <- read_rds("seqtab.nochim.rds")

nochim.seqs <- getSequences(seqtab.nochim)

merger_to_df <- function(x){
  
  merge_df <- dplyr::bind_rows(x, .id = "Sample")
  merge_df$seqlen <- nchar(merge_df$sequence)
  
  merge_df <- merge_df %>% as_tibble() 
  
  return(merge_df)
}


merge_df <- vector("list", length(mergers))

names(merge_df) <- names(mergers)

for (i in 1:length(mergers)) {
  merge_df[[i]] <- merger_to_df(mergers[[i]])
}

# 
# workbook <-
#   mapply(`[<-`, mergers, 'sp', value = names(mergers), SIMPLIFY = FALSE)
# 
# workbook <- do.call(rbind, workbook)

do.call(bind_rows, merge_df) %>% as_tibble() -> merge_df 

merge_df %>%
  mutate(Sample = sapply(strsplit(Sample, "-"), `[`, 1)) %>%
  count(Sample)
library("ggExtra")

col1 = "#d8e1cf" 
col2 = "#438484"

targetLength <- seq(100,150)

merge_df %>% summarise(m = min(seqlen), M = max(seqlen))

p1 <- merge_df %>%
  mutate(Sample = sapply(strsplit(Sample, "-"), `[`, 1)) %>%
  # mutate(trimmed = ifelse(seqlen %in% targetLength, 'Used', "Trimmed")) %>%
  mutate(trimmed = ifelse(sequence %in% nochim.seqs, 'nochim', "chim")) %>%
  group_by(Sample) %>% 
  mutate(abundance = abundance/sum(abundance)) %>%
  # summarise(min(abundance), max(abundance), mean(abundance))
  ggplot(aes(x = nmatch, y = seqlen)) +
  geom_point(aes(size = abundance, color = trimmed), alpha = 3/5) +
  scale_color_manual(name = NULL, values = c("#de2d26", "black")) + 
  # xlim(0,150) +
  scale_size(name = "Abundancia\nRelativa", 
             # range = c(0, 5),
             # breaks = c(0.1, 0.5, 1),
             labels = scales::percent_format(scale = 1)) +
  theme_bw(base_family='GillSans', base_size = 12) +
  theme(legend.position = "bottom") +
  xlab("Región del empalme (nt)") +
  ylab("Tamaño del amplicón (nt)")

# Marginal histogram plot

empalmePlot <- ggMarginal(p1, type = "density",
                          fill = col1, color = 'black')

empalmePlot

ggsave(empalmePlot,
       filename = 'empalme.png',
       # path = output_path,
       width = 7, height = 7,
       dpi =  300 )
