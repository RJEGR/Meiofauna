
mergers <- read_rds("mergers.rds")

merger_to_df <- function(merger){
  merge_df <- bind_rows(merger, .id = "Sample")
  merge_df$seqlen <- nchar(merge_df$sequence)
  
  return(merge_df)
}

dim(merge_df <- merger_to_df(mergers))

library("ggExtra")

col1 = "#d8e1cf" 
col2 = "#438484"

targetLength <- seq(300,320)

nochim.seqs <- getSequences(seqtab.nochim)

p1 <- merge_df %>%
  mutate(Sample = sapply(strsplit(Sample, "[_]"), `[`, 1)) %>%
  #mutate(trimmed = ifelse(seqlen %in% targetLength, 'Used', "Trimmed")) %>%
  mutate(trimmed = ifelse(sequence %in% nochim.seqs, 'nochim', "chim")) %>%
  group_by(Sample) %>% 
  # mutate(abundance = raTrans(abundance)) %>% 
  # summarise(min(abundance), max(abundance), mean(abundance))
  ggplot(aes(x = nmatch, y = seqlen)) +
  geom_point(aes(size = abundance, color = trimmed), alpha = 3/5) +
  scale_color_manual(name = NULL, values = c("#de2d26", "black")) + 
  xlim(0,150) +
  scale_size(name = "Abundancia\nRelativa", 
             range = c(0, 5),
             breaks = c(1, 5, 10, 20, 30),
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
       path = output_path,
       width = 7, height = 7,
       dpi =  300 )
