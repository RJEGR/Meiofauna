
# path <- "C:/Users/Israel V/Documents/MEIOFAUNA//raw-seqs/"
path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/"

setwd(path)

library(tidyverse)
library(dada2)

mergers <- read_rds("mergers.rds")

seqtab.nochim <- read_rds("seqtab.nochim.rds")

nochim.seqs <- getSequences(seqtab.nochim)

merger_to_df <- function(x){
  
  test_if_emtpy_df <-  function(x){
    
    if(nrow(x)>0){
      test <- TRUE
    }
    else{
      test<- FALSE
    }
    
    return(test)
  }
  
  x <- x[unlist(lapply(x, test_if_emtpy_df))]
  
  merge_df <- dplyr::bind_rows(x, .id = "Sample")
  
  merge_df$seqlen <- nchar(merge_df$sequence)
  
  merge_df <- merge_df %>% as_tibble() 
  
  return(merge_df)
}

# merger_to_df(mergers$PE1)
  
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
  mutate(Cruise = sapply(strsplit(Sample, "-"), `[`, 1)) %>%
  group_by(Cruise) %>%
  summarise(n_asvs = n(), Tabundance = sum(abundance), n_sam = length(unique(Sample))) 

library(ggdensity)

merge_df %>%
  mutate(trimmed = ifelse(sequence %in% nochim.seqs, 'nochim', "chim")) %>%
  ggplot(aes(y = nmatch, x = seqlen)) +
  geom_hdr_lines() +
  theme_bw(base_family='GillSans', base_size = 16) +
  theme(legend.position = "bottom") +
  ylab("Merged region (nt)") +
  xlab("Amplicon Size (nt)") +
  ylim(50, 160)

library("ggExtra")

col1 = "#d8e1cf" 
col2 = "#438484"

targetLength <- seq(0,0)

merge_df %>% summarise(m = min(seqlen), M = max(seqlen))

library(ggdensity)

p1 <- merge_df %>%
  mutate(Sample = sapply(strsplit(Sample, "_"), `[`, 1)) %>%
  # mutate(trimmed = ifelse(seqlen %in% targetLength, 'Used', "Trimmed")) %>%
  mutate(trimmed = ifelse(sequence %in% nochim.seqs, 'nochim', "chim")) %>%
  group_by(Sample) %>% 
  mutate(abundance = abundance/sum(abundance)) %>%
  # summarise(min(abundance), max(abundance), mean(abundance))
  ggplot(aes(x = nmatch, y = seqlen)) +
  # geom_hdr() +
  # geom_point(aes(size = abundance, color = trimmed), alpha = 3/5) +
  geom_point(aes( color = trimmed), 
    shape = 21, size = 1.5, stroke = 1.2, fill = "white") +
  scale_color_manual(name = NULL, values = c("#de2d26", "black"))  +
  # xlim(0,150) +
  scale_size(name = "Relative\nAb.",
             # range = c(0, 5),
             # breaks = c(0.1, 0.5, 1),
             labels = scales::label_percent()
    ) +
  theme_bw(base_family='GillSans', base_size = 12) +
  theme(legend.position = "bottom") +
  xlab("Merged region (nt)") +
  ylab("Amplicon Size (nt)")

# Marginal histogram plot

empalmePlot <- ggMarginal(p1, type = "density",
                          fill = col1, color = 'black')

empalmePlot

ggsave(empalmePlot,
       filename = 'multi_run_empalme_apr9.png',
       # path = output_path,
  device = png,
       width = 5, height = 5,
       dpi =  300 )

# TRACK PLOT -----

out <- read_tsv(list.files(pattern = "_filterAndTrim.tsv"))

out %>% mutate(Sample = sapply(strsplit(Sample, "-"), `[`, 1)) %>% group_by(Sample) %>% tally(reads.out) %>% view()

getN <- function(x) sum(getUniques(x))

# sapply(mergers[[1]], getN) %>% dplyr::bind_rows() %>% pivot_longer(names(.))
       
merge_out <- vector("list", length(mergers))

names(merge_out) <- names(mergers)

for (i in 1:length(mergers)) {
  merge_out[[i]] <- sapply(mergers[[i]], getN) %>% dplyr::bind_rows() %>% 
    pivot_longer(names(.), names_to = "Sample", values_to = "merged") %>%
    mutate(Sample = sapply(strsplit(Sample, "[_]"), `[`, 1))
  
}

do.call(bind_rows, merge_out) -> merge_out 


# track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))

seqtab <- read_rds("seqtab.collapsed.rds")

track <- data.frame(tabled = rowSums(seqtab), nonchim = rowSums(seqtab.nochim)) %>% 
  dplyr::as_tibble(rownames = "Sample") %>%
  mutate(Sample = sapply(strsplit(Sample, "[_]"), `[`, 1)) 

track %>% mutate(Sample = sapply(strsplit(Sample, "-"), `[`, 1)) %>% group_by(Sample) %>% tally(nonchim) %>% view()


track <- out %>% left_join(merge_out) %>% left_join(track)

colnames(track) <- c("Sample", "input", "filtered", "merged", "tabled", "nonchim")

track <- track %>% mutate_if(is.numeric, function(x) {x / track$input})

track <- track %>% pivot_longer(-Sample, values_to = "Reads", names_to = "Step")

track$Step <- factor(track$Step, 
                          levels = c("input", "filtered", "merged", "tabled", "nonchim"),
                          ordered=T)


track <- track %>% mutate(batch = sapply(strsplit(Sample, "-"), `[`, 1))

# order samples by size
track$Sample <- factor(track$Sample, levels = out[order(-out$reads.in),]$Sample,
                         ordered = T) 

track_plot <- ggplot(track, aes(x=Sample, y=Reads, fill=Step)) + 
  geom_col(position = position_identity(), width = 0.8, aes(fill = Step), alpha = 0.7) +
  scale_fill_brewer(palette = "BrBG", direction = 1) +
  labs(title=paste0(run, ". % Reads count through processing"), x="") +
  theme_bw(base_family='GillSans', base_size = 12) +
  theme(legend.position = "top", 
        # axis.text.x = element_text(angle = 45, hjust=1, size = 5)
        axis.text.x = element_blank())
  # coord_flip() 

track_plot

track_plot <- track_plot + facet_grid(~ batch , scales = "free", space = "free")

ggsave(track_plot,
       device = png,
       filename = 'multi_run_track-bar.png',
       # path = output_path,
       width = 20, height = 4,
       dpi =  300 )

p <- track %>%
  ggplot(aes(y = Step, x =  Sample, fill = Reads)) +
  geom_tile(color = 'white', linewidth = 0) +
  scale_fill_viridis_c(option = "B", name = "% RA", direction = -1, na.value = "white") +
  # ggsci::scale_fill_material("deep-orange")
  facet_grid(~ batch , scales = "free", space = "free") +
  theme_bw(base_family='GillSans', base_size = 12) +
  theme(legend.position = "top", 
        # axis.text.x = element_text(angle = 45, hjust=1, size = 5)
        axis.text.x = element_blank())



ggsave(p,
       device = png,
       filename = 'multi_run_track.png',
       # path = output_path,
       width = 20, height = 4,
       dpi =  300 )

