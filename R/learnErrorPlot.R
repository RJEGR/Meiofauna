
# Learn error viz

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

work_path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/raw-seqs-bkp//filtN/cutadapt/Illumina/filterAndTrim/"

# setwd(path)

str(rdata <- sort(list.files(work_path, pattern="learnErrors.RData", full.names = T)))

read_rds_l <- function(rdata) {
  
  rds <- readr::read_rds(rdata)
  
  B <- sapply(strsplit(basename(rdata), "[_]"), `[`, 1)
  
  errF <- transdf_q(rds$errF)
  errR <- transdf_q(rds$errR)
  
  out_df <- transdf_q(rds$errF) %>% rbind(., transdf_q(rds$errF)) %>% mutate(Batch = B)
  
  return(out_df)
}

err_df <- lapply(rdata, read_rds_l)

do.call(bind_rows, err_df) -> err_df 

nti = c("A", "C", "G", "T")
ntj = c("A", "C", "G", "T")

err_df  %>%
  filter(from %in% nti & to %in% ntj) %>%
  ggplot(aes(x = Qual)) +
  geom_point(aes(y = Observed, color = Batch), 
    na.rm = TRUE, alpha = 1/5) + 
  geom_line(aes(y = Estimated, color = Batch), 
    linetype = "solid", size = 1, alpha = 4/5) +
  geom_line(aes(y = Nominal), 
    linetype = "dashed", size =0.5, alpha = 4/5, color = "red") +
  scale_y_log10() + 
  facet_wrap(~Transition, nrow = length(nti)) +
  labs(x = "Consensus quality score",
    y = "Error frequency (log10)",
    caption = 'Estimate (Lines) and Observed (Points) rates',
    title = 'Transitions error rate') +
  theme_bw() +
  see::scale_color_metro_d() +
  theme(legend.position = "top")

library(rstatix)

err_df %>%
  rstatix::cor_mat(vars = c("Estimated", "Input", "Nominal")) %>%
  cor_reorder() %>%
  pull_lower_triangle() %>%
  cor_plot(label = TRUE)


err_df %>% group_by(Batch) %>% cor_test(vars = "Nominal", vars2 = c("Input"))
