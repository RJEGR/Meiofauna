
# PERFORM ANCOM TEST

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "~/MEIOFAUNA/INPUTS/"

library(phyloseq)
library(microbiome)
library(tidyverse)

out_df <- function(out_ancombc) {
  
  res = out_ancombc$res
  
  df_fig0 <- res %>% select_at(vars(starts_with(c("taxon", "diff")))) %>%
    pivot_longer(-taxon, names_to = "group", values_to = "diff") %>%
    mutate(group = gsub("diff_", "", group), diff = as.integer(diff))
  
  df_fig1 <- res %>% 
    select_at(vars(starts_with(c("taxon","lfc_")))) %>% 
    pivot_longer(-taxon, names_to = "group", values_to = "logFC") %>%
    mutate(group = gsub("lfc_", "", group))
  
  df_fig2 <- res %>% 
    select_at(vars(starts_with(c("taxon","se_")))) %>% 
    pivot_longer(-taxon, names_to = "group", values_to = "SE") %>%
    mutate(group = gsub("se_", "", group))
  
  # res %>% 
  #   select_at(vars(starts_with(c("taxon","W_")))) %>% 
  #   pivot_longer(-taxon, names_to = "group", values_to = "W") %>%
  #   mutate(group = gsub("W_", "", group))
  
  df_fig3 <- res %>% 
    select_at(vars(starts_with(c("taxon","q_")))) %>% 
    pivot_longer(-taxon, names_to = "group", values_to = "q_val") %>%
    mutate(group = gsub("q_", "", group)) %>%
    mutate(star = ifelse(q_val <.001, "***", 
                         ifelse(q_val <.01, "**",
                                ifelse(q_val <.05, "*", ""))))
  
  df_fig0 %>% left_join(df_fig1) %>% left_join(df_fig2) %>% left_join(df_fig3)
}

phyloseq <- read_rds(paste0(wd, '/phyloseq.rds'))

# set the contrast (look at the contrasts)
sample_data(phyloseq) %>% with(., table(Region))

ancombc_data1 <-  phyloseq %>%
  # subset_samples(Tissue=="Foregut") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .) # %>%
  # aggregate_taxa(., "p")

formula <- "Region"

# .out1 <- ANCOMBC::ancombc(data = ancombc_data1, formula = formula,
#                          p_adj_method = "holm",
#                          prv_cut = 0.2, # zero_cut
#                          lib_cut = 100,
#                          group =  formula, struc_zero = TRUE, neg_lb = TRUE,
#                          tol = 1e-5, max_iter = 1000, conserve = TRUE,
#                          alpha = 0.05, global = FALSE)
# 
# out_df(.out1)


out1 <- ANCOMBC::ancombc2(data = ancombc_data1, fix_formula = formula,
                         group =  formula, tax_level = "p", pairwise = T)

recode_to <- c("(Intercept)" = "Yucatan")

bar_df <- out_df(out1) %>% 
  filter(diff == 1) %>%
  mutate(group = gsub("Region", "", group)) %>%
  # mutate(wrap = paste0("Yucatan-", group)) %>%
  mutate(group = dplyr::recode_factor(group, !!!recode_to)) %>%
  mutate(
    y_star = logFC + (0.2+SE)*sign(logFC),
    ymin = (abs(logFC) - SE) * sign(logFC),
    ymax = (abs(logFC) + SE) * sign(logFC)) 

bar_df %>%
  ggplot(data = ., 
    aes(x = reorder(taxon, desc(taxon)), y = logFC, fill = group)) + 
  geom_bar(stat = "identity", width = 0.5, 
    position = position_dodge(width = 0.4)) +
  ggh4x::facet_nested(  ~ Region, scales = "free", space = "free", switch = "y") +
  coord_flip() +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2,
    position = position_dodge(0.05), color = "black") + 
  geom_text(aes(y = y_star, label=star), 
    vjust=.7, color="black", position=position_dodge(width = .5)) +
  ggsci::scale_fill_aaas() 


ancombc_data1 %>% aggregate_taxa(., "p") %>%  
  transform_sample_counts(function(x) sqrt(x / sum(x))) %>% 
  phyloseq::plot_heatmap() +
  facet_grid(~ Region, scales = "free_x")

   