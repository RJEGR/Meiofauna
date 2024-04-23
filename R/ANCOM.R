
# PERFORM ANCOM TEST

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/"

library(phyloseq)
library(microbiome)
library(tidyverse)
library(ANCOMBC)

out_df <- function(out_ancombc) {
  
  res_pair = out1$res_pair
  
  
  res_pair %>% select_at(vars(starts_with(c("diff"))))
  res_pair %>% filter_at(vars(starts_with("diff")), any_vars(. == 1))
  
  df_fig0 <- res_pair %>% select_at(vars(starts_with(c("taxon", "diff")))) %>%
    pivot_longer(-taxon, names_to = "group", values_to = "diff") %>%
    mutate(group = gsub("diff_", "", group), diff = as.integer(diff))
  
  df_fig1 <- res_pair %>% 
    select_at(vars(starts_with(c("taxon","lfc_")))) %>% 
    pivot_longer(-taxon, names_to = "group", values_to = "logFC") %>%
    mutate(group = gsub("lfc_", "", group))
  
  df_fig2 <- res_pair %>% 
    select_at(vars(starts_with(c("taxon","se_")))) %>% 
    pivot_longer(-taxon, names_to = "group", values_to = "SE") %>%
    mutate(group = gsub("se_", "", group))
  
  df_fig3 <- res_pair %>% 
    select_at(vars(starts_with(c("taxon","q_")))) %>% 
    pivot_longer(-taxon, names_to = "group", values_to = "q_val") %>%
    mutate(group = gsub("q_", "", group)) %>%
    mutate(star = ifelse(q_val <.001, "***", 
      ifelse(q_val <.01, "**",
        ifelse(q_val <.05, "*", ""))))
  
  df <- df_fig0 %>% left_join(df_fig1) %>% left_join(df_fig2) %>% left_join(df_fig3)
  
  df <- df %>% dplyr::mutate(label = ifelse(diff == 1, round(logFC, 2), NA)) %>% dplyr::arrange(taxon)
  
  return(df)
}

phyloseq <- read_rds(paste0(wd, '/phyloseq.rds'))

# Generate the true abundances
# abn_data = sim_plnm(abn_table = as(otu_table(phyloseq), "matrix"), taxa_are_rows = FALSE, prv_cut = 0.05, 
#   n = 160, lib_mean = 1e8, disp = 0.5)

# set the contrast (look at the contrasts)
sample_data(phyloseq) %>% with(., table(Region))

reg_levels <- c("Deep-sea", "NW Slope", "NW Shelf", "Yucatan")


reframe <- as(sample_data(phyloseq), "data.frame") %>% mutate(Region = factor(Region, levels = reg_levels))

sample_data(phyloseq) <- reframe

ancombc_data1 <-  phyloseq %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)  %>%
  microbiome::aggregate_taxa(., "p") # as to wgcna, ancom perform good resulsts to p

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
                         p_adj_method = "holm",
                         global = TRUE, dunnet = TRUE, 
                         trend = FALSE, trend_control = NULL,
                         group =  formula, tax_level = NULL, pairwise = T)



recode_to <- c("Yucatan_NW Shelf" = "Yucatan-NW Shelf ", 
               "Yucatan_NW Slope"  = "Yucatan-NW Slope",
               "NW Shelf_NW Slope" = "NW Shelf-NW Slope")

out_df(out1) %>% count(group)

bar_df <- out_df(out1) %>% 
  # filter(diff == 1) %>%
  mutate(group = gsub("Region", "", group)) %>%
  mutate(wrap = paste0("Deep-sea:", group)) %>%
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
  ggh4x::facet_nested(  ~ group, scales = "free", space = "free", switch = "y") +
  coord_flip() +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2,
    position = position_dodge(0.05), color = "black") + 
  geom_text(aes(y = y_star, label=star), 
    vjust=.7, color="black", position=position_dodge(width = .5)) +
  ggsci::scale_fill_aaas() +
  theme_classic(base_size = 14, base_family = "GillSans")



lo = floor(min(df$logFC))
up = ceiling(max(df$logFC))
mid = (lo + up)/2

bar_df %>%
  ggplot(aes(x = group, y = taxon, fill = logFC)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up),
    name = NULL) +
  geom_text(aes(group, taxon, label = label), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to multiple c. subjects") +
  theme_minimal(base_size = 14, base_family = "GillSans") +
  theme(plot.title = element_text(hjust = 0.5))



# ancombc_data1 %>% aggregate_taxa(., "c") %>%  
#   transform_sample_counts(function(x) sqrt(x / sum(x))) %>% 
#   phyloseq::plot_heatmap() +
#   facet_grid(~ Region, scales = "free_x")

   