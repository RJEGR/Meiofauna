
# PERFORM ANCOM TEST

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "~/MEIOFAUNA/INPUTS/"

library(phyloseq)
library(microbiome)


out_df <- function(out_ancombc) {
  
  res = out_ancombc$res
  
  df_fig1 = data.frame(res$W * res$diff_abn, check.names = FALSE) %>% 
    rownames_to_column("taxon_id") %>%
    pivot_longer(-taxon_id, names_to = "group", values_to = "logFC")
  
  df_fig2 = data.frame(res$se * res$diff_abn, check.names = FALSE) %>% 
    rownames_to_column("taxon_id") %>% 
    pivot_longer(-taxon_id, names_to = "group", values_to = "SE")
  
  df_fig3 = data.frame(res$q_val * res$diff_abn, check.names = FALSE) %>% 
    rownames_to_column("taxon_id") %>%
    pivot_longer(-taxon_id, names_to = "group", values_to = "q_val") %>%
    mutate(star = ifelse(q_val <.001, "***", 
                         ifelse(q_val <.01, "**",
                                ifelse(q_val <.05, "*", ""))))
  
  # table(df_fig3$star)
  # colnames(df_fig2)[-1] = paste0(colnames(df_fig2)[-1], "SD")
  
  df_fig1 %>% left_join(df_fig2) %>% left_join(df_fig3)
}


phyloseq <- read_rds(paste0(wd, '/phyloseq.rds'))

ancombc_data1 <-  phyloseq %>%
  # subset_samples(Tissue=="Foregut") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .) # %>%
  # aggregate_taxa(., "p")

formula <- "Region" 

.out1 <- ANCOMBC::ancombc(data = ancombc_data1, formula = formula,
                         p_adj_method = "holm",
                         prv_cut = 0.2, # zero_cut
                         lib_cut = 100,
                         group =  formula, struc_zero = TRUE, neg_lb = TRUE,
                         tol = 1e-5, max_iter = 1000, conserve = TRUE,
                         alpha = 0.05, global = FALSE)

out_df(.out1)


out1 <- ANCOMBC::ancombc2(data = ancombc_data1, fix_formula = formula,
                         group =  formula, tax_level = "p")

out_ancombc <- out1

