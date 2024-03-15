library(dplyr)
library(tidyr)

group_1 <- c("GTGATCAAGTAACATG.1", "TGACAACGTTTGCCCT.1", "ACTATTCTCTATACCT.1")
group_2 <- c("TGCATGATCGAGAACG.1", "TTTGAGGAGAGGAATG.1", "ACAGACTTCCACGCTT.1")
group_3 <- c("TTAGCGAAGTCTAGAA.1", "TGGTCCTCACATAAAG.1", "GGAGGATGTTCCGCGA.1")

AF_metrics <- models[["counts_and_theta"]][["metrics"]]
AF_imp_metrics <- models[["counts_and_theta_imp"]][["metrics"]]

AF_g1_metrics <- list()
AF_g2_metrics <- list()
AF_g3_metrics <- list()

AF_imp_g1_metrics <- list()
AF_imp_g2_metrics <- list()
AF_imp_g3_metrics <- list()

for (i in names(AF_metrics)){
  AF_g1_metrics[[i]] <- AF_metrics[[i]][rownames(AF_metrics[[i]]) %in% group_1,]
  AF_g2_metrics[[i]] <- AF_metrics[[i]][rownames(AF_metrics[[i]]) %in% group_2,]
  AF_g3_metrics[[i]] <- AF_metrics[[i]][rownames(AF_metrics[[i]]) %in% group_3,]
  
  AF_imp_g1_metrics[[i]] <- AF_imp_metrics[[i]][rownames(AF_imp_metrics[[i]]) %in% group_1,]
  AF_imp_g2_metrics[[i]] <- AF_imp_metrics[[i]][rownames(AF_imp_metrics[[i]]) %in% group_2,]
  AF_imp_g3_metrics[[i]] <- AF_imp_metrics[[i]][rownames(AF_imp_metrics[[i]]) %in% group_3,]
}
