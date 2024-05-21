#this script plots posterior distributions for models

library(depmixS4)

#set the directories
input_dir <- "C:\\Users\\liber\\Desktop\\Study\\LMU\\Thesis Project - MariaCT's Lab\\Data\\HMM outputs\\Alleloscope_500kb_100nSNP\\3_state_model"
output_dir <- ".\\posterior_plots\\"

#set the working directory
setwd(input_dir)
dir.create(output_dir)

models <- readRDS("models.RDS")

#iterate through each model and cell
pdf(paste(output_dir, "posterior_plots.pdf"))
par(mfrow = c(3, 1))
for (model in names(models)){
  for (cell in names(models[[model]][["fitted_models"]])){
    fm <- models[[model]][["fitted_models"]][[cell]]
    
    response_mtx <- summary(fm, which = "response")
    
    cnv_state_vector <- c()
    cnv_state_vector["loss"] <- which(response_mtx[,"Re1.(Intercept)"] == min(response_mtx[,"Re1.(Intercept)"]))
    cnv_state_vector["gain"] <- which(response_mtx[,"Re1.(Intercept)"] == max(response_mtx[,"Re1.(Intercept)"]))
    cnv_state_vector["base"] <- which(!(c(1,2,3) %in% cnv_state_vector))
    
    posterior_df <- models[[model]][["fitted_models"]][[cell]]@posterior
    
    for (cnv_state in c("loss", "base", "gain")){
      hist(posterior_df[,cnv_state_vector[cnv_state]+1],
                main = paste(cnv_state, " state posterior probabilities distribution in ",
                        model, " model (cell: ", cell, ")", sep = ""),
                xlab = "Posterior Probability",
                ylab = "Number of Regions")
      # plot_list[[paste(model, cell, cnv_state, sep = "_")]] <- p
    }
  }
}
dev.off()
