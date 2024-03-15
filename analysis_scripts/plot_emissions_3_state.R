#this script plots emission distributions for models

library(depmixS4)

#set the directories
input_dir <- "C:\\Users\\liber\\Desktop\\Study\\LMU\\Thesis Project - MariaCT's Lab\\Data\\HMM outputs\\Alleloscope_500kb_100nSNP\\3_state_model"
output_dir <- ".\\emission_plots\\"

#set the working directory
setwd(input_dir)
dir.create(output_dir)

models <- readRDS("models.RDS")

color_set <- c("#1362ce", "#13ce42", "#ce1313")
names(color_set) <- c("loss", "base", "gain")

#iterate through each model and cell
for (model in names(models)){
  pdf(paste(output_dir, paste(model, "_emission_plots.pdf")))
  if (model != "counts"){
    par(mfrow = c(2, 1))
  }
  for (cell in names(models[[model]][["fitted_models"]])){
    fm <- models[[model]][["fitted_models"]][[cell]]
    
    response_mtx <- summary(fm, which = "response")
    
    #determine which class is which CNV
    cnv_state_vector <- c()
    cnv_state_vector["loss"] <- which(response_mtx[,"Re1.(Intercept)"] == min(response_mtx[,"Re1.(Intercept)"]))
    cnv_state_vector["gain"] <- which(response_mtx[,"Re1.(Intercept)"] == max(response_mtx[,"Re1.(Intercept)"]))
    cnv_state_vector["base"] <- which(!(c(1,2,3) %in% cnv_state_vector))
    
    #extract response parameters
    theta_data <- c()
    counts_data <- models[[model]][["data"]]$counts
    if (model != "counts"){
      theta_data <- models[[model]][["data"]]$theta_hat
      if (model == "counts_and_theta"){
        theta_data <- theta_data[!is.na(theta_data)]
      }
    }
    
    #plot the distributions for counts
    density_counts <- density(counts_data, )
    y_max <- max(density_counts$y)
    plot(density_counts, lwd = 2, col = "white", ylim = c(0, y_max*2.5),
         main = paste("Counts emission distribution in ",
                      model, " model\n(cell: ", cell, ")", sep = ""),
         xlab = "Count",
         ylab = "Portion of Regions")
    
    for (cnv_state in c("loss", "base", "gain")){
      curve(dnorm(x, mean = response_mtx[cnv_state_vector[cnv_state],1], 
                  sd = response_mtx[cnv_state_vector[cnv_state],2]), 
            add = TRUE, col = color_set[cnv_state], lwd = 2,
            ylab = "Portion of Regions", xlab = "Count",
            main = paste("Counts emission distributions for ", model, " model\n(cell: ",
                         cell, ")", sep = ""))
    }
    
    #plot the distributions for theta if applicable
    if (model != "counts"){
      density_theta <- density(theta_data)
      y_max <- max(density_theta$y)
      plot(density_theta, lwd = 2, col = "white", ylim = c(0, y_max*2.5),
           main = paste("AF emission distribution in ",
                        model, " model\n(cell: ", cell, ")", sep = ""),
           xlab = "AF",
           ylab = "Portion of Regions")
      
      for (cnv_state in c("loss", "base", "gain")){
        curve(dnorm(x, mean = response_mtx[cnv_state_vector[cnv_state],3], 
                    sd = response_mtx[cnv_state_vector[cnv_state],4]), 
              add = TRUE, col = color_set[cnv_state], lwd = 2,
              ylab = "Portion of Regions", xlab = "AF",
              main = paste("AF emission distributions for ", 
                           model, " model\n(cell: ", cell, ")", sep = ""))
      }
    }
  }
  dev.off()
}