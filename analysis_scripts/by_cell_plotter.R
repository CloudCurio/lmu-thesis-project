per_cell_summary <- function(models, output_dir){
  ##############################################################################
  # loading libraries
  ##############################################################################
  
  library(depmixS4)
  
  library(dplyr)
  library(tidyr)
  
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(gridGraphics)
  
  ##############################################################################
  # Initiating important variables
  ##############################################################################
  
  #extact unique cell names into a vector
  cells_vector <- unique(models[["counts"]][["data"]]$cell)
  
  #create an object to save intermediary results into
  output_list <- list()
  
  ##############################################################################
  # create a heatmap of CNV predictions
  ##############################################################################
  
  #iterate through cells to create individual heatmaps
  for (cell in cells_vector){
    #create a cell entry in the output_list
    output_list[[cell]] <- list()
    
    #Create a dataframe with input data
    HMM_preds <- NULL
    
    #iterate through models, adding prediction data to the HMM_preds
    for (model in names(models)){
      #subset data by model and cell
      subset_df <- models[[model]][['predictions']]
      subset_df <- subset_df[which(subset_df$cell == cell),]
      subset_df$model <- model
      
      #split the region column into *chr* label and region starting position
      subset_df$chr <- gsub("chr([0-9]+):.*", "chr\\1", subset_df$region)
      subset_df$fragment <- as.numeric(gsub("chr[0-9]+:([0-9]+)", "\\1", subset_df$region))
      
      #merge data for same cell but different models
      if (is.null(HMM_preds)){
        HMM_preds <- subset_df
      } else {
        HMM_preds <- rbind(HMM_preds, subset_df)
      }
    }
    #add WGS data as a subset to the HMM_preds
    WGS_info <- subset_df
    WGS_info$HMM_pred <- WGS_info$WGS_score
    WGS_info$model <- "WGS"
    HMM_preds <- rbind(HMM_preds, WGS_info)
    
    #order the dataframe properly
    #make chr column a factor
    HMM_preds$chr <- factor(HMM_preds$chr, levels = paste0("chr", 1:22))
    
    #reorder the dataframe
    HMM_preds <- HMM_preds[order(HMM_preds$chr, HMM_preds$fragment), ]
    
    #add X and Y positions and record positions for labels
    HMM_preds <- HMM_preds %>%
      mutate(X_pos = as.numeric(factor(region, levels = unique(region))))
    
    chromosome_label_positions <- HMM_preds %>%
      filter(!duplicated(chr)) %>%
      pull(X_pos)
    
    #create a heatmap of predictions and WGS info
    p <- ggplot(HMM_preds, aes(x = X_pos, y = model, fill = factor(HMM_pred))) +
      geom_tile() +
      geom_vline(xintercept = chromosome_label_positions - 0.5, color = "black", linetype = "dotted", size = 0.5) +
      geom_hline(yintercept = 3.5, color = "black", linetype = "solid", size = 1) +
      labs(title = paste("HMM CNV predictions (cell: ", cell, ")", sep = ""),
           x = "Region",
           y = "Model",
           fill = "CNV class") +
      scale_x_continuous(breaks = chromosome_label_positions, labels = unique(HMM_preds$chr),
                         expand = c(0.05, 0)) +
      scale_fill_manual(values = c("#2171B5", "#1fb424", "#b4421f"),
                        labels = c("loss", "base", "gain")) +
      theme_minimal() +
      theme(legend.position = "top", 
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) 
    
    #save the heatmap in the corresponding cell's entry
    output_list[[cell]][["heatmap"]] <- p
  }
  
  ##############################################################################
  # create a table with performance metrics
  ##############################################################################
  
  for (cell in cells_vector){
    #create an empty dataframe for metrics storage
    metrics_df <- as.data.frame(matrix(ncol = 3, nrow = 4))
    colnames(metrics_df) <- names(models)
    rownames(metrics_df) <- c("Overall_Acc", "Loss_Acc", "Base_Acc", "Gain_Acc")
    
    #iterate through models and extract the metrics for a particular cell
    for (model in names(models)){
      Overall_Acc <- round(models[[model]][['metrics']][['overall']][cell, "overall_Acc"], 3)
      Loss_Acc <- round(models[[model]][['metrics']][['loss']][cell, "Balanced Accuracy"], 3)
      Base_Acc <- round(models[[model]][['metrics']][['base']][cell, "Balanced Accuracy"], 3)
      Gain_Acc <- round(models[[model]][['metrics']][['gain']][cell, "Balanced Accuracy"], 3)
      
      metrics_vector <- c(Overall_Acc, Loss_Acc, Base_Acc, Gain_Acc)
      metrics_df[,model] = metrics_vector
    }
    
    #save the metrics dataframe in the corresponding cell's entry
    output_list[[cell]][["metrics"]] <- metrics_df
  }
  
  ##############################################################################
  # create emission distribution plots
  ##############################################################################
  
  #designate colors to use
  color_set <- c("#1362ce", "#13ce42", "#ce1313")
  names(color_set) <- c("loss", "base", "gain")
  
  for (cell in cells_vector){
    #create an entry for the distribution plots
    output_list[[cell]][["emission_plots"]] <- list()
    
    for (model in names(models)){
      if (model == "counts"){
        next
      }
      
      #extract the fitted HMM model
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
      theta_data <- models[[model]][["data"]]$theta_hat
      if (model == "counts_and_theta"){
        theta_data <- theta_data[!is.na(theta_data)]
      }
      
      #plot counts distributions
      counts_plot <- ggplot(data = data.frame(x = c(0:max(counts_data))), aes(x = x, colour = cnv_tag)) +
        stat_function(fun = dnorm, n = 10001, 
                      args = list(mean = response_mtx[cnv_state_vector["loss"],1], 
                                  sd = response_mtx[cnv_state_vector["loss"],2]),
                      colour = color_set["loss"]) + ylab("") +
        stat_function(fun = dnorm, n = 10001, 
                      args = list(mean = response_mtx[cnv_state_vector["base"],1], 
                                  sd = response_mtx[cnv_state_vector["base"],2]),
                      colour = color_set["base"]) + ylab("") +
        stat_function(fun = dnorm, n = 10001, 
                      args = list(mean = response_mtx[cnv_state_vector["gain"],1], 
                                  sd = response_mtx[cnv_state_vector["gain"],2]),
                      colour = color_set["gain"]) + ylab("") +
        theme(axis.line = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              axis.title = element_text(size = 12)) +
        labs(x = "Count", y = "Portion of Regions", 
             title = paste("Counts emission distribution\n(model: ", model, ")", sep = ""))
      
      #plot AF distributions
      AF_plot <- ggplot(data = data.frame(x = c(0:1)), aes(x)) +
        stat_function(fun = dnorm, n = 10001, 
                      args = list(mean = response_mtx[cnv_state_vector["loss"],3], 
                                  sd = response_mtx[cnv_state_vector["loss"],4]),
                      colour = color_set["loss"]) + ylab("") +
        stat_function(fun = dnorm, n = 10001, 
                      args = list(mean = response_mtx[cnv_state_vector["base"],3], 
                                  sd = response_mtx[cnv_state_vector["base"],4]),
                      colour = color_set["base"]) + ylab("") +
        stat_function(fun = dnorm, n = 10001, 
                      args = list(mean = response_mtx[cnv_state_vector["gain"],3], 
                                  sd = response_mtx[cnv_state_vector["gain"],4]),
                      colour = color_set["gain"]) + ylab("") +
        theme(axis.line = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              axis.title = element_text(size = 12)) +
        labs(x = "Allele Frequency", y = "Portion of Regions", 
             title = paste("AF emission distribution\n(model: ", model, ")", sep = ""))
      
      #create model storage lists and save the plots into them
      output_list[[cell]][["emission_plots"]][[model]][["counts"]] <- counts_plot
      output_list[[cell]][["emission_plots"]][[model]][["AF"]] <- AF_plot
    }
  }
  ##############################################################################
  # combine data in an output file
  ##############################################################################
  
  #open a pdf device
  pdf(paste(output_dir, "per_cell_summary.pdf"), width = 8.3, height = 11.7)
  
  #make a separate page per cell
  for (cell in cells_vector){
    
    #merge emission distribution plots into a patch
    emission_patch <- output_list[[cell]][["emission_plots"]][[1]][["counts"]] + 
      output_list[[cell]][["emission_plots"]][[2]][["counts"]] + 
      output_list[[cell]][["emission_plots"]][[1]][["AF"]] +
      output_list[[cell]][["emission_plots"]][[2]][["AF"]] + 
      plot_layout(nrow = 2, byrow = TRUE)
    
    #join everything together
    final_patch <- output_list[[cell]][["heatmap"]] / 
      gridExtra::tableGrob(output_list[[cell]][["metrics"]]) / 
      emission_patch +
      plot_layout(heights = c(0.75,1,1.8))
    
    print(final_patch)
  }
  #close the device
  dev.off()
  
}

