# Plot time series per chain
figure1B <- function(data){
  # Extract the epi and sequence data
  dt_cases <- data[[1]]
  df_seq <- data[[2]]
  nchain <- length(unique(dt_cases$chain))
  # Compute the weekly number of cases
  week_dates <- lubridate::floor_date(dt_cases$date_ons, "week")
  # Convert week_dates to a factor 
  all_dates <- factor(week_dates, levels = as.character(
    seq(min(week_dates), as.Date("2022-02-07"), "week")))
  month_dates <- c(lubridate::floor_date(c(as.Date(as.character(levels(all_dates)))), "month"))
  ymax <- max(table(all_dates))
  # For each chain, plot the weekly number of cases
  par(mfrow = c(nchain, 1), mar = c(2, 4, 0, 2), las = 1, bty = "l", 
      cex.axis = 1.1, cex.main = 1.3, cex.lab = 1.3, oma = c(0,0,0,0))
  cols <- c("lightblue", "lightgreen", "darkred", "gold", "green", "pink", "darkblue", "purple")
  for(i in seq_along(unique(dt_cases[order(date_ons), chain]))){
    # Extract the onset dates of the cases in chain_i
    chain_i <- unique(dt_cases[order(date_ons, decreasing = T),chain])[i]
    dates_i <- all_dates[dt_cases$chain == chain_i]
    dates_i_seq <- all_dates[dt_cases$chain == chain_i & df_seq$seq != ""]
    b <- barplot(table(dates_i) ~ as.Date(names(table(dates_i))), ylim = c(0, ymax), 
                 col = NA, xlab = "", ylab = "", xaxt = "n")
    barplot(table(dates_i_seq) ~ as.Date(names(table(dates_i_seq))), 
            ylim = c(0, ymax), col = cols[i], add = T, xaxt = "n")
    # Add title and axis
    title(adj = 1, paste0("Chain ", nchain + 1 - i), line = -2)
    num_dates <- b[,1]
    index <- c(seq(1, length(unique(month_dates)), 3))
    axis(1, num_dates[!duplicated(month_dates)][index], rep(" ", length(index)))
  }
  
  axis(1, num_dates[!duplicated(month_dates)][seq(1, length(num_dates), 3)],
       unique(format(unique(month_dates), "%b %Y"))[seq(1, length(num_dates), 3)])
  title(xlab = "Date of Symptom Onset")
  title(ylab = "Weekly Case Count", line = - 1.5, outer = T)
  
}

# Generate alpha plots
figure3 <- function(out_wt, out_comb, data, burnin){
  # Define colour scheme
  cols <- c("#88CCEE","#44AA99","#AA4499","#999933","#117733",
            "#882255","#6699CC")
  # Select the cases from the epi dataset that were sequenced 
  dt_cases <- data[[1]][data[[2]]$seq != ""][order(chain, date_ons)]
  # Re-order chain from 1 to length(unique(dt_cases$chain))
  dt_cases$chain <- as.numeric(factor(dt_cases$chain))
  chains <- unique(dt_cases$chain)
  # Generate epi and combined ggplot (without plotting)
  plot_WT <- plot(out_wt, type = "alpha", burnin = burnin)
  plot_O2 <- plot(out_comb, type = "alpha", burnin = burnin)
  
  # Horizontal and vertical lines corresponding to each chan
  v_lines <- c(0,cumsum(table(dt_cases$chain)) + .5)
  
  # Convert columns of plot_WT and plot_O2
  data_plot_WT <- convert_data_plot(plot_WT, dt_cases, chains)
  data_plot_O2 <- convert_data_plot(plot_O2, dt_cases, chains)

  data_plot_WT$Type <- "A"
  data_plot_O2$Type <- "B"
  
  # Merge data_plot_WT and data_plot_O2
  data_plot <- rbind(data_plot_O2, data_plot_WT)
  data_plot$Type <- factor(data_plot$Type, levels = c("A", "B"))
  
  # Alpha plots
  plot_tot <- ggplot() + geom_point(data = data_plot, 
                                    aes(x = Infectee, y = Infector, color = Chain,
                                        size = Frequency)) + 
    facet_grid(~Type) + 
    scale_color_manual(values = c(cols, "grey")) + theme_classic() + 
    theme(text = element_text(size = 15),
          panel.spacing = unit(2, "lines"),
          legend.title = element_text(size = 20), 
          axis.title = element_text(size = 20),
          strip.text = element_text(hjust = -0.01, size = 20),
          strip.background = element_blank()
    ) + xlab("Infectee ID") + ylab("Infector ID") + 
    guides(colour = guide_legend(override.aes = list(size=5))) +
    geom_vline(xintercept = v_lines, lty = 2, col = "black") + 
    geom_hline(yintercept = v_lines, lty = 2, col = "black") + 
    xlim(0, max(data_plot$Infectee)) + ylim(0, max(data_plot$Infector))
  
  plot_tot
  
}

# Function used in figure3 to re-organise the data plot dataframe
convert_data_plot <- function(data_plot, dt_cases, chains_nb){
  # Convert "from" and "to" to numeric values
  num_from <- as.numeric(levels(data_plot$data$from)[data_plot$data$from])
  num_to <- as.numeric(levels(data_plot$data$to)[data_plot$data$to])
  # Set importation to NA
  num_from[num_from==0] <- NA
  # Set the links between different chains to "grey"
  plot_cols <- rep("other", nrow(data_plot$data))
  same_chain <- which(!is.na(dt_cases[num_from, chain]) &
                        !is.na(dt_cases[num_to, chain]) &
                        dt_cases[num_from, chain] == dt_cases[num_to, chain])
  same_chain <- same_chain[is.element(dt_cases[num_from, chain][same_chain], 
                                      chains_nb)]
  plot_cols[same_chain] <- as.character(dt_cases[num_from, chain][same_chain])
  plot_cols <- paste0("chain ", plot_cols)
  plot_cols[plot_cols == "chain other"] <- "other chains"
  
  new_data_plot <- data_plot$data
  # Re-arrange the format / name of each column of new_data_plot
  new_data_plot$Chain <- factor(plot_cols, levels = c(paste0("chain ", chains_nb), "other chains"))
  new_data_plot$Infector <- as.numeric(levels(new_data_plot$from)[new_data_plot$from])
  new_data_plot$Infectee <- as.numeric(levels(new_data_plot$to)[new_data_plot$to])
  new_data_plot$Frequency <- new_data_plot$frequency
  return(new_data_plot)
  
}

# Plot consensus trees
figure4 <- function(out_wt, out_comb, data, burnin){
  # Define colour scheme
  cols <- c("#88CCEE","#44AA99","#AA4499","#999933","#117733",
            "#882255","#6699CC")
  # Select the cases from the epi dataset that were sequenced 
  dt_cases <- data[[1]][data[[2]]$seq != ""][order(chain, date_ons)]
  # Re-order chain from 1 to length(unique(dt_cases$chain))
  dt_cases$chain <- as.numeric(factor(dt_cases$chain))
  # Generate consensus tree of out_wt and out_comb
  graph_wt <- out_to_graph(out_wt, dt_cases, cols, burnin)
  graph_comb <- out_to_graph(out_comb, dt_cases, cols, burnin)
  # Generate plot
  par(mfrow = c(1,2), mar = c(0,0,0,0), omi = c(0,0,0,0), cex.axis = 1.5, 
      cex.lab = 2, cex = 1.5)
  set.seed(1)
  ## Epi approach
  plot.igraph(graph_wt, edge.arrow.size=0.4, edge.arrow.width=0.3,
              vertex.label.cex=0, edge.label.cex=1, vertex.size=7.5,
              vertex.label.font=2)
  mtext("A", side = 3, line = -2, adj = 0.1, cex = 2)
  set.seed(1)
  ## Combined approach
  plot.igraph(graph_comb, edge.arrow.size=0.4, edge.arrow.width=0.3,
              vertex.label.cex=0, edge.label.cex=1, vertex.size=7.5,
              vertex.label.font=2)
  mtext("B", side = 3, line = -2, adj = 0.1, cex = 2)
  par(fig = c(0, 1, 0, 1),omi = c(0,0,0,0), mar = c(0,0,0,0), new = T)
  plot.new()
  ## Add legend
  legend('center', legend = unique(dt_cases$chain), fill = c(cols), 
         border="black", bty = "n")
  
}

# Function used to compute the consensus trees (used in figure4())
out_to_graph <- function(out_o2, dt_cases, cols, burnin){
  # Compute the maximum number of importations accross the different iterations
  n_imports <- out_o2[, grep("alpha", colnames(out_o2))] %>%
    is.na %>% rowSums() %>% max
  # Extract the most likely transmission tree, and convert to a data table
  trans_tree <- data.table(summary(out_o2, burnin = burnin)$tree)
  # If there are extra imports in trans_tree (i.e. if there are more imports in 
  # trans_tree than the maximum number across iterations), add the "best indexes"
  # as secondary cases
  if(n_imports < nrow(trans_tree[is.na(from),])){
    # Find the import with weakest support
    min_freq_na <- (trans_tree[is.na(from),support] %>% sort(decreasing = T))[n_imports]
    # The new indexes will be the importations with support worse than min_freq_na
    which_new_index <- which(is.na(trans_tree$from) & 
                               trans_tree$support < min_freq_na)
    # Best index for each case
    best_index <- 
      apply(out_o2[, grep("alpha", colnames(out_o2))], 2, 
            function(X){ 
              if(all(is.na(X))) return(0)
              return(names(table(X))[which.max(table(X))])}) %>% 
      as.numeric
    # Best generation
    best_gen <- apply(out_o2[, grep("kappa", colnames(out_o2))], 2, 
                      function(X){ 
                        if(all(is.na(X))) return(0)
                        return(names(table(X))[which.max(table(X))])}) %>% 
      as.numeric
    # Support for best index
    best_prop <- (apply(out_o2[, grep("alpha", colnames(out_o2))], 2, 
                        function(X){ 
                          if(all(is.na(X))) return(0)
                          return(names(table(X))[which.max(table(X))])}) %>% 
                    as.numeric)/nrow(out_o2)
    # Replace weakest importations by the most supported index
    trans_tree[which_new_index, from := best_index[which_new_index]]
    trans_tree[which_new_index, generations := best_gen[which_new_index]]
    trans_tree[which_new_index, support := best_prop[which_new_index]]
  }
  # Test if there are any loop (i.e. i infects j which infects i)
  for(i in seq_len(nrow(trans_tree))){
    if(!is.na(trans_tree[to == i, from]) &&
       !is.na(trans_tree[to == trans_tree[to == i, from], from]) &&
       trans_tree[to == trans_tree[to == i, from], from] == i){
      # Extract i's best ancestor
      index <- trans_tree[to == i, from]
      # Find all potential infectors of "index" 
      all_ances <- table(as.numeric(out_o2[, paste0("alpha_", index)]))
      all_ances <- all_ances[!is.element(names(all_ances), 
                                         trans_tree[from == index, to])]
      # Replace the infector of j by the next best ancestor
      if(length(all_ances) > 0) 
        trans_tree[to == index, 
                   from := as.numeric(names(sort(all_ances, decreasing = T))[1])]
    }
  }
  # Add the chain of the infectee to trans_tree 
  trans_tree[, chain_final := dt_cases[to, chain]]
  trans_tree[is.na(from), from := 0]
  # Create graph object
  Ebola_graph <- graph.data.frame(trans_tree)
  
  E(Ebola_graph)$label <- trans_tree[,from]
  # Remove edges for importations
  Ebola_graph <- delete_edges(Ebola_graph, which(E(Ebola_graph)$label==0))
  
  E(Ebola_graph)$label <- ""
  Ebola_graph <- delete_vertices(Ebola_graph, "0")
  V(Ebola_graph)$label <- labels(V(Ebola_graph))
  # Add colours depending on the infectee's chain
  V(Ebola_graph)$color <-  "grey"
  V(Ebola_graph)$color <- sapply(V(Ebola_graph)$label, function(X)  {
    if(is.na(trans_tree[to==X, chain_final]))
      return(V(Ebola_graph)$color[V(Ebola_graph)$name==X])
    else if(trans_tree[to==X, chain_final] == 1) return(cols[1])
    else if(trans_tree[to==X, chain_final] == 2) return(cols[2])
    else if(trans_tree[to==X, chain_final] == 3) return(cols[3])
    else if(trans_tree[to==X, chain_final] == 4) return(cols[4])
    else if(trans_tree[to==X, chain_final] == 5) return(cols[5])
    else if(trans_tree[to==X, chain_final] == 6) return(cols[6])
    else if(trans_tree[to==X, chain_final] == 7) return(cols[7])
    else return(V(Ebola_graph)$color[V(Ebola_graph)$name==X])})
  V(Ebola_graph)$label <- ""
  return(Ebola_graph)
}

# Alpha plot by chain
figure5 <- function(out_wt, out_comb, data, burnin){
  # Define colour scheme
  cols <- c("#88CCEE","#44AA99","#AA4499","#999933","#117733",
            "#882255","#6699CC")
  # Select the cases from the epi dataset that were sequenced 
  dt_cases <- data[[1]][data[[2]]$seq != ""][order(chain, date_ons)]
  dt_cases$chain <- as.numeric(factor(dt_cases$chain))
  chains <- unique(dt_cases$chain)
  # Generate probability of transition between chains in out_comb and out_wt
  chain_prop_o2 <- alpha_to_chain(out_comb, dt_cases, burnin, chains)
  chain_prop_wt <- alpha_to_chain(out_wt, dt_cases, burnin, chains)

  # Move chain_prop_o2 and chain_prop_wt to long format 
  data_comb <- data.table(From = rep(colnames(chain_prop_o2), 
                                     ncol(chain_prop_o2)),
                          To = rep(colnames(chain_prop_o2), 
                                   each = ncol(chain_prop_o2)),
                          Frequency = c(chain_prop_o2))
  data_comb$From <- factor(data_comb$From, levels = chains)
  data_comb$To <- factor(data_comb$To, levels = chains)
  
  data_wt <- data.table(From = rep(colnames(chain_prop_wt), 
                                     ncol(chain_prop_wt)),
                          To = rep(colnames(chain_prop_wt), 
                                   each = ncol(chain_prop_wt)),
                          Frequency = c(chain_prop_wt))
  data_wt$From <- factor(data_wt$From, levels = chains)
  data_wt$To <- factor(data_wt$To, levels = chains)
  
  # Remove connections where the probability of connection is below 5%
  data_comb <- data_comb[Frequency > .05]
  data_wt <- data_wt[Frequency > .05]
  
  data_comb$Type <- "B"
  data_wt$Type <- "A"
  # Create the variable "clust" (irrelevant for simulated data, was used to 
  # highlight the connection between the different chains)
  data_wt$clust <- "Others"
  
  data_comb$clust <- "Others"
  data_comb$clust[data_comb$From == data_comb$To] <- 
    paste0("Group ", data_comb$From[data_comb$From == data_comb$To])
  # Merge data_comb and data_wt
  data_plot <- rbind(data_comb, data_wt)
  data_plot$Type <- factor(data_plot$Type, levels = c("A", "B"))
  # Generate plot
  plot_tot <- ggplot() + geom_point(data = data_plot, 
                                    aes(x = To, y = From, color = clust,
                                        size = Frequency)) + 
    facet_grid(~Type) + 
    scale_color_manual(values = c(cols, "grey")) + theme_classic()  + 
    theme(text = element_text(size = 15),
          panel.spacing = unit(2, "lines"),
          legend.title = element_text(size = 20), 
          axis.title = element_text(size = 20),
          strip.text = element_text(hjust = -0.01, size = 20),
          strip.background = element_blank()
    ) + xlab("Infectee chain") + ylab("Infector chain") + 
    guides(colour = guide_legend(override.aes = list(size=5)))
  plot_tot
  
}

# Function used to compute the connectivity between chains (used to generate the 
# chain alpha plots)
alpha_to_chain <- function(out, dt_cases, burnin = 2e4, chains){
  # Select columns containing the infectors for each case
  alpha <- out[, grep("alpha", colnames(out))]
  # Remove burnin phase
  alpha <- alpha[out$step >= burnin,]
  # Initialise matrix containing the chain of each infector
  chain_mat <- alpha * 0
  for (j in seq_len(ncol(alpha))){
    chain_mat[,j] <- dt_cases[alpha[,j], chain]
  }
  # Initialise matrix indicating what chain each case is connected to,
  chain_tab <- matrix(nrow = length(chains), ncol = ncol(alpha))
  rownames(chain_tab) <- chains
  for (j in seq_len(ncol(alpha))){
    # Convert the chains to factor 
    vals <- factor(chain_mat[,j], levels = chains)
    # Proportion of each chain
    chain_tab[,j] <- table(vals) / sum(!is.na(alpha[,j]))
    if(all(is.na(alpha[,j]))) chain_tab[,j] <- 0
  }
  chain_tab <- chain_tab[, !is.na(chains)]
  # Aggregate chain_tab using the infectee's chain
  chain_prop <- 
    aggregate(t(chain_tab), by = list(dt_cases$chain[!is.na(dt_cases$chain)]), mean)[, -1]
  rownames(chain_prop) <- colnames(chain_prop)
  return(chain_prop %>% t)
}

# Plot the match between the data and inferred trees
figure6 <- function(out_wt, out_comb, data, burnin){
  # Define colour scheme
  cols <- c("#88CCEE","#44AA99","#AA4499","#999933","#117733",
            "#882255","#6699CC")
  # Select the cases from the epi dataset that were sequenced 
  dt_cases <- data[[1]][data[[2]]$seq != ""][order(chain, date_ons)]
  dt_cases$chain <- as.numeric(factor(dt_cases$chain))
  chains <- unique(dt_cases$chain)
  ## Find the closest sequened infector for each case
  dt_cases[, infector_ind := infector]
  # Extract the epi dataset containing all cases
  dt_cases_all <- data[[1]][order(chain, date_ons)]
  vec_infector <- dt_cases_all$infector
  names(vec_infector) <- dt_cases_all$ID
  # While at least one entry of dt_cases$infector_ind that is not included in dt_cases$ID
  # Set infector_ind to the infector of infector_ind using vec_infector.
  while(any(!is.element(paste0("case_", dt_cases$infector_ind), dt_cases$ID) & 
            !is.na(dt_cases$infector_ind))){
    cases_indirect <- 
      which(!is.element(paste0("case_", dt_cases$infector_ind), dt_cases$ID) & 
              !is.na(dt_cases$infector_ind))
    dt_cases[cases_indirect, 
             infector_ind := vec_infector[paste0("case_", infector_ind)]]
  }
  ## For each case, compute the proportion of iterations where the case is connected
  ## to the right infector
  case_comb <- inference_case(out_comb, dt_cases, burnin)
  case_wt <- inference_case(out_wt, dt_cases, burnin)
  
  fact_chains <- factor(dt_cases$chain)
  ## Plot the proportion of iterations where cases are linked to the right
  ## infector. 
  par(mfrow = c(1,2), oma = c(4,4,0,0), mar = c(0,0,1,1) + 1.3, bty = "l",
      cex.axis = 1.5, cex.lab = 2)
  ## In epi model
  plot(case_wt, type = "h", ylim = c(0,1), pch = 19, cex = 1.2, lwd = 8,
       lend = "butt", col = cols[fact_chains[as.numeric(names(case_wt))]])
  legend('top', legend = paste0("chain ", 1:8), cex = 1.5,
         fill = c(cols), border="black", bty = "n", ncol = 2)
  mtext("A", side = 3, line = 0, adj = 0.1, cex = 2)
  ## In combined model
  plot(case_comb, type = "h", ylim = c(0,1), pch = 19, cex = 1.2, lwd = 8,
       lend = "butt", col = cols[fact_chains[as.numeric(names(case_comb))]])
  mtext("B", side = 3, line = 0, adj = 0.1, cex = 2)
  title(xlab = "Case", outer = T, line = 2)
  title(ylab = "Proportion", outer = T, line = 2)
}

# Function used to compute the proportion of iteration where cases are 
# linked to the right infector
inference_case <- function(out, dt_cases, burnin){
  # Extract infector for each iteration
  alpha <- out[, grep("alpha", colnames(out))]
  alpha <- alpha[out$step >= burnin,]
  # Extract (indirect) infector for each case
  vec_index <- dt_cases$infector_ind
  vec_index[!is.element(paste0("case_", vec_index), dt_cases$ID)] <- NA
  # Initialise the proportion of iteration where each case is connected to the right 
  # infector
  prop <- numeric(sum(!is.na(vec_index)))
  names(prop) <- which(!is.na(vec_index))
  for (j in seq_along(which(!is.na(vec_index)))){
    # Extract the position j corresponds to in vec_index
    i <- which(!is.na(vec_index))[j]
    # Extract all infectors of i in the matrix alpha
    infectors <- alpha[,i]
    # Identify the actual infector of i in vec_index
    case_i <- which(dt_cases$ID == paste0("case_", vec_index[i]))
    # Compute the number of iterations (and the proportion) where the right 
    # infector was identified by the model
    sum_i <- (infectors == case_i) %>% sum(na.rm = T)
    prop[j] <- sum_i / length(infectors[!is.na(infectors)])
  }
  # Remove imported cases
  prop <- prop[is.finite(prop)]
  return(prop)
}

# Plot the posterior trace in each model
supp_fig2 <- function(out_wt, out_comb, burnin){
  par(mfrow = c(2,2), oma = c(4,4,0,0), mar = c(0,0,1,1) + 1.3, bty = "l",
      cex.axis = 1.5, cex.lab = 2)
  plot(out_wt$post ~ out_wt$step, type = "l", 
       main = "\nEpidemiological approach")
  plot(out_wt[out_wt$step > burnin,]$post ~ 
         out_wt[out_wt$step > burnin,]$step, type = "l", 
       main = paste0("\nEpidemiological approach, burnin = ", burnin))
  plot(out_comb$post ~ out_comb$step, type = "l", 
       main = "\nCombined approach")
  plot(out_comb[out_comb$step > burnin,]$post ~ 
         out_comb[out_comb$step > burnin,]$step, type = "l", 
       main = paste0("\nCombined approach, burnin = ", burnin))
  title(xlab = "Iterations", ylab = "Posterior", outer = T, line = 2)
  
}