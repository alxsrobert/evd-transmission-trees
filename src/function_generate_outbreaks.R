generate_outbreak <- function(n_import, date_range, max_cases, seq_ref, prop_var,
                              w_dens, f_dens, mut_rate, r0, prob_seq_rep){
  # Initialise data frame with sequences
  df_seq <- data.frame()
  # Draw which elements of the genetic sequence can be changed
  which_change <- sample(x = seq_along(seq_ref), 
                         size = round(prop_var * length(seq_ref)))
  # Initialise count of the number of cases, and the number of importations
  count <- 1
  n_cases <- 0
  # Loop over the maximum number of importations
  for(i in seq_len(n_import)){
    # Draw the infection date
    import_date_inf <- sample(size = 1, x = seq(date_range[1], date_range[2], "days"))
    # Use the incubation period to draw the onset date
    import_date_ons <- import_date_inf + 
      sample(size = 1, x = seq_along(f_dens), prob = f_dens)
    # Draw changes to the genetic sequence
    import_seq <- seq_ref
    import_seq[which_change] <- sample(x = c("a", "t", "c", "g"), 
                                       size = length(which_change), replace = T)
    # If first case, initialise dt_cases (epi information), and df_seq (genetic sequences)
    if(count == 1){
      dt_cases <- data.frame(ID = "case_1", date_inf = import_date_inf, 
                             date_ons = import_date_ons, infector = NA, chain = i) 
      df_seq <- data.frame(ID = paste0("case_", count), 
                           seq = toupper(paste0(import_seq, collapse = "")))
    }else{
      # Otherwise, add the new importation to the epi and sequence datasets
      dt_cases <- rbind.data.frame(dt_cases, c(paste0("case_", count),
                                               as.character(import_date_inf), 
                                               as.character(import_date_ons), NA, i))
      df_seq <- rbind.data.frame(df_seq, c(paste0("case_", count),
                                           toupper(paste0(import_seq, collapse = ""))))
    }
    # Increment n_cases and count
    n_cases <- n_cases + 1
    count <- count + 1
    # Stop the loop if the number of cases is above max_cases
    if(nrow(dt_cases) < max_cases){
      j <- n_cases
      # For each case in the chain (stop chain if there are more than 40 cases)
      while(j <= n_cases & sum(dt_cases$chain == i) < 40){
        # Extract the infection date and sequence of j
        date_inf <- dt_cases$date_inf[j]
        seq_inf <- as.character(df_seq[j,]$seq)
        # Draw the number of secondary cases using a negative binomial distribution
        n_sec <- rnbinom(1, mu = r0, size = .7)
        # Add n_sec to n_cases
        n_cases <- n_cases + n_sec
        
        # Generate the characteristics of the new cases:
        if(n_sec > 0){
          for(k in seq_len(n_sec)){
            # Use the function generate_case to get the onset date, infection
            # date and genetic sequence of k
            new_case <- 
              generate_case(w_dens, f_dens, mut_rate, date_inf, seq_inf, which_change)
            row_new_case <- c(paste0("case_", count),
                              as.character(new_case$date[1]), 
                              as.character(new_case$date[2]), j, i)
            # Add the new characteristics to dt_cases
            dt_cases <- rbind.data.frame(dt_cases, row_new_case)
            # Add the new sequence to df_seq
            df_seq <- 
              rbind.data.frame(df_seq, c(paste0("case_", count),
                                         toupper(paste0(new_case$sequence, 
                                                        collapse = ""))))
            # Update case count
            count <- count + 1
          }
        }
        j <- j + 1
      }
    }
  }
  # Draw the number of reported sequences using prob_seq_rep
  nb_reported_seq <- round(nrow(df_seq) * prob_seq_rep)
  keep_seq <- sample(seq_len(nrow(df_seq)), size = nb_reported_seq)
  # Remove unreported sequences
  df_seq[-keep_seq, "seq"] <- ""
  # Return the simulated epi data and genetic sequences
  return(list(dt = as.data.table(dt_cases), seq = df_seq))
}


generate_case <- function(w_dens, f_dens, mut_rate, date_inf, seq_inf, which_change){
  # Initialise the new sequence as a chain of characters
  new_seq <- strsplit(seq_inf, "")[[1]]
  # Draw the infection date using the serial interval and the infector's infection date
  inf_date <- date_inf + sample(size = 1, x = seq_along(w_dens), prob = w_dens)
  # Draw the onset date
  ons_date <- inf_date + sample(size = 1, x = seq_along(f_dens), prob = f_dens)
  # Draw the number of mutation, using the mutation rate and the serial interval
  n_mut <- rpois(1, mut_rate * as.numeric(inf_date - date_inf))
  # Draw the positions that will mutate
  changes <- sample(x = which_change, size = n_mut, replace = F)
  new_seq[changes] <- sample(x = c("a", "t", "c", "g"), 
                                size = length(changes), replace = T)
  
  return(list(date = c(inf_date, ons_date), 
              sequence = paste0(new_seq, collapse = "")))
}
