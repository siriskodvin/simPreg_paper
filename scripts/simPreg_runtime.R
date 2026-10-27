
##---------------------------------------------------------------##
#                                                                 #
#  Script for benchmarking the runtime of the simPreg functions.  #
#  Runtime will vary depending on the machine specifications.     #
#                                                                 #
##---------------------------------------------------------------##

library(simPreg)
library(ggplot2)
library(dplyr)

# Initiate results list
res_l <- list()

# Number of repetitions
j_max <- 1000

# Loop on varying sample size
for(i in 9:1){
  
  # Print status
  print(paste0("i = ", i))

  # Loop repetitions
  for(j in 1:j_max){
    
    # Time proportions
    t <- Sys.time()
    df_prop <- simPregProp()
    t <- as.numeric(Sys.time() - t)
    
    # Time sampling
    v <- Sys.time()
    df_samp <- simPregSamp(df_prop, n = 10^i, expand = FALSE)
    v <- as.numeric(Sys.time() - v)
    
    # Save resuls
    res_l[[length(res_l)+1]] <- data.frame(i = i,
                                           t = t,
                                           v = v)
  }
}

# Collapse result list into dataframe
df_res <- do.call(rbind, res_l)

# Adjust to plot
df_pl <- df_res
df_pl$tv <- df_pl$t + df_pl$v

# Aggregate results
df_pl_tab <- aggregate(cbind(t, v, tv) ~ i,
                       data = df_pl,
                       FUN = function(x) c(mean = mean(x), sd = sd(x)))
df_pl_tab <- data.frame(
  df_pl_tab[c("i")],
  setNames(
    as.data.frame(do.call(cbind, df_pl_tab[c("t", "v", "tv")])),
    c("t.mean", "t.sd", "v.mean", "v.sd", "tv.mean", "tv.sd")))

# Create dataset for stacking
df_pl_tab_stc <- rbind(
  data.frame(
    df_pl_tab[c("i")],
    lower = 0,
    upper   = df_pl_tab$t.mean,
    func = "prop"),
  data.frame(
    df_pl_tab[c("i")],
    lower = df_pl_tab$t.mean,
    upper   = df_pl_tab$tv.mean,
    func = "samp"))

# Plot t and v together for expand = FALSE
df_pl_sub <- df_pl_tab_stc
pl <- ggplot(df_pl_sub, aes(x = i)) +#, color = hr_fac, shape = hr_fac)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(func)), alpha = 0.65) +
  scale_x_continuous(
    breaks = 1:9,
    labels = function(x) parse(text = paste0("10^", x))) +
  scale_fill_manual(values = c("prop" = "#6796A7",
                               "samp" = "#824F73"),
                    breaks = c("samp", "prop"),
                    labels = c("prop" = "simPregProp",
                               "samp" = "simPregSamp"),
                    guide = guide_legend(order = 1)) +
  labs(title = NULL, x = "Sample size", y = "Mean runtime (seconds)", fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.border =  element_rect(color = "black", fill = NA, linewidth = 0.5))


