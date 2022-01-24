library(tidyverse)
library(gganimate)
library(cbastylr)
library(forcats)
library(ggpubr)

niter <- 20000
probs <- c(0.2, 0.25, 0.26, 0.3)
exp_prop <- 0.1
status_table <- data.frame(iteraton = numeric(), Arm = numeric(), n_a = numeric(), 
                           n_b = numeric(), n_c = numeric(), n_d = numeric(), 
                           s_a = numeric(), s_b = numeric(), s_c = numeric(), 
                           s_d = numeric(), r_a = numeric(), r_b = numeric(),
                           r_c = numeric(), r_d = numeric())

# Initialise with 20 samples, 5 from each group
n_counter <- rep(0, length(probs))
s_counter <- rep(0, length(probs))

for(i in 1:20){
  selected_arm <- (1:20 %% 4 + 1)[i]
  n_counter[selected_arm] <- n_counter[selected_arm] + 1
  s_counter[selected_arm] <- s_counter[selected_arm] + rbinom(1, 1, probs[selected_arm])
  work_output <- c(i, selected_arm, n_counter, s_counter, s_counter/n_counter)
  status_table <- rbind(status_table, work_output)
}

colnames(status_table) <- c('iteration', 'Arm', 'n_a', 'n_b', 'n_c', 'n_d', 's_a', 's_b',
                            's_c', 's_d', 'r_a', 'r_b', 'r_c', 'r_d')

for(i in 21:niter){
  best <- which.max(s_counter/n_counter)
  selected_arm <- ifelse(runif(1) <= exp_prop,
                         (matrix(c(1,2,3,4), nrow = 1) %*% rmultinom(1, 1, c(1,1,1,1)))[[1,1]],
                         best)
  n_counter[selected_arm] <- n_counter[selected_arm] + 1
  s_counter[selected_arm] <- s_counter[selected_arm] + rbinom(1, 1, probs[selected_arm])
  work_output <- c(i, selected_arm, n_counter, s_counter, s_counter/n_counter)
  status_table <- rbind(status_table, work_output)
}

## Calculate estimates and error bars for each arm at each iteration
status_unstacked <- status_table %>% 
  select(iteration, n_a:n_d) %>% 
  gather(group, sample_size, -iteration) %>% 
  mutate(group = gsub('n_','',group)) %>% 
  inner_join(status_table %>% 
               select(iteration, r_a:r_d) %>% 
               gather(group, p_est, -iteration) %>% 
               mutate(group = gsub('r_','',group))) %>% 
  rowwise() %>% 
  mutate(lci = max(p_est - 1.96 * sqrt(p_est * (1 - p_est) / sample_size),0.),
         uci = min(p_est + 1.96 * sqrt(p_est * (1 - p_est) / sample_size),1.),
         group = factor(group, levels = c('a', 'b', 'c', 'd'))
  )

iterationplotlist <- c(seq(20, 100, 5), seq(200, 1150, 50), seq(1200, niter, 200))

pdf(paste('./Output/EpsGreedy.pdf'), width = 12, height = 7)

for(iter in iterationplotlist){
  plot_greedyr50 <- status_unstacked %>% 
    filter(iteration == iter) %>%
    mutate(iteration = as.character(iteration)) %>% 
    ggplot(aes(x = forcats::fct_rev(group), y = p_est, fill = forcats::fct_rev(group))) + 
    geom_bar(stat = 'identity') + 
    geom_errorbar(aes(ymin = lci, ymax = uci)) + theme_classic(base_size = 14) +
    coord_flip() + scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    scale_fill_discrete(guide = guide_legend(reverse=TRUE)) +
    labs(title = "Distribution of Response Rates",
         subtitle = paste0("Total Sample Size = ", iter),
         x = "Group", y = "Estimated Response Rate", fill = "Group")
  
  plot_greedyn50 <- status_unstacked %>% 
    filter(iteration == iter) %>%
    mutate(prop = sample_size / iteration,
           iteration = as.character(iteration)) %>% 
    ggplot(aes(x = 1, y = prop, fill = fct_rev(group))) + 
    geom_bar(stat = 'identity') + theme_classic(base_size = 14) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_x_continuous(breaks = NULL) + coord_flip() +
    scale_fill_discrete(guide = guide_legend(reverse=TRUE)) +
    labs(title = "Distribution of Sample Sizes",
         x = " ", y = "Proportion of Sample")
  
  combined <- ggarrange(plot_greedyr50, plot_greedyn50,
                        ncol = 1, nrow = 2, heights = c(4, 1), 
                        legend = 'bottom', common.legend = T)
  
  combined <- 
    annotate_figure(combined,
                    top = text_grob("Epsilon Greedy: Sample Sizes and Respose Estimates",
                                    face = "bold", size = 18),
                    bottom = text_grob("Response rates - Arm A = 20%, Arm B = 25%, Arm C = 26%, Arm D = 30%", 
                                       face = "italic", size = 14),
    )
  print(combined) 
}
dev.off()

