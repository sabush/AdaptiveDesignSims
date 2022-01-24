library(magick)
library(tidyverse)
library(gganimate)
library(cbastylr)
library(plotly)

niter <- 20000
probs <- c(0.2, 0.25, 0.26, 0.3)
exp_prop <- 0.1
status_table <- data.frame(iteraton = numeric(), Arm = numeric(), n_a = numeric(), n_b = numeric(),
                           n_c = numeric(), n_d = numeric(), s_a = numeric(), s_b = numeric(),
                           s_c = numeric(), s_d = numeric(), r_a = numeric(), r_b = numeric(),
                           r_c = numeric(), r_d = numeric())

# Initialise with 10 random samples
n_counter <- rep(0, length(probs))
s_counter <- rep(0, length(probs))

for(i in 1:niter){
  UCB <- s_counter/n_counter + sqrt(2 * log(i) / n_counter)
  UCB[is.nan(UCB)] <- Inf
  selected_arm <- which.max(UCB)
  n_counter[selected_arm] <- n_counter[selected_arm] + 1
  s_counter[selected_arm] <- s_counter[selected_arm] + rbinom(1, 1, probs[selected_arm])
  work_output <- c(i, selected_arm, n_counter, s_counter, s_counter/n_counter)
  status_table <- rbind(status_table, work_output)
}

colnames(status_table) <- c('iteration', 'Arm', 'n_a', 'n_b', 'n_c', 'n_d', 's_a', 's_b',
                            's_c', 's_d', 'r_a', 'r_b', 'r_c', 'r_d')


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

pdf(paste('./Output/UCB1.pdf'), width = 12, height = 7)

for(iter in iterationplotlist){
  plot_greedyr50 <- status_unstacked %>% 
    filter(iteration == iter) %>%
    mutate(iteration = as.character(iteration)) %>% 
    ggplot(aes(x = fct_rev(group), y = p_est, fill = fct_rev(group))) + geom_bar(stat = 'identity') + 
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
                    top = text_grob("Upper Confidence Bound: Sample Sizes and Respose Estimates",
                                    face = "bold", size = 18),
                    bottom = text_grob("Response rates - Arm A = 20%, Arm B = 25%, Arm C = 26%, Arm D = 30%", 
                                       face = "italic", size = 14),
    )
  print(combined) 
}
dev.off()


# # Generate pdfs for plotting
# 
# iterationplotlist <- c(5:20, seq(20, 100, 5), seq(200, 1000, 50), seq(1200, niter, 200))
# 
# 
# 
# 
# 
# plot_UCBn <- status_table %>%
#   
#   dplyr::select(iteration, n_a, n_b, n_c, n_d) %>%
#   
#   gather(group, sample_size, -iteration) %>%
#   
#   mutate(group = gsub("n_", "", group)) %>%
#   
#   left_join(status_table %>%
#               
#               dplyr::select(iteration, n_a, n_b, n_c, n_d) %>%
#               
#               gather(group, sample_size, -iteration) %>%
#               
#               mutate(group = gsub("n_", "", group)) %>%
#               
#               group_by(iteration) %>%
#               
#               summarise(total = sum(sample_size))) %>%
#   
#   mutate(prop = sample_size / total) %>%
#   
#   filter(iteration < 10) %>%
#   
#   ggplot(aes(x = group, y = prop)) +
#   
#   geom_bar(stat = 'identity', fill = cbastylr_corp()[1]) + theme_classic(base_size = 16) +
#   
#   scale_y_continuous(labels = scales::percent_format()) +
#   
#   labs(title = "Distribution of Sample Sizes",
#        
#        caption = "Response rates - Arm A = 20%, Arm B = 25%, Arm C = 26%, Arm D = 30%",
#        
#        # subtitle = "Total Sample Size: {frame_time}",
#        
#        x = "Group", y = "Proportion of Sample") +
#   
#   transition_states(
#     
#     states = iteration,
#     
#     transition_length = 0.1,
#     
#     state_length = 0.1
#     
#   ) +
#   
#   enter_fade() +
#   
#   exit_shrink()
# 
# 
# 
# 
# 
# UCBn_gif <- animate(plot_UCBn, width = 600, height = 400, renderer = magick_renderer())
# 
# UCBn_gif
# 
# 
# 
# 
# 
# plot_UCBr <- status_table %>%
#   
#   dplyr::select(iteration, r_a, r_b, r_c, r_d) %>%
#   
#   gather(group, response_rate, -iteration) %>%
#   
#   mutate(group = gsub("n_", "", group)) %>%
#   
#   ggplot(aes(x = group, y = response_rate)) +
#   
#   geom_bar(stat = 'identity', fill = cbastylr_corp()[1]) + theme_classic(base_size = 16) +
#   
#   scale_y_continuous(labels = scales::percent_format()) +
#   
#   labs(title = "Distribution of Response Rates",
#        
#        # subtitle = "Total Sample Size: {frame_time}",
#        
#        caption = "Response Rates - Arm A = 20%, Arm B = 25%, Arm C = 26%, Arm D = 30%",
#        
#        x = "Group", y = "Response Rate") +
#   
#   transition_states(
#     
#     states = iteration,
#     
#     transition_length = 0.1,
#     
#     state_length = 0.1
#     
#   ) +
#   
#   enter_fade() +
#   
#   exit_shrink()
# 
# 
# 
# UCBr_gif <- animate(plot_UCBr, width = 600, height = 400, renderer = magick_renderer())
# 
# UCBr_gif
# 
# 
# 
# plot_UCBss <- data.frame(iteration = 1:niter, x = c('a'), y = 1:niter) %>%
#   
#   ggplot(aes(x = x, y = y)) + geom_bar(stat = 'identity', fill = cbastylr_corp()[1]) +
#   
#   coord_flip() + theme_classic(base_size = 16) + xlim(0, niter) +
#   
#   theme(axis.text.y = element_blank(), axis.ticks.x = element_blank()) +
#   
#   labs(title = "Total Sample Size",
#        
#        y = "Sample Size") +
#   
#   transition_states(
#     
#     states = iteration,
#     
#     transition_length = 0.1,
#     
#     state_length = 0.1
#     
#   ) +
#   
#   enter_fade() +
#   
#   exit_shrink()
# 
# 
# 
# UCBss_gif <- animate(plot_UCBss, width = 600, height = 400, renderer = magick_renderer())
# 
# UCBss_gif
# 
# 
# 
# 
# 
# UCB_gif <- image_append(c(UCBss_gif[1], UCBn_gif[1], UCBr_gif[1]))
# 
# for(i in 2:niter){
#   
#   combined <- image_append(c(UCBss_gif[i], UCBn_gif[i], UCBr_gif[i]))
#   
#   UCB_gif <- c(UCB_gif, combined)
#   
# }
# 
# 
# 
# UCB_gif
# 
# image_write(UCB_gif, "UCB.gif")