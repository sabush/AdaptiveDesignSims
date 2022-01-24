library(gganimate)
library(gifski)
library(tidyverse)
library(magick)
library(transformr)

##### Equal Sampling - pvalues -------------------------------------------------

rr1 <- 0.3
rr2 <- 0.3
niter <- 100
nsims <- 100

cumulative_equal <-
  tibble(
    simulation = numeric(), iteration = numeric(), T1N = numeric(),
    T1S = numeric(), T2N = numeric(), T2S = numeric()
  )

for (cycle in 1:nsims) {
  T1_n <- 0
  T1_succ <- 0
  T2_n <- 0
  T2_succ <- 0
  for (i in 1:niter) {
    T1_succ <- T1_succ + ifelse(rbernoulli(1, rr1), 1, 0)
    T2_succ <- T2_succ + ifelse(rbernoulli(1, rr2), 1, 0)
    cumulative_equal <- cumulative_equal %>%
      add_row(
        simulation = cycle, iteration = i, T1N = i, 
        T1S = T1_succ, T2N = i, T2S = T2_succ
      )
  }
}

cumulative_equal <- cumulative_equal %>%
  rowwise() %>%
  mutate(
    pval = prop.test(c(T1S, T2S), c(T1N, T2N))$p.value,
    pval = if_else(!is.na(pval), pval, 1),
    Simulation = factor(simulation)
  ) %>%
  group_by(simulation) %>%
  mutate(numsig = cumsum(pval <= 0.05),
         numsig2 = cumsum(numsig)) %>%
  filter(numsig2 <= 1)

write.csv(cumulative_equal, 'D:/Dropbox/Dropbox/RProjects/animations/Output/pvalue_evolution_equal.csv')
cumulative_equal <- read_csv('D:/Dropbox/Dropbox/RProjects/animations/Output/pvalue_evolution_equal.csv') %>%
  select(-X1) %>%
  mutate(numPos = cumsum(numsig),
         percPos = round(numPos / if_else(Simulation == 1, 1, Simulation - 1) * 100, digits = 0)) %>%
  mutate(plot_time = row_number()) %>%
  mutate(Simulation = as.character(Simulation))

add_point <- cumulative_equal %>%
  filter(numsig == 1 | iteration == 100) %>%
  dplyr::select(plot_time, iteration, pval, Simulation) %>% 
  mutate(Simulation = as.numeric(Simulation))

old_plot <- cumulative_equal %>%
  rename(iteration2 = iteration)

add_point2 <- add_point %>%
  mutate(iteration2 = iteration) %>%
  select(-iteration)

plotlist <- cumulative_equal %>%
  filter(Simulation == 1) %>%
  ggplot(aes(x = iteration, y = pval, colour = Simulation, group = Simulation)) +
  geom_line(stat = 'identity', size = 2, colour = 'blue') + 
  theme_classic(base_size = 14) +
  theme(legend.position = "none") + ylim(0, 1) + xlim(0, 100) +
  labs(
    title = "Evolution of p-values: 0% Experiments stopped early (0 simulations)",
    subtitle = "Two treatments, both with a 30% chance of success (no difference)",
    x = "Sample Size (each group)", y = "p-value") +
  geom_hline(yintercept = 0.05, color = "green", size = 1.5) +
  transition_reveal(iteration)

giflist <-
  animate(
    plotlist, width = 600, height = 300, renderer = magick_renderer(), fps = 100,
    nframes = max(3,floor(cumulative_equal %>% filter(Simulation == 1) %>% nrow / 20))
  )

for(iter in 2:10) {
  print(iter)
  workdata <- cumulative_equal %>% 
    filter(Simulation == iter)
  
  perc <- workdata %>% 
    top_n(1, iteration) %>% 
    select(percPos) %>% 
    unlist()
  
  plotlist_app <- old_plot %>%
    filter(Simulation %in% 1:(iter-1)) %>%
    ggplot(aes(x = iteration2, y = pval, group = Simulation)) +
    geom_line(size = 0.5, color = 'gray90') +
    geom_line(
      data = cumulative_equal %>% filter(Simulation == iter),
      aes( x = iteration, y = pval, colour = Simulation, group = Simulation ),
      stat = 'identity', size = 2, colour = 'blue') + 
    theme_classic(base_size = 14) + theme(legend.position = "none") +
    geom_point(
      data = add_point2 %>% filter(Simulation < iter),
      aes(x = iteration2, y = pval, group = Simulation),
      colour = 'red', size = 3
    ) +
    ylim(0, 1) + xlim(0, 100) +
    labs(
      title = paste0("Evolution of p-values: ",perc,"% of Experiments stopped early (",
                     iter - 1," simulations)"),
      subtitle = "Two treatments, both with a 30% chance of success (no difference)",
      x = "Sample Size (each group)", y = "p-value"
    ) +
    geom_hline(yintercept = 0.05, color = "green", size = 1.5) +
    transition_reveal(iteration)
  
  gif_new <-
    animate(
      plotlist_app, width = 600, height = 300, renderer = magick_renderer(), fps = 100,
      nframes = max(3,floor(cumulative_equal %>% filter(Simulation == iter) %>% nrow / 20))
    )
  giflist <- c(giflist, gif_new)
}

for(iter in 11:nsims) {
  print(iter)
  workdata <- cumulative_equal %>% 
    filter(Simulation == iter)
  
  perc <- workdata %>% 
    top_n(1, iteration-1) %>% 
    select(percPos) %>% 
    unlist()
  
  plotlist_app <- old_plot %>%
    filter(Simulation %in% 1:(iter-1)) %>%
    ggplot(aes(x = iteration2, y = pval, group = Simulation)) +
    geom_line(size = 0.5, color = 'gray90') +
    geom_line(
      data = cumulative_equal %>% filter(Simulation == iter),
      aes( x = iteration, y = pval, colour = Simulation, group = Simulation ),
      stat = 'identity', size = 2, colour = 'blue') + 
    theme_classic(base_size = 14) + theme(legend.position = "none") +
    geom_point(
      data = add_point2 %>% filter(Simulation < iter),
      aes(x = iteration2, y = pval, group = Simulation),
      colour = 'red', size = 3
    ) +
    ylim(0, 1) + xlim(0, 100) +
    labs(
      title = paste0("Evolution of p-values: ",perc,"% of Experiments stopped early (",
                     iter - 1," simulations)"),
      subtitle = "Two treatments, both with a 30% chance of success (no difference)",
      x = "Sample Size (each group)", y = "p-value"
    ) +
    geom_hline(yintercept = 0.05, color = "green", size = 1.5) +
    transition_reveal(iteration)
  
  gif_new <-
    animate(
      plotlist_app, width = 600, height = 300, renderer = magick_renderer(), fps = 100,
      nframes = 2)
  giflist <- c(giflist, gif_new)
}

for(i in 1:30) {
  print(iter)
  workdata <- cumulative_equal %>% 
    filter(Simulation == iter)
  
  perc <- workdata %>% 
    top_n(1, iteration) %>% 
    select(percPos) %>% 
    unlist()
  
  plotlist_app <- old_plot %>%
    filter(Simulation %in% 1:(iter-1)) %>%
    ggplot(aes(x = iteration2, y = pval, group = Simulation)) +
    geom_line(size = 0.5, color = 'gray90') +
    geom_line(
      data = cumulative_equal %>% filter(Simulation == iter),
      aes( x = iteration, y = pval, colour = Simulation, group = Simulation ),
      stat = 'identity', size = 2, colour = 'blue') + 
    theme_classic(base_size = 14) + theme(legend.position = "none") +
    geom_point(
      data = add_point2 %>% filter(Simulation < iter),
      aes(x = iteration2, y = pval, group = Simulation),
      colour = 'red', size = 3
    ) +
    ylim(0, 1) + xlim(0, 100) +
    labs(
      title = paste0("Evolution of p-values: ",perc,"% of Experiments stopped early (",
                     iter," simulations)"),
      subtitle = "Two treatments, both with a 30% chance of success (no difference)",
      x = "Sample Size (each group)", y = "p-value"
    ) +
    geom_hline(yintercept = 0.05, color = "green", size = 1.5) +
    transition_null()
  
  gif_new <-
    animate(
      plotlist_app, width = 600, height = 300, renderer = magick_renderer(), fps = 100,
      nframes = 2)
  giflist <- c(giflist, gif_new)
}

# giflist
image_write(giflist, "D:/Dropbox/Dropbox/RProjects/animations/Output/pvalues.gif")
