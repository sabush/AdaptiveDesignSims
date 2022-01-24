##### Adaptive Distributions ---------------------------------------------------

library(gganimate)
library(gifski)
library(tidyverse)
library(magick)
library(transformr)
library(plotly)

# 1 Two arms -------------------------------------------------------------------
## 1.1 Run Simulations ---------------------------------------------------------

niter <- 1000
probs <- c(0.2, 0.3)
alist <- c(1,1)
blist <- c(1,1)
result_table <- tibble(iteration = numeric(), Arm = numeric(), aval = numeric(), 
                       bval = numeric())

for(i in 1:niter){
  bestb = NA
  maxsample = -1
  allsamples = list()
  for(x in 1:length(probs)){
    sample_b = rbeta(1, alist[x], blist[x])
    allsamples = append(allsamples, sample_b)
    
    if(sample_b > maxsample){
      maxsample = sample_b
      bestb = x
    }
  }
  
  # Pull the arm
  outcome <- if_else(rbernoulli(1, probs[bestb]), 1, 0)
  
  # Update priors
  alist[bestb] <- alist[bestb] + outcome
  blist[bestb] <- blist[bestb] + 1 - outcome
  updaterows <- data.frame(matrix(c(rep(i, length(probs)), 1:length(probs), 
                                    alist, blist), nrow=length(probs)))
  colnames(updaterows) <- c("iteration", "Arm", "aval", "bval")
  result_table <- result_table %>% rbind(updaterows)
}

# 1.2 Generate pdfs for plotting -------------------------------------------------

iterationplotlist <- c(1:10, seq(20, 100, 5), seq(200, niter, 50))
# iterationplotlist <- 1:niter
result_table_plot <- expand.grid(iterationplotlist, 1:length(probs), x = seq(0, 1, 0.001))
colnames(result_table_plot) <-  c("iteration", "Arm", "xvalue")

result_table_plot <- result_table_plot %>%
  inner_join(result_table) %>%
  mutate(density = dbeta(xvalue, aval, bval),
         Arm = factor(Arm))

plot_mab <- result_table_plot %>%
  ggplot(aes(x = xvalue, y = density, group = Arm, colour = Arm)) +
  geom_line(aes(frame = iteration)) + theme_classic(base_size = 16) +
  # scale_colour_manual(values = cbastylr_corp()) +
  scale_x_continuous(breaks = seq(0,1,0.1), limits = c(0, 1)) +
  labs(title = "Multi Armed Bandit - Distribution of response rates",
       subtitle = "Response rates - Arm 1 = 20%, Arm 2 = 30%",
       x = "Response Rate", y = "Density")

plot_mab <- result_table_plot %>%
  ggplot(aes(x = xvalue, y = density, group = Arm, colour = Arm)) +
  geom_line(aes(frame = iteration)) + theme_classic(base_size = 16) +
  # scale_colour_manual(values = cbastylr_corp()) +
  scale_x_continuous(breaks = seq(0,1,0.1), limits = c(0, 1)) +
  labs(title = "Multi Armed Bandit - Distribution of response rates",
       subtitle = "Response rates - Arm 1 = 20%, Arm 2 = 30%",
       x = "Response Rate", y = "Density")

plotly_mab <- ggplotly(plot_mab)
htmlwidgets::saveWidget(plotly_mab, file = "plotly_mab.html")

result_table_ss <- result_table %>%
  filter(iteration %in% iterationplotlist) %>%
  mutate(ss = aval + bval - 2) %>%
  select(-aval, -bval) %>%
  accumulate_by(~iteration)

plot_mab_ss <- result_table_ss %>%
  # mutate(Arm = factor(Arm)) %>%
  ggplot(aes(x = iteration, y = ss, group = Arm, colour = Arm)) +
  geom_line(aes(frame = iteration)) +
  theme_classic(base_size = 16) +
  # scale_colour_manual(values = cbastylr_corp()) +
  # scale_x_log10() +
  labs(title = "Multi Armed Bandit - Cumulative Sample Size",
       subtitle = "Response rates - Arm 1 = 20%, Arm 2 = 30%",
       x = "Iteration", y = "Sample Size")

plotly_mab_ss <- ggplotly(plot_mab_ss)

htmlwidgets::saveWidget(plotly_mab_ss, file = "plotly_mab_ss.html")

# Try four arms ----------------------------------------------------------------
# Run Simulations

niter <- 5000
probs <- c(0.2, 0.25, 0.26, 0.3)
alist <- rep(1, length(probs))
blist <- rep(1, length(probs))
result_table2 <- tibble(iteration = numeric(), Arm = numeric(), 
                        aval = numeric(), bval = numeric())

for(i in 1:niter){
  bestb = NA
  maxsample = -1
  allsamples = list()
  for(x in 1:length(probs)){
    sample_b = rbeta(1, alist[x], blist[x])
    allsamples = append(allsamples, sample_b)
    
    if(sample_b > maxsample){
      maxsample = sample_b
      bestb = x
    }
  }

    # Pull the arm
  outcome <- if_else(rbernoulli(1, probs[bestb]), 1, 0)
  
  # Update priors
  alist[bestb] <- alist[bestb] + outcome
  blist[bestb] <- blist[bestb] + 1 - outcome
  updaterows <- data.frame(matrix(c(rep(i, length(probs)), 1:length(probs), 
                                    alist, blist), nrow=length(probs)))
  colnames(updaterows) <- c("iteration", "Arm", "aval", "bval")
  result_table2 <- result_table2 %>% rbind(updaterows)
}

# Generate pdfs for plotting
iterationplotlist <- c(5:20, seq(20, 100, 5), seq(200, 1000, 50), seq(1200, niter, 200))
# iterationplotlist <- 1:niter
result_table_plot2 <- expand.grid(iterationplotlist, 1:length(probs), x = seq(0, 1, 0.001))
colnames(result_table_plot2) <-  c("iteration", "Arm", "xvalue")
result_table_plot2 <- result_table_plot2 %>%
  inner_join(result_table2) %>%
  mutate(density = dbeta(xvalue, aval, bval),
         Arm = factor(Arm)) %>%
  arrange(iteration, xvalue, Arm)

plot_mab2 <- result_table_plot2 %>%
  ggplot(aes(x = xvalue, y = density, group = Arm, colour = Arm)) +
  geom_line(aes(frame = iteration)) + theme_classic(base_size = 16) +
  # scale_colour_manual(values = cbastylr_main()) +
  scale_x_continuous(breaks = seq(0,1,0.1), limits = c(0, 1)) +
  geom_vline(xintercept = probs, 
             # color = cbastylr_main()[1:4], 
             linetype = 'dashed') +
  # geom_text(x = 0.8, y = 50, label = paste("N = ",iteration)) +
  labs(title = "Multi Armed Bandit - Distribution of response rates",
       caption = "Response rates - Arm 1 = 20%, Arm 2 = 25%, Arm 3 = 26%, Arm 4 = 30%",
       x = "Response Rate", y = "Density")

plotly_mab2 <- ggplotly(plot_mab2)
htmlwidgets::saveWidget(plotly_mab2, file = "plotly_mab2.html")
