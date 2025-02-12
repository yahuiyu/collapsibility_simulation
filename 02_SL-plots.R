###############################################################################
## Project: Collapsibility 
## Program 2:  create plots for SL coefficients distributions

## Created by: Ya-Hui Yu (Feb 13, 2025)
################################################################################

#...............................................................................
# Plots for SL coefficients distribution (DGM with interaction)

coef2 <- read_csv(here("output","coef2_inter.csv"))
coef3<- read_csv(here("output","coef3_inter.csv"))

# Scenario 9 (2 SL algorithms)
coef2.9<- coef2 %>%
  filter(scenario == 9) %>%
  select(SL.lm_All,SL.glm_All) 

n<-nrow(coef2.9)
coef2.9$id<-1:n

coef2.9 <- coef2.9 %>%
  arrange(desc(SL.lm_All)) %>%
  mutate(id = factor(id, levels = rev(id))) %>%
  mutate(cum_pct = (row_number() / n()) * 100)  

coef2.9_long <- coef2.9 %>%
  pivot_longer(cols = c(SL.lm_All,SL.glm_All), 
               names_to = "Algorithm", values_to = "Coefficient")

# Plot stacked horizontal bar chart
p_coef2.9 <- ggplot(coef2.9_long, aes(x = Coefficient, y = id, fill = Algorithm)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("SL.lm_All" = "lightblue1", 
                               "SL.glm_All" = "grey65")) +  
  labs(x = "Proportion", y = "Observations")+
  theme(text = element_text(size = 12),
        axis.text.y = element_blank()) +  
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "lightpink2", size = 0.5) +
  geom_hline(yintercept = 0.5*n, linetype = "dashed", 
             color = "darkblue", size = 0.3)+
  guides(fill = guide_legend(ncol = 2))

#...............................................................................
# Scenario 10 (2 SL algorithms)
coef2.10<- coef2 %>%
  filter(scenario == 10) %>%
  select(SL.lm_All,SL.glm_All) 

n<-nrow(coef2.10)
coef2.10$id<-1:n

coef2.10 <- coef2.10 %>%
  arrange(desc(SL.lm_All)) %>%
  mutate(id = factor(id, levels = rev(id))) %>%
  mutate(cum_pct = (row_number() / n()) * 100)  

coef2.10_long <- coef2.10 %>%
  pivot_longer(cols = c(SL.lm_All,SL.glm_All), 
               names_to = "Algorithm", values_to = "Coefficient")

# Plot stacked horizontal bar chart
p_coef2.10 <- ggplot(coef2.10_long, 
                     aes(x = Coefficient, y = id, fill = Algorithm)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("SL.lm_All" = "lightblue1", 
                               "SL.glm_All" = "grey65")) +  
  labs(x = "Proportion", y = "Observations")+
  theme(text = element_text(size = 12),
        axis.text.y = element_blank()) +  
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "lightpink2", size = 0.5) +
  geom_hline(yintercept = 0.5*n, linetype = "dashed", 
             color = "darkblue", size = 0.3)+
  guides(fill = guide_legend(ncol = 2))

#...............................................................................
# Scenario 11 (3 SL algorithms)
coef3.11<- coef3 %>%
  filter(scenario == 11) %>%
  select(SL.lm_All,SL.glm_All,SL.glm.interaction_All) 

n<-nrow(coef3.11)
coef3.11$id<-1:n

coef3.11 <- coef3.11 %>%
  arrange(desc(SL.lm_All), desc(SL.glm_All), desc(SL.glm.interaction_All)) %>%
  mutate(id = factor(id, levels = rev(id))) %>%
  mutate(cum_pct = (row_number() / n()) * 100)  

coef3.11_long <- coef3.11 %>%
  pivot_longer(cols = c(SL.lm_All,SL.glm_All,SL.glm.interaction_All), 
               names_to = "Algorithm", values_to = "Coefficient")

coef3.11_long$Algorithm <- factor(coef3.11_long$Algorithm, 
                                  levels = c("SL.glm.interaction_All",
                                             "SL.glm_All", "SL.lm_All"))

# Plot stacked horizontal bar chart
p_coef3.11 <- ggplot(coef3.11_long, aes(x = Coefficient, y = id, fill = Algorithm)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("SL.lm_All" = "lightblue1", 
                               "SL.glm_All" = "grey65",
                               "SL.glm.interaction_All" = "cornsilk2")) +  
  labs(x = "Proportion", y = "Observations")+
  theme(text = element_text(size = 12),
        axis.text.y = element_blank()) +  
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "lightpink2", size = 0.5) +
  geom_hline(yintercept = 0.5*n, linetype = "dashed", 
             color = "darkblue", size = 0.3)+
  guides(fill = guide_legend(ncol = 2))

#...............................................................................
# Scenario 12 (3 SL algorithms)
coef3.12<- coef3 %>%
  filter(scenario == 12) %>%
  select(SL.lm_All,SL.glm_All,SL.glm.interaction_All) 

n<-nrow(coef3.12)
coef3.12$id<-1:n

coef3.12 <- coef3.12 %>%
  arrange(desc(SL.lm_All), desc(SL.glm_All), desc(SL.glm.interaction_All)) %>%
  mutate(id = factor(id, levels = rev(id))) %>%
  mutate(cum_pct = (row_number() / n()) * 100)  

coef3.12_long <- coef3.12 %>%
  pivot_longer(cols = c(SL.lm_All,SL.glm_All,SL.glm.interaction_All), 
               names_to = "Algorithm", values_to = "Coefficient")

coef3.12_long$Algorithm <- factor(coef3.12_long$Algorithm, 
                                  levels = c("SL.glm.interaction_All",
                                             "SL.glm_All", "SL.lm_All"))

# Plot stacked horizontal bar chart
p_coef3.12 <- ggplot(coef3.12_long, aes(x = Coefficient, y = id, fill = Algorithm)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("SL.lm_All" = "lightblue1", 
                               "SL.glm_All" = "grey65",
                               "SL.glm.interaction_All" = "cornsilk2")) +  
  labs(x = "Proportion", y = "Observations")+
  theme(text = element_text(size = 12),
        axis.text.y = element_blank()) +  
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "lightpink2", size = 0.5) +
  geom_hline(yintercept = 0.5*n, linetype = "dashed", 
             color = "darkblue", size = 0.3)+
  guides(fill = guide_legend(ncol = 2))

#...............................................................................
## overall figure:
p_combined_inter <-  (p_coef2.9 + labs(title = "a) Linear")) /
  (p_coef3.11 + labs(title = "c) Linear")) |
  (p_coef2.10 + labs(title = "b) Logit"))/
  (p_coef3.12 + labs(title = "d) Logit"))

p_combined_inter <- p_combined_inter + plot_annotation(title = "DGM with interaction")

print(p_combined_inter)

ggsave(here("output","sl_coef_inter.png"))
ggsave(here("output","sl_coef_inter.pdf"))

#-------------------------------------------------------------------------------

# Plots for SL coefficients distribution (DGM without interaction)
coef2 <- read_csv(here("output","coef2_nointer.csv"))
coef3<- read_csv(here("output","coef3_nointer.csv"))

# Scenario 9 (2 SL algorithms)
coef2.9<- coef2 %>%
  filter(scenario == 9) %>%
  select(SL.lm_All,SL.glm_All) 

n<-nrow(coef2.9)
coef2.9$id<-1:n

coef2.9 <- coef2.9 %>%
  arrange(desc(SL.lm_All)) %>%
  mutate(id = factor(id, levels = rev(id))) %>%
  mutate(cum_pct = (row_number() / n()) * 100)  

coef2.9_long <- coef2.9 %>%
  pivot_longer(cols = c(SL.lm_All,SL.glm_All), 
               names_to = "Algorithm", values_to = "Coefficient")

# Plot stacked horizontal bar chart
p_coef2.9 <- ggplot(coef2.9_long, aes(x = Coefficient, y = id, fill = Algorithm)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("SL.lm_All" = "lightblue1", 
                               "SL.glm_All" = "grey65")) +  
  labs(x = "Proportion", y = "Observations")+
  theme(text = element_text(size = 12),
        axis.text.y = element_blank()) +  
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "lightpink2", size = 0.5) +
  geom_hline(yintercept = 0.5*n, linetype = "dashed", 
             color = "darkblue", size = 0.3)+
  guides(fill = guide_legend(ncol = 2))

#...............................................................................
# Scenario 10 (2 SL algorithms)
coef2.10<- coef2 %>%
  filter(scenario == 10) %>%
  select(SL.lm_All,SL.glm_All) 

n<-nrow(coef2.10)
coef2.10$id<-1:n

coef2.10 <- coef2.10 %>%
  arrange(desc(SL.lm_All)) %>%
  mutate(id = factor(id, levels = rev(id))) %>%
  mutate(cum_pct = (row_number() / n()) * 100)  

coef2.10_long <- coef2.10 %>%
  pivot_longer(cols = c(SL.lm_All,SL.glm_All), 
               names_to = "Algorithm", values_to = "Coefficient")

# Plot stacked horizontal bar chart
p_coef2.10 <- ggplot(coef2.10_long, 
                     aes(x = Coefficient, y = id, fill = Algorithm)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("SL.lm_All" = "lightblue1", 
                               "SL.glm_All" = "grey65")) +  
  labs(x = "Proportion", y = "Observations")+
  theme(text = element_text(size = 12),
        axis.text.y = element_blank()) +  
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "lightpink2", size = 0.5) +
  geom_hline(yintercept = 0.5*n, linetype = "dashed", 
             color = "darkblue", size = 0.3)+
  guides(fill = guide_legend(ncol = 2))

#...............................................................................
# Scenario 11 (3 SL algorithms)
coef3.11<- coef3 %>%
  filter(scenario == 11) %>%
  select(SL.lm_All,SL.glm_All,SL.glm.interaction_All) 

n<-nrow(coef3.11)
coef3.11$id<-1:n

coef3.11 <- coef3.11 %>%
  arrange(desc(SL.lm_All), desc(SL.glm_All), desc(SL.glm.interaction_All)) %>%
  mutate(id = factor(id, levels = rev(id))) %>%
  mutate(cum_pct = (row_number() / n()) * 100)  

coef3.11_long <- coef3.11 %>%
  pivot_longer(cols = c(SL.lm_All,SL.glm_All,SL.glm.interaction_All), 
               names_to = "Algorithm", values_to = "Coefficient")

coef3.11_long$Algorithm <- factor(coef3.11_long$Algorithm, 
                                  levels = c("SL.glm.interaction_All",
                                             "SL.glm_All", "SL.lm_All"))

# Plot stacked horizontal bar chart
p_coef3.11 <- ggplot(coef3.11_long, aes(x = Coefficient, y = id, fill = Algorithm)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("SL.lm_All" = "lightblue1", 
                               "SL.glm_All" = "grey65",
                               "SL.glm.interaction_All" = "cornsilk2")) +  
  labs(x = "Proportion", y = "Observations")+
  theme(text = element_text(size = 12),
        axis.text.y = element_blank()) +  
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "lightpink2", size = 0.5) +
  geom_hline(yintercept = 0.5*n, linetype = "dashed", 
             color = "darkblue", size = 0.3)+
  guides(fill = guide_legend(ncol = 2))

#...............................................................................
# Scenario 12 (3 SL algorithms)
coef3.12<- coef3 %>%
  filter(scenario == 12) %>%
  select(SL.lm_All,SL.glm_All,SL.glm.interaction_All) 

n<-nrow(coef3.12)
coef3.12$id<-1:n

coef3.12 <- coef3.12 %>%
  arrange(desc(SL.lm_All), desc(SL.glm_All), desc(SL.glm.interaction_All)) %>%
  mutate(id = factor(id, levels = rev(id))) %>%
  mutate(cum_pct = (row_number() / n()) * 100)  

coef3.12_long <- coef3.12 %>%
  pivot_longer(cols = c(SL.lm_All,SL.glm_All,SL.glm.interaction_All), 
               names_to = "Algorithm", values_to = "Coefficient")

coef3.12_long$Algorithm <- factor(coef3.12_long$Algorithm, 
                                  levels = c("SL.glm.interaction_All",
                                             "SL.glm_All", "SL.lm_All"))

# Plot stacked horizontal bar chart
p_coef3.12 <- ggplot(coef3.12_long, aes(x = Coefficient, y = id, fill = Algorithm)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("SL.lm_All" = "lightblue1", 
                               "SL.glm_All" = "grey65",
                               "SL.glm.interaction_All" = "cornsilk2")) +  
  labs(x = "Proportion", y = "Observations")+
  theme(text = element_text(size = 12),
        axis.text.y = element_blank()) +  
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "lightpink2", size = 0.5) +
  geom_hline(yintercept = 0.5*n, linetype = "dashed", 
             color = "darkblue", size = 0.3)+
  guides(fill = guide_legend(ncol = 2))

#...............................................................................
## overall figure:
p_combined_nointer <-  (p_coef2.9 + labs(title = "a) Linear")) /
  (p_coef3.11 + labs(title = "c) Linear")) |
  (p_coef2.10 + labs(title = "b) Logit"))/
  (p_coef3.12 + labs(title = "d) Logit"))


p_combined_nointer <- p_combined_nointer + 
                    plot_annotation(title = "DGM without interaction")

print(p_combined_nointer)

ggsave(here("output","sl_coef_nointer.png"))
ggsave(here("output","sl_coef_nointer.pdf"))
