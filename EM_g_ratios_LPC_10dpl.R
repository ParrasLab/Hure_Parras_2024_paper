setwd("...")
library(ggplot2)
library(ggpubr)
library(dplyr)
library(multcomp)
data <- read.table("EM_gratios.txt",sep="\t", header=TRUE,  fill=TRUE, quote = "")
head(data)
tail(data)
dim(data) #753   3
dat <- subset(data, axon_diameter < 2.7) #exclude 5 outliers to the largest axonal diameter
head(dat)
dim(dat) # 748   3



dat$condition <- factor(dat$condition, levels = c("Vehicle","Clem", "Sm5", "Sm11"))
   colors<- c("gray80","gray40", "lightskyblue", "dodgerblue3")
   conditions <- c("Vehicle","Clem", "Sm5", "Sm11")
   
   #To convert into factors the character (<chr>) columns, specifying the order given by treatment and stage vectors
   dat$condition <- factor(dat$condition, levels =conditions)

  
   library(ggpubr)
    ggplot(dat, aes(x = axon_diameter, y = g_ratio, color = condition)) +
     geom_point(size = 2, alpha = 0.9) +
     scale_color_manual(values = setNames(colors, conditions)) +
     labs(x = "Axon diameter", y = "G-ratio", color = "Condition") +
     scale_x_continuous(breaks = seq(0.5, 3, by = 0.5)) +
     scale_y_continuous(breaks = seq(0.65, 0.95, by = 0.05)) +
     geom_smooth(method = "lm", se = FALSE) +
     theme(panel.grid.major = element_line(color = "gray70"),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           axis.line = element_line(color = "gray70"),
           legend.position = "right") +
     theme_classic() +
     guides(color = guide_legend(reverse = F)) #+
     #stat_cor(label.x = 1.9, label.y= c(0.66, 0.64, 0.62, 0.60))
   ggsave("scatter_plot_lesion_dat.pdf", width = 10, height = 7)
   
 
   
   #############################################################################
   ###### Axons from 1 to 1.5 > p.val = 0.0164 * with Sm11/Leurcovorin
   ############################################################################# 
   
   filtered_data <- dat %>%
     filter(axon_diameter >= 1 & axon_diameter <= 1.5)
   
   
   violin_plot <- ggplot(filtered_data, aes(x = condition, y = g_ratio, color = condition)) +
     geom_violin(trim = FALSE, alpha = 0.5, size = 1.1) +  # Increase size here
     geom_jitter(width = 0.3, size =3, alpha = 0.5) +
     labs(x = "Condition", y = "G-ratio", color = "Condition") +
     ggtitle("Violin Plot of g-ratios from axon diameters from 1 to 1.5") +
     scale_color_manual(values = setNames(colors, conditions)) +
     theme_classic() + # Change theme to classic
     
     theme(
       plot.title = element_text(hjust = 0.5, margin = margin(0, 0, 10, 0))
     )
   print(violin_plot)
   ggsave("violin_plot_1to15_dim.pdf", width = 8, height = 7)
   
  # Perform ANOVA. the code uses a one-way ANOVA to compare the cell counts across different treatments. 
   outaov <- aov(g_ratio ~ condition, data = filtered_data)
    summary(outaov) # F(3,257)=4.60, p=0.004
    #                 Df Sum Sq  Mean Sq F value  Pr(>F)   
    # condition     3 0.0418 0.013941   4.599 0.00373 **
    #   Residuals   257 0.7790 0.003031 
    outaov <- aov(g_ratio ~ condition, data = filtered_data)
    
    # Perform Dunnett's test to compare each treatment with Vehicle
    # Specify "Vehicle" as the control group
    dunnett_test <- multcomp::glht(outaov, linfct = mcp(condition = "Dunnett"))
   
    # Get summary of the Dunnett test
    summary(dunnett_test)
    # 
    # Simultaneous Tests for General Linear Hypotheses
    # 
    # Multiple Comparisons of Means: Dunnett Contrasts
    # 
    # 
    # Fit: aov(formula = g_ratio ~ condition, data = filtered_data)
    # 
    # Linear Hypotheses:
    #                      Estimate Std. Error t value Pr(>|t|)  
    #                      Estimate Std. Error t value Pr(>|t|)  
    # Clem - Vehicle == 0 -0.003355   0.009668  -0.347   0.9718  
    # Sm5 - Vehicle == 0   0.006365   0.009929   0.641   0.8557  
    # Sm11 - Vehicle == 0 -0.026655   0.009602  -2.776   0.0164 *
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # (Adjusted p values reported -- single-step method)
    

   
  ###########################################################################
  ###### Axons from larger than 1 p= 0.00555 for with Sm11/Leurcovorin
  ########################################################################### 
  
  filtered_data <- dat %>%
    filter(axon_diameter >= 1 )
  
  
  violin_plot <- ggplot(filtered_data, aes(x = condition, y = g_ratio, color = condition)) +
    geom_violin(trim = FALSE, alpha = 0.5, size = 1.1) +  # Increase size here
    geom_jitter(width = 0.3, size =3, alpha = 0.5) +
    labs(x = "Condition", y = "G-ratio", color = "Condition") +
    ggtitle("Violin Plot of g-ratios from axon diameters from 1 to 1.5") +
    scale_color_manual(values = setNames(colors, conditions)) +
    theme_classic() + # Change theme to classic
    
    theme(
      plot.title = element_text(hjust = 0.5, margin = margin(0, 0, 10, 0))
    )
  print(violin_plot)
  ggsave("violin_plot_more_than_1_dim.pdf", width = 8, height = 7)
  
  # Perform ANOVA. the code uses a one-way ANOVA to compare the cell counts across different treatments. 
  outaov <- aov(g_ratio ~ condition, data = filtered_data)
  summary(outaov) # F(3,257)=4.60, p=0.004
  #                Df Sum Sq  Mean Sq F value  Pr(>F)   
  # condition       3 0.0461 0.015365   4.611 0.00349 **
  #   Residuals   389 1.2963 0.003332  
  outaov <- aov(g_ratio ~ condition, data = filtered_data)
  
  # Perform Dunnett's test to compare each treatment with Vehicle
  # Specify "Vehicle" as the control group
  dunnett_test <- multcomp::glht(outaov, linfct = mcp(condition = "Dunnett"))
  
  # Get summary of the Dunnett test
  summary(dunnett_test)
  
  # Fit: aov(formula = g_ratio ~ condition, data = filtered_data)
  # 
  # Linear Hypotheses:
  #   Estimate Std. Error t value Pr(>|t|)   
  # Clem - Vehicle == 0 -0.005817   0.008271  -0.703  0.82112   
  # Sm5 - Vehicle == 0   0.000865   0.008566   0.101  0.99925   
  # Sm11 - Vehicle == 0 -0.025181   0.008068  -3.121  0.00555 **
  #   ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # (Adjusted p values reported -- single-step method)

  
 