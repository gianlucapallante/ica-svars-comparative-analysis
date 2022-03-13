## Code to replicate Figure 3 of the paper
## Load packege here for managing workflow

if(length(which(installed.packages() %in% c("rstudioapi") == T)) == 0){
  install.packages("rstudioapi")
  library(rstudioapi)
}else{library(rstudioapi)}

if(length(which(installed.packages() %in% c("here") == T)) == 0){
  install.packages("here")
  library(here)
}else{library(here)}

## Open ica-svars-comparative-analysis.Rproj
if(!base::grepl(x = here(),pattern = "ica-svars-comparative-analysis")){
  file_path <- file.choose()
}

if(!base::grepl(x = here(),pattern = "ica-svars-comparative-analysis")){
  openProject(file_path)
}

## Check that working directory contains the ica-svars-comparative-analysis.Rproj file
library(here)
dr_here()
here()
## Independent Component Analysis with different methods (MC analysis)
source(paste0(here(),"/Rpackages.R"),local = T)


## Colors
palette <- c('CvM'     = '#66c2a5',
             'fastICA' = '#fc8d62',
             'DCov'    = '#8da0cb',
             'PML'     = '#e78ac3')


## Parameters of the p-generalized distribution
## Number of shape parameter values (in the specific assessment and for inference table we only focuse on 4 distributional scenario - See paper)
seq1   <- seq(from = 0.5, to = 3.5, length.out = 15) ## Sequence of p-Generalized parameter for general evaluation exercise (1st part)
seq2   <- seq(from = 4, to = 100, length.out = 5) ## Sequence of p-Generalized parameter for general evaluation exercise (2nd part)
#SEQ     <- c(0.5,1.57,2.43,100) ## Sequence of p-Generalized parameter for statistical inference exercise
SEQ    <- as.numeric(c(seq1, seq2))

rm(seq1,seq2)


## Load Databases with simulated data
average_performance_spec <- read_csv(file = here("final_databases","average_performance_specific.csv")) %>% 
  dplyr::select(-"X1")

errors <- read_csv(file =  here("final_databases","p_value_errors.csv")) %>% 
  dplyr::select(-"X1")


errors %>% group_by(dimension,column,scenario,p_shape) %>% 
  left_join(average_performance_spec,.) %>% 
  filter(scenario != "n=200") %>% 
  rename(MDI = "Minimum Distance") %>% 
  group_by(estimator,p_shape,scenario,dimension) %>% 
  summarise_at(vars(MDI,p.value), list(mean = ~ mean(.),
                                              sd   = ~ sd(.))) %>% 
  ungroup() %>% 
  mutate(p_shape = as.factor(p_shape)) %>% 
  ggplot(.,aes(x = p_shape, group = estimator))+
  geom_line(aes(y = MDI_mean,color = estimator), size = 1)+
  geom_point(aes(y = MDI_mean,color = estimator, shape = estimator), size = 3)+
  geom_ribbon(data = . %>% filter(p.value_mean>0.05), 
              aes(xmin = p_shape[min(p.value_mean) & estimator == estimator], xmax = p_shape[max(p.value_mean)& estimator == estimator],
                  ymin = -Inf, ymax = Inf),
              alpha=0.3, fill = "grey" )+
  facet_grid(dimension~scenario)+
  scale_color_manual(values = palette, breaks = c( "CvM","DCov","fastICA","PML"))+
  scale_fill_manual(values = palette, breaks = c( "CvM","DCov","fastICA","PML"))+
  scale_x_discrete(breaks = c("0.5","1.57","2","2.43","4","52","100"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "white"))+
  labs(title = "Specific Assessment", y = "MDI")
