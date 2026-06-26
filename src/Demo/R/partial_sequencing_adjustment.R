library(mgcv)
library(ggsci)

# read in sequencing proportion
# ToDO: Should get 2020 pcrpositives_sequenced time series.

get_partial_sequencing_adjustment <- function()
{
  proportion_of_pcrpositives_sequenced <- read_csv("K:/SSI_sekvensering_corona/christian/Data_from_jacob/Sequencing_proportions/proportion_of_pcrpositives_sequenced.csv", 
                                                   col_names = FALSE)
  proportion_of_pcrpositives_sequenced$date <- as.Date('2020-12-31') + days(x = seq.int(1,365))
  proportion_of_pcrpositives_sequenced <- proportion_of_pcrpositives_sequenced %>% rename(pct=X1) %>% dplyr::select(date,pct)
  
  # smooth with GAM and/or emwa (try both)
  x        <- matrix(c(1:365,(1:365 %% 7)),ncol=2)
  dat      <- as_tibble(data.frame(y=proportion_of_pcrpositives_sequenced$pct,x1=x[,1],x2=as.factor(x[,2])))
  gam_prop <- gam(y~s(x1,bs='cr',k=50)+s(x2,bs='fs',k=7),data=dat)
  
  gam_pred <- predict(gam_prop,data=dat)
  proportion_of_pcrpositives_sequenced$gam_pred <- gam_pred
  
  lambda <- 0.2
  proportion_of_pcrpositives_sequenced %>% 
    mutate(pct_ewma = accumulate(pct, ~ lambda * .y + (1-lambda)*.x)) %>%
    pivot_longer(-date,values_to = 'pct',names_to = 'smooth') %>%
    ggplot(aes(x=date, y=pct,col=smooth)) + geom_line(aes(linewidth=smooth)) +
    scale_discrete_manual('linewidth',values = c(2,0.5,1))
  
  # read in iar's
  serial.interval <- read_csv('K:/SSI_sekvensering_corona/christian/R_projects/danish_npi_project/data/serial_interval.csv')
  si_short        <- c(serial.interval$fit[1:29],serial.interval$fit[30:100] |> sum())
  rt_iar_file     <- read_csv('K:/SSI_sekvensering_corona/christian/R_projects/danish_npi_project/data/denmark_rt_iar.csv') %>% 
    mutate(Rc    = rollapply(Rt_adj,width=30,w=si_short,FUN=weighted.mean,fill=NA,align = 'left'),
           Rc_lb = rollapply(Rt_adj_0.025,width=30,w=si_short,FUN=weighted.mean,fill=NA,align = 'left'),
           Rc_ub = rollapply(Rt_adj_0.975,width=30,w=si_short,FUN=weighted.mean,fill=NA,align = 'left'),
           iar_shifted = rollapply(iar_0.5,width=30,w=si_short,FUN=weighted.mean,fill=NA,align = 'left'),
           iar_lb_shifted = rollapply(iar_0.025,width=30,w=si_short,FUN=weighted.mean,fill=NA,align = 'left'),
           iar_ub_shifted = rollapply(iar_0.975,width=30,w=si_short,FUN=weighted.mean,fill=NA,align = 'left'))
  
  rt_iar_file <- rt_iar_file %>% left_join(proportion_of_pcrpositives_sequenced,by=c('date')) %>%
    mutate(gam_pred_adj = replace_na(gam_pred,1),
           partial_sequencing = gam_pred_adj*iar_0.5,
           partial_sequencing_shifted = gam_pred_adj*iar_shifted)
  
  return(rt_iar_file)
}



# combine both time series
#proportion_of_pcrpositives_sequenced <- proportion_of_pcrpositives_sequenced %>% left_join(rt_iar_file %>% dplyr::select(date,contains('iar'))) %>%
#  mutate(partial_sequencing=gam_pred*iar_0.5 )

if(FALSE)
{
  pp_seq <- get_partial_sequencing_adjustment() %>% filter(date>as.Date('2020-09-01') & date < as.Date('2022-01-01'))
  pp_plt <- pp_seq  %>% dplyr::select(date,iar,pct,gam_pred) %>%
    rename(`Sequenced proportion of\n infection (smoothed)`=gam_pred,`Infection assertainment\n ratio`=iar, `% of cases\n sequenced`=pct )%>%
    pivot_longer(-date,values_to = 'pct',names_to = 'Measure') %>%
    ggplot(aes(x=date, y=pct,col=Measure)) + geom_line(aes(linewidth=Measure)) + scale_color_lancet() + 
    scale_discrete_manual('linewidth',values = c(0.5,1,1.5)) + 
    geom_ribbon(data=pp_seq,aes(x=date,ymin=iar_0.025,ymax=iar_0.975,col='Infection assertainment\n ratio',fill='Infection assertainment\n ratio'), alpha=0.2, linetype = "dotted") +
    theme_light() + xlab('') + ylab('Percentage (%)') + guides(fill='none')
  ggsave('K://SSI_sekvensering_corona/christian/python_scripts/seq_prop_plot.pdf',pp_plt,width=11,height=6)
  
  geom_ribbon(aes(x=SampleDate,ymin=iar_0.025, ymax=iar_0.975), alpha=0.2, linetype = "dotted")
  
  proportion_of_pcrpositives_sequenced %>% dplyr::select(!contains('iar')) %>%
    mutate(pct_ewma = accumulate(pct, ~ lambda * .y + (1-lambda)*.x)) %>%
    pivot_longer(-date,values_to = 'pct',names_to = 'smooth') %>%
    ggplot(aes(x=date, y=pct,col=smooth)) + geom_line(aes(linewidth=smooth)) + scale_color_lancet() + 
    scale_discrete_manual('linewidth',values = c(2,1.5,0.5,1)) + theme_light() 
  
  proportion_of_pcrpositives_sequenced %>% dplyr::select(date,pct,gam_pred,iar,partial_sequencing) %>%
    mutate(log_ps     = log(partial_sequencing),
           log_one_ps = log(1-partial_sequencing))
  
  # Check how many clusters we have per week + percentage of clusters with outdegree 0 (that is no onwards transmission)
  # reg_data %>% group_by(isoweek) %>% summarise(n=n(),singleton=sum(out_degree==0)/n) %>% View()
  
  # run simulations to see how much adjustment would be
  get_log_proba_size_cluster <- function(cluster_size,R,k,p=1)
  {
    return(lgamma(k*cluster_size+cluster_size-1) -
             lgamma(k*cluster_size) -
             lgamma(cluster_size+1) +
             (cluster_size-1)*log(p*R/k) -
             (k*cluster_size+cluster_size -1) * log(1+p*R/k))
  }
  
  get_proba_size_cluster <- function(cluster_size,R,k,p=1)
  {
    return(gamma(k*cluster_size+cluster_size-1) / 
             (gamma(k*cluster_size)*
                gamma(cluster_size+1) ) *
             (p*R/k)^(cluster_size-1) *
             (1+p*R/k)^(1-k*cluster_size-cluster_size) )
  }
  
  l <- 20
  j <- 1:20
  p_detect <- as_tibble(data.frame(p_detect=c(0.01,0.1,0.5))) 
  true_R   <- as_tibble(data.frame(true_R=seq(0.5,1.75,by=0.25))) 
  true_k   <- as_tibble(data.frame(true_k=c(0.01,0.1,0.5,1,2,5))) 
  
  vec_proba_clust_size <- as_tibble(data.frame(j=j)) %>% cross_join(true_R) %>% cross_join(true_k) %>%
    mutate(vec_log_proba_clust_size = get_log_proba_size_cluster(j,true_R,true_k),
           vec_proba_clust_size     = exp(vec_log_proba_clust_size),
           blah                     = get_proba_size_cluster(j,true_R,true_k))
  
  log_proba_helper <- as_tibble(data.frame(j=j)) %>%
    mutate(l          = l,
           l_choose_j = lchoose(l,j),
           l_minus_j  = l-j ) %>%
    dplyr::select(l,j,l_choose_j,l_minus_j) %>% cross_join(p_detect) %>%
    mutate(log_p_detect = log(p_detect),
           log_1_minus_p_detect = log(1-p_detect)) %>%
    mutate(sum_of_terms = l_choose_j + j*log_p_detect + l_minus_j*log_1_minus_p_detect,
           prob = exp(sum_of_terms),
           prob_ub = exp(log(1)+sum_of_terms)) # this assumes prob of cluster of this size is 1
  
  
  # prob of no onward transmission:
  prob_obs_zero_elements <- vec_proba_clust_size %>% filter(j==1) %>% group_by(true_R,true_k) %>% summarize(prob_obs_zero_elements=0)
  for(i in j)
  {
    tmp <- vec_proba_clust_size %>% filter(j==i) %>% group_by(true_R,true_k) %>%
      summarize(tmp = blah*(1-0.1)^j) 
    
    prob_obs_zero_elements$prob_obs_zero_elements <- prob_obs_zero_elements$prob_obs_zero_elements + tmp$tmp
  }
  
}


