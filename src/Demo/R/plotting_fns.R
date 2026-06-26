library(rmutil)

plot_fits_vs_data_household <- function(data,isoweek_in='all',random_id_in='ML',label='A')
{
  setting = 'address'
  params = read_csv('data/Rc_param_estimates_ML_random_trees_prioritise_settings_by_setting_full_period_v4.csv')
  
  if(isoweek_in != 'all')
    data      <- data %>% filter(isoweek==isoweek_in)
  
  
  Rc          <- data %>% pull(paste('out_degree',setting,sep="_"))
  data_list   <- list(N=length(Rc),y=Rc,hsize=data$hsize_large)
  params      <- params %>% filter(isoweek==isoweek_in & str_ends(random_id, paste(random_id_in, setting, sep='_') ) )
  
  hist_nb <- as_tibble(list(out_degree=rnbinom(data_list$N,mu=params %>% filter(model=='negative_binomial' & parameter == 'R') %>% pull(mean),
                                               size=params %>% filter(model=='negative_binomial' & parameter == 'k') %>% pull(mean)))) %>% 
    group_by(out_degree) %>%
    summarise(count=n())
  hist_poisson <- as_tibble(list(out_degree=rpois(data_list$N,params %>% filter(model=='poisson' & parameter == 'lambda') %>% pull(mean)))) %>% 
     group_by(out_degree) %>%
     summarise(count=n())
  alpha = params %>% filter(model=='beta_binomial_partial_obs' & parameter == 'alpha') %>% pull(mean)
  beta  = params %>% filter(model=='beta_binomial_partial_obs' & parameter == 'beta') %>% pull(mean)
  hist_bb <- as_tibble(list(out_degree=c(rbetabinom( n=data_list$N, size=2.5,m=alpha/(alpha+beta),s=alpha+beta)))) %>%
    group_by(out_degree) %>%
    summarise(count=n())
  hist_zinb <- as_tibble(list(out_degree=c(rnbinom(data_list$N*(1-params %>% filter(model=='zero_inflated_negative_binomial' & parameter == 'theta') %>% pull(mean)),
                                                   mu=params %>% filter(model=='zero_inflated_negative_binomial' & parameter == 'R') %>% pull(mean),
                                                   size=params %>% filter(model=='zero_inflated_negative_binomial' & parameter == 'k') %>% pull(mean)),
                                           rep(0,data_list$N*params %>% filter(model=='zero_inflated_negative_binomial' & parameter == 'theta') %>% pull(mean))))) %>%
    group_by(out_degree) %>%
    summarise(count=n())
  
  plt <- data %>% 
    mutate(out_degree = !!sym(paste('out_degree',setting, sep='_')) ) %>%
    ggplot(aes(x=out_degree)) + geom_histogram(fill='lightgrey') +
    geom_line(data=hist_nb,aes(y=count,colour='NegBinomial')) + 
    geom_point(data=hist_nb,aes(y=count),col='red',shape=18) + 
    geom_line(data=hist_poisson,aes(y=count,colour='Poisson')) +
    geom_point(data=hist_poisson,aes(y=count),col='blue',shape=4) +
    geom_line(data=hist_bb,aes(y=count,colour='Beta Binomial')) +
    geom_point(data=hist_bb,aes(y=count),col='darkgreen',shape=14) +
    geom_line(data=hist_zinb,aes(y=count,colour='ZI NegBinomial')) +
    geom_point(data=hist_zinb,aes(y=count),col='black',shape=16) +
    scale_color_manual('',breaks=c('NegBinomial','Poisson','Beta Binomial','ZI NegBinomial'),values=c('red','blue','darkgreen','black')) +
    theme_light() + theme(legend.position = 'bottom') + ggtitle(paste0('(', label, ') ', 'Estimates for week ',isoweek_in, " and id=", random_id_in ))
  
  
  return(plt)
}

plot_fits_vs_data_setting <- function(data,isoweek_in='all',random_id_in='ML',setting='school', label='A')
{
  params = read_csv('data/Rc_param_estimates_ML_random_trees_prioritise_settings_by_setting_full_period_v4.csv')
  
  if(isoweek_in != 'all')
    data      <- data %>% filter(isoweek==isoweek_in)
  
  
  Rc          <- data %>% pull(paste('out_degree',setting,sep="_"))
  data_list   <- list(N=length(Rc),y=Rc,hsize=data$hsize_large)
  params      <- params %>% filter(isoweek==isoweek_in & str_ends(random_id, paste(random_id_in, setting, sep='_') ) )
  
  hist_nb <- as_tibble(list(out_degree=rnbinom(data_list$N,mu=params %>% filter(model=='negative_binomial' & parameter == 'R') %>% pull(mean),
                                               size=params %>% filter(model=='negative_binomial' & parameter == 'k') %>% pull(mean)))) %>% 
    group_by(out_degree) %>%
    summarise(count=n())
  hist_poisson <- as_tibble(list(out_degree=rpois(data_list$N,params %>% filter(model=='poisson' & parameter == 'lambda') %>% pull(mean)))) %>% 
    group_by(out_degree) %>%
    summarise(count=n())
  alpha = params %>% filter(model=='beta_binomial_partial_obs' & parameter == 'alpha') %>% pull(mean)
  beta  = params %>% filter(model=='beta_binomial_partial_obs' & parameter == 'beta') %>% pull(mean)
  hist_zinb <- as_tibble(list(out_degree=c(rnbinom(data_list$N*(1-params %>% filter(model=='zero_inflated_negative_binomial' & parameter == 'theta') %>% pull(mean)),
                                                   mu=params %>% filter(model=='zero_inflated_negative_binomial' & parameter == 'R') %>% pull(mean),
                                                   size=params %>% filter(model=='zero_inflated_negative_binomial' & parameter == 'k') %>% pull(mean)),
                                           rep(0,data_list$N*params %>% filter(model=='zero_inflated_negative_binomial' & parameter == 'theta') %>% pull(mean))))) %>%
    group_by(out_degree) %>%
    summarise(count=n())
  
  plt <- data %>% 
    mutate(out_degree = !!sym(paste('out_degree',setting, sep='_')) ) %>%
    ggplot(aes(x=out_degree)) + geom_histogram(fill='lightgrey') +
    geom_line(data=hist_nb,aes(y=count,colour='NegBinomial')) + 
    geom_point(data=hist_nb,aes(y=count),col='red',shape=18) + 
    geom_line(data=hist_poisson,aes(y=count,colour='Poisson')) +
    geom_point(data=hist_poisson,aes(y=count),col='blue',shape=4) +
    geom_line(data=hist_zinb,aes(y=count,colour='ZI NegBinomial')) +
    geom_point(data=hist_zinb,aes(y=count),col='black',shape=16) +
    scale_color_manual('',breaks=c('NegBinomial','Poisson','Beta Binomial','ZI NegBinomial'),values=c('red','blue','darkgreen','black')) +
    theme_light() + theme(legend.position = 'bottom') + ggtitle(subtitle=paste0('(', label, ') ', 'Estimates for week ',isoweek_in, " and id=", random_id_in ))
  
  
  return(plt)
}

plot_fits_vs_data_setting <- function(params=NULL, data=node_attrs, isoweek_in='all',random_id_in='ML', setting = 'school',label='A')
{
  if(is.null(params))
    params = read_csv('data/Rc_param_estimates_ML_random_trees_prioritise_settings_by_setting_full_period_v3.csv')
  
  if(isoweek_in!='all')
    data      <- data %>% filter(isoweek==isoweek_in)
  
  
  Rc          <- data %>% pull(paste('out_degree',setting,sep="_"))
  data_list   <- list(N=length(Rc),y=Rc)
  params      <- params %>% filter(isoweek==isoweek_in & str_ends( random_id, paste(random_id_in, setting, sep='_') ) )
  
  hist_nb <- as_tibble(list(out_degree=rnbinom(data_list$N,mu=params %>% filter(model=='negative_binomial' & parameter == 'R') %>% pull(mean),
                                               size=params %>% filter(model=='negative_binomial' & parameter == 'k') %>% pull(mean)))) %>% 
    group_by(out_degree) %>%
    summarise(count=n())
  # hist_poisson <- as_tibble(list(out_degree=rpois(data_list$N,params %>% filter(model=='poisson' & parameter == 'lambda') %>% pull(mean)))) %>% 
  #   group_by(out_degree) %>%
  #   summarise(count=n())
  # hist_zip <- as_tibble(list(out_degree=c(rpois(data_list$N*(1-params %>% filter(model=='zero_inflated_poisson' & parameter == 'theta') %>% pull(mean)),
  #                                               params %>% filter(model=='zero_inflated_poisson' & parameter == 'lambda') %>% pull(mean)),
  #                                         rep(0,data_list$N*params %>% filter(model=='zero_inflated_poisson' & parameter == 'theta') %>% pull(mean))))) %>% 
  #   group_by(out_degree) %>%
  #   summarise(count=n())
  hist_zinb <- as_tibble(list(out_degree=c(rnbinom(data_list$N*(1-params %>% filter(model=='zero_inflated_negative_binomial' & parameter == 'theta') %>% pull(mean)),
                                                   mu=params %>% filter(model=='zero_inflated_negative_binomial' & parameter == 'R') %>% pull(mean),
                                                   size=params %>% filter(model=='zero_inflated_negative_binomial' & parameter == 'k') %>% pull(mean)),
                                           rep(0,data_list$N*params %>% filter(model=='zero_inflated_negative_binomial' & parameter == 'theta') %>% pull(mean))))) %>%
    group_by(out_degree) %>%
    summarise(count=n())
  
  plt <- data %>% 
    mutate(out_degree = !!sym(paste('out_degree',setting, sep='_')) ) %>%
    ggplot(aes(x=out_degree)) + geom_histogram(fill='lightgrey') +
    geom_line(data=hist_nb,aes(y=count,colour='NegBinomial')) + 
    geom_point(data=hist_nb,aes(y=count),col='red',shape=18) + 
    #geom_line(data=hist_poisson,aes(y=count,colour='Poisson')) +
    #geom_point(data=hist_poisson,aes(y=count),col='blue',shape=4) +
    #geom_line(data=hist_zip,aes(y=count,colour='ZI Poisson')) +
    #geom_point(data=hist_zip,aes(y=count),col='darkgreen',shape=14) +
    geom_line(data=hist_zinb,aes(y=count,colour='ZI NegBinomial')) +
    geom_point(data=hist_zinb,aes(y=count),col='black',shape=16) +
    scale_color_manual('',breaks=c('NegBinomial','Poisson','ZI Poisson','ZI NegBinomial'),values=c('red','blue','darkgreen','black')) +
    theme_light() + theme(legend.position = 'bottom') + ggtitle(paste0('(', label, ') ', 'Estimates for week ',isoweek_in, " and id=", random_id_in ))
  
  
  return(plt)
}


plot_fits_vs_data <- function(params=k_results, data=node_attrs, isoweek_in='all',random_id_in='ML',label='A')
{
  if(isoweek_in!='all')
    data      <- data %>% filter(isoweek==isoweek_in)
  
  
  Rc          <- data %>% pull(out_degree)
  data_list   <- list(N=length(Rc),y=Rc)
  params      <- params %>% filter(isoweek==isoweek_in & random_id == random_id_in )
  #params      <- params %>% filter(isoweek==isoweek_in & str_detect(random_id,random_id_in) & str_detect(random_id,'_all') )
  
  hist_nb <- as_tibble(list(out_degree=rnbinom(data_list$N,mu=params %>% filter(model=='negative_binomial' & parameter == 'R') %>% pull(mean),
                                               size=params %>% filter(model=='negative_binomial' & parameter == 'k') %>% pull(mean)))) %>% 
    group_by(out_degree) %>%
    summarise(count=n())
  hist_poisson <- as_tibble(list(out_degree=rpois(data_list$N,params %>% filter(model=='poisson' & parameter == 'lambda') %>% pull(mean)))) %>% 
    group_by(out_degree) %>%
    summarise(count=n())
  hist_zip <- as_tibble(list(out_degree=c(rpois(data_list$N*(1-params %>% filter(model=='zero_inflated_poisson' & parameter == 'theta') %>% pull(mean)),
                                                params %>% filter(model=='zero_inflated_poisson' & parameter == 'lambda') %>% pull(mean)),
                                          rep(0,data_list$N*params %>% filter(model=='zero_inflated_poisson' & parameter == 'theta') %>% pull(mean))))) %>% 
    group_by(out_degree) %>%
    summarise(count=n())
  hist_zinb <- as_tibble(list(out_degree=c(rnbinom(data_list$N*(1-params %>% filter(model=='zero_inflated_negative_binomial' & parameter == 'theta') %>% pull(mean) |> mean() ),
                                                   mu=params %>% filter(model=='zero_inflated_negative_binomial' & parameter == 'R') %>% pull(mean) |> mean(),
                                                   size=params %>% filter(model=='zero_inflated_negative_binomial' & parameter == 'k') %>% pull(mean)|> mean()),
                                           rep(0,data_list$N*params %>% filter(model=='zero_inflated_negative_binomial' & parameter == 'theta') %>% pull(mean)|> mean() )))) %>%
    group_by(out_degree) %>%
    summarise(count=n())
  
  plt <- data %>% ggplot(aes(x=out_degree)) + geom_histogram(fill='lightgrey') +
    geom_line(data=hist_nb,aes(y=count,colour='NegBinomial')) + 
    geom_point(data=hist_nb,aes(y=count),col='red',shape=18) + 
    geom_line(data=hist_poisson,aes(y=count,colour='Poisson')) +
    geom_point(data=hist_poisson,aes(y=count),col='blue',shape=4) +
    geom_line(data=hist_zip,aes(y=count,colour='ZI Poisson')) +
    geom_point(data=hist_zip,aes(y=count),col='darkgreen',shape=14) +
    geom_line(data=hist_zinb,aes(y=count,colour='ZI NegBinomial')) +
    geom_point(data=hist_zinb,aes(y=count),col='black',shape=16) +
    scale_color_manual('',breaks=c('NegBinomial','Poisson','ZI Poisson','ZI NegBinomial'),values=c('red','blue','darkgreen','black')) +
    theme_light() + theme(legend.position = 'bottom') + ggtitle(paste0('(', label, ') ', 'Estimates for week ',isoweek_in, " and id=", random_id_in ))
  

  return(plt)
}


# Only want to plot Rc as the Rt we can convert to Rc as we have daily data
plot_Rc <- function(model_in='negative_binomial',node_isoweek_dates,Rc_agg=NULL,max_date=as.Date('2021-10-01'),min_date=as.Date('2020-10-01'),sim_data=NULL) 
{
  # plot population level rather an output of model.
  k_results <- read_csv('stan_output/Rc_param_estimates.csv')
  
  # 1 for ML tree
  r_plot_data <- k_results %>% filter(random_id=='ML' & isoweek!='all' & model == model_in) %>% 
    filter(parameter!='lp__') %>% 
    dplyr::select(parameter,mean,`2.5%`,`97.5%`,isoweek) %>%
    pivot_wider(names_from = parameter, values_from = c(mean,`2.5%`,`97.5%`)) %>%
    pivot_longer(
      -isoweek,
      names_to = c("quantile", "parameter"),
      names_pattern = "(.*)_(.*)",
      values_to = "value"
    ) %>% 
    filter(parameter %in% c('R')) %>%
    pivot_wider(names_from = quantile,values_from = value) %>%
    left_join(node_isoweek_dates%>% mutate(isoweek=as.character(isoweek)),by=c('isoweek')) %>%
    left_join(Rc_agg,by=c('date_max'='date'))
  
  plot <- r_plot_data %>% group_by(parameter) %>% filter(date_max<max_date&date_max>min_date) %>%
    ggplot(aes(x=date_mid,y=mean,col=parameter,group=parameter)) + geom_line(lwd=1) + 
    geom_ribbon(aes(x=date_mid,ymin=`2.5%`, ymax=`97.5%`, fill = parameter,  
                    color = parameter), alpha=0.2, linetype = "dotted") +
    geom_line(aes(y=Rc,col='Rc aggregate'), linetype='dashed') +
    geom_ribbon(aes(x=date_mid,ymin=Rc_lb,ymax=Rc_ub,col='Rc aggregate',fill='Rc aggregate'),alpha=0.2, linetype='dashed') +
    geom_point(aes(x=date_mid,y=outdegree_average),shape=8) + 
    scale_colour_manual("",breaks = c("R",'Rc aggregate', 'R sim data'), values=c('red4','lightblue3','darkgreen')) +
    scale_fill_manual("",breaks = c("R",'Rc aggregate', 'R sim data'), values=c('red4','lightblue3','darkgreen')) +
    ylim(c(0,2.5)) + theme_light() + ylab('') + xlab('') + ggtitle('Reproduction Number Estimates')
  
  if(!is.null(sim_data))
  {
    sim_data <- sim_data %>% filter(isoweek!='all') %>%
      mutate(isoweek=as.double(isoweek)) %>%
      left_join(node_isoweek_dates,by=c('isoweek')) 
    
    plot <- plot + geom_line(data=sim_data %>% filter(parameter=='R'),aes(y=mean,col='R sim data'), lwd=1) +
      geom_ribbon(aes(x=date_mid,ymin=`2.5%`, ymax=`97.5%`,col='R sim data',fill='R sim data'),alpha=0.2, linetype='dashed',data=sim_data%>% filter(parameter=='R')) 
  }
  # 2 show variation across random trees
  
  # 
  # if(!is.null(Rc_agg))
  # {
  #   plot <- Rc_data %>% ggplot(aes(x=date_max,y=Estimate+Intercept,col='Rc_ind')) + geom_line() +
  #     geom_ribbon(aes(x=date_max,ymin=Q2.5+Intercept_Q2.5, ymax=Q97.5+Intercept_Q97.5), 
  #                 alpha=0.2, linetype = "dotted") +
  #     geom_line(aes(x=date_max,y=Rc, col='Rc_agg')) + 
  #     theme_light()
  # } else {
  #   plot <- Rc_data %>% ggplot(aes(x=date_max,y=Estimate+Intercept,col='Rc_ind')) + geom_line() +
  #     geom_ribbon(aes(x=date_max,ymin=Q2.5+Intercept_Q2.5, ymax=Q97.5+Intercept_Q97.5), 
  #                 alpha=0.2, linetype = "dotted") +
  #     theme_light()
  # }

  return(plot)
}

# plot Rc for random trees to compare with ML estimate
plot_Rc_random_trees <- function(model_in='negative_binomial',node_isoweek_dates,Rc_agg=NULL,max_date=as.Date('2021-08-01'),sim_data=NULL) 
{
  # plot population level rather an output of model.
  k_results            <- read_csv('./stan_output/Rc_param_estimates_combined_ML_trees.csv')
  k_results_random     <- read_csv('data/Rc_param_estimates_combined_random_trees.csv') %>% unique()

  # 1 for ML tree
  r_plot_data <-k_results_random %>% filter(isoweek!='all' & model == model_in) %>% filter(parameter!='lp__') %>% 
    dplyr::select(parameter,mean,`2.5%`,`97.5%`,isoweek,random_id) %>%
    group_by(isoweek,random_id)  %>%
    pivot_wider(names_from = parameter, values_from = c(mean,`2.5%`,`97.5%`)) %>%
    pivot_longer(
      -c(isoweek,random_id),
      names_to = c("quantile", "parameter"),
      names_pattern = "(.*)_(.*)",
      values_to = "value"
    ) %>% 
    filter(parameter %in% c('R')) %>%
    pivot_wider(names_from = quantile,values_from = value) %>%
    left_join(node_isoweek_dates%>% mutate(isoweek=as.character(isoweek)),by=c('isoweek')) %>%
    left_join(Rc_agg,by=c('date_max'='date')) %>%
    left_join(tmp)
  
  plot <- r_plot_data %>% filter(date_max<max_date) %>% 
    ggplot(aes(x=date_mid,y=mean,col=random_id)) + geom_line() +
         #geom_ribbon(aes(x=date_max,ymin=`2.5%`, ymax=`97.5%`), 
         #             alpha=0.2, linetype = "dotted") +
         geom_line(aes(x=date_mid,y=Rc, col='green')) + ylim(c(0,2.25)) +
         theme_light() + theme(legend.position = 'none')
  
  plot <- r_plot_data %>% group_by(parameter) %>% filter(date_max<max_date, date_max > as.Date('2020-10-01')) %>%
    ggplot(aes(x=date_mid,y=mean,col=parameter,group=parameter)) + geom_line(lwd=1) + 
    geom_ribbon(aes(x=date_mid,ymin=`2.5%`, ymax=`97.5%`, fill = parameter,  
                    color = parameter), alpha=0.2, linetype = "dotted") +
    geom_line(aes(y=Rc,col='Rc aggregate'), linetype='dashed') +
    geom_ribbon(aes(x=date_mid,ymin=Rc_lb,ymax=Rc_ub,col='Rc aggregate',fill='Rc aggregate'),alpha=0.2, linetype='dashed') +
    #geom_point(aes(x=date_mid,y=outdegree_average),shape=8) + 
    scale_colour_manual("",breaks = c("R",'Rc aggregate', 'R sim data'), values=c('red4','lightblue3','darkgreen')) +
    scale_fill_manual("",breaks = c("R",'Rc aggregate', 'R sim data'), values=c('red4','lightblue3','darkgreen')) +
    ylim(c(0,2.5)) + theme_light() + ylab('') + xlab('') + ggtitle('Reproduction Number Estimates')
  
  if(!is.null(sim_data))
  {
    sim_data <- sim_data %>% filter(isoweek!='all') %>%
      mutate(isoweek=as.double(isoweek)) %>%
      left_join(node_isoweek_dates,by=c('isoweek')) 
    
    plot <- plot + geom_line(data=sim_data %>% filter(parameter=='R'),aes(y=mean,col='R sim data'), lwd=1) +
      geom_ribbon(aes(x=date_mid,ymin=`2.5%`, ymax=`97.5%`,col='R sim data',fill='R sim data'),alpha=0.2, linetype='dashed',data=sim_data%>% filter(parameter=='R')) 
  }
  # 2 show variation across random trees
  
  # 
  # if(!is.null(Rc_agg))
  # {
  #   plot <- Rc_data %>% ggplot(aes(x=date_max,y=Estimate+Intercept,col='Rc_ind')) + geom_line() +
  #     geom_ribbon(aes(x=date_max,ymin=Q2.5+Intercept_Q2.5, ymax=Q97.5+Intercept_Q97.5), 
  #                 alpha=0.2, linetype = "dotted") +
  #     geom_line(aes(x=date_max,y=Rc, col='Rc_agg')) + 
  #     theme_light()
  # } else {
  #   plot <- Rc_data %>% ggplot(aes(x=date_max,y=Estimate+Intercept,col='Rc_ind')) + geom_line() +
  #     geom_ribbon(aes(x=date_max,ymin=Q2.5+Intercept_Q2.5, ymax=Q97.5+Intercept_Q97.5), 
  #                 alpha=0.2, linetype = "dotted") +
  #     theme_light()
  # }
  
  return(plot)
}

plot_overdispersion_random_trees <- function(model_in='negative_binomial',node_isoweek_dates,Rc_agg=NULL,max_date=as.Date('2021-08-01'),sim_data=NULL) 
{
  # plot population level rather an output of model.
  k_results            <- read_csv('data/Rc_param_estimates_combined_ML_trees.csv')
  k_results_random     <- read_csv('data/Rc_param_estimates_combined_random_trees.csv') %>% unique()
  
  # 1 for ML tree
  k_plot_data <- k_results_random %>% filter(isoweek!='all' & model == model_in) %>% filter(parameter!='lp__') %>% 
    dplyr::select(parameter,mean,`2.5%`,`97.5%`,isoweek,random_id) %>%
    group_by(isoweek,random_id) %>%
    pivot_wider(names_from = parameter, values_from = c(mean,`2.5%`,`97.5%`)) %>%
    pivot_longer(
      -c(isoweek,random_id),
      names_to = c("quantile", "parameter"),
      names_pattern = "(.*)_(.*)",
      values_to = "value"
    ) %>% 
    filter(parameter %in% c('k')) %>%
    pivot_wider(names_from = quantile,values_from = value) %>%
    left_join(node_isoweek_dates%>% mutate(isoweek=as.character(isoweek)),by=c('isoweek')) %>%
    left_join(Rc_agg,by=c('date_max'='date'))
  
  plot <- k_plot_data %>% filter(date_max<max_date) %>% 
    ggplot(aes(x=date_mid,y=mean,col=random_id)) + geom_line() +
    #geom_ribbon(aes(x=date_max,ymin=`2.5%`, ymax=`97.5%`), 
    #             alpha=0.2, linetype = "dotted") +
    theme_light() + theme(legend.position = 'none')
  
  plot <- r_plot_data %>% group_by(parameter) %>% filter(date_max<max_date) %>%
    ggplot(aes(x=date_mid,y=mean,col=parameter,group=parameter)) + geom_line(lwd=1) + 
    geom_ribbon(aes(x=date_mid,ymin=`2.5%`, ymax=`97.5%`, fill = parameter,  
                    color = parameter), alpha=0.2, linetype = "dotted") +
    geom_line(aes(y=Rc,col='Rc aggregate'), linetype='dashed') +
    geom_ribbon(aes(x=date_mid,ymin=Rc_lb,ymax=Rc_ub,col='Rc aggregate',fill='Rc aggregate'),alpha=0.2, linetype='dashed') +
    geom_point(aes(x=date_mid,y=outdegree_average),shape=8) + 
    scale_colour_manual("",breaks = c("R",'Rc aggregate', 'R sim data'), values=c('red4','lightblue3','darkgreen')) +
    scale_fill_manual("",breaks = c("R",'Rc aggregate', 'R sim data'), values=c('red4','lightblue3','darkgreen')) +
    ylim(c(0,2.5)) + theme_light() + ylab('') + xlab('') + ggtitle('Reproduction Number Estimates')
  
  if(!is.null(sim_data))
  {
    sim_data <- sim_data %>% filter(isoweek!='all') %>%
      mutate(isoweek=as.double(isoweek)) %>%
      left_join(node_isoweek_dates,by=c('isoweek')) 
    
    plot <- plot + geom_line(data=sim_data %>% filter(parameter=='R'),aes(y=mean,col='R sim data'), lwd=1) +
      geom_ribbon(aes(x=date_mid,ymin=`2.5%`, ymax=`97.5%`,col='R sim data',fill='R sim data'),alpha=0.2, linetype='dashed',data=sim_data%>% filter(parameter=='R')) 
  }
  
  return(plot)
}

# Plot of overdispersion parameter
plot_overdispersion <- function(model_in='negative_binomial',node_isoweek_dates, sim_data=NULL)
{
  # plot population level rather an output of model.
  k_results <- read_csv('data/Rc_param_estimates.csv')
  
  # 1 for ML tree
  k_plot_data <- k_results %>% filter(random_id=='ML' & isoweek!='all' & model == model_in) %>% 
    filter(parameter!='lp__') %>% 
    dplyr::select(parameter,mean,`2.5%`,`97.5%`,isoweek) %>%
    pivot_wider(names_from = parameter, values_from = c(mean,`2.5%`,`97.5%`)) %>%
    pivot_longer(
      -isoweek,
      names_to = c("quantile", "parameter"),
      names_pattern = "(.*)_(.*)",
      values_to = "value"
    ) %>% 
    filter(parameter %in% c('k')) %>%
    pivot_wider(names_from = quantile,values_from = value) %>%
    left_join(node_isoweek_dates%>% mutate(isoweek=as.character(isoweek)),by=c('isoweek')) 
  
  plot <- k_plot_data %>% group_by(parameter) %>% filter(date_max<as.Date('2021-08-01')) %>%
    ggplot(aes(x=date_mid,y=mean,col=parameter,group=parameter)) + geom_line(lwd=1) + 
    geom_ribbon(aes(x=date_mid,ymin=`2.5%`, ymax=`97.5%`, fill = parameter,  
                    color = parameter), alpha=0.2, linetype = "dotted") +
    theme_light() + ylab('') + xlab('') + ggtitle('Overdispersion Estimates')
  
  if(!is.null(sim_data))
  {
    sim_data <- sim_data %>% filter(isoweek!='all') %>%
      mutate(isoweek=as.double(isoweek)) %>% 
      left_join(node_isoweek_dates,by=c('isoweek')) 
    
    plot <- plot + geom_line(data=sim_data %>% filter(parameter=='k'),aes(y=mean,col='k sim data'), lwd=1) +
      geom_ribbon(aes(x=date_mid,ymin=`2.5%`, ymax=`97.5%`,col='k sim data',fill='k sim data'),alpha=0.2, linetype='dashed',data=sim_data%>% filter(parameter=='k')) 
  }
  # 2 show variation across random trees
  
  # 
  # if(!is.null(Rc_agg))
  # {
  #   plot <- Rc_data %>% ggplot(aes(x=date_max,y=Estimate+Intercept,col='Rc_ind')) + geom_line() +
  #     geom_ribbon(aes(x=date_max,ymin=Q2.5+Intercept_Q2.5, ymax=Q97.5+Intercept_Q97.5), 
  #                 alpha=0.2, linetype = "dotted") +
  #     geom_line(aes(x=date_max,y=Rc, col='Rc_agg')) + 
  #     theme_light()
  # } else {
  #   plot <- Rc_data %>% ggplot(aes(x=date_max,y=Estimate+Intercept,col='Rc_ind')) + geom_line() +
  #     geom_ribbon(aes(x=date_max,ymin=Q2.5+Intercept_Q2.5, ymax=Q97.5+Intercept_Q97.5), 
  #                 alpha=0.2, linetype = "dotted") +
  #     theme_light()
  # }
  
  return(plot)
}

# Plot of NPI effect sizes on transmission
plot_NPI <- function(model)
{
  model_coefs             <- as_tibble(fixef(model))
  model_coefs$covariate   <- rownames(fixef(model))
  
  npi_data                <- model_coefs %>% filter(str_detect(covariate,'combined_numeric')) %>%
    mutate(covariate = str_replace(covariate,'_combined_numeric',''),
           NPI       = case_when(covariate == 'C1' ~ 'Schools Closing',
                                 covariate == 'C2' ~ 'Workplace Closure',
                                 covariate == 'C4' ~ 'Restrictions on gatherings',
                                 covariate == 'H6' ~ 'Facial Coverings',
                                 covariate == 'H8' ~ 'Protection of elderly people',
                                 TRUE ~ 'other NPI')) %>%
    mutate(col = case_when(Q97.5 < 0 ~ 'red',
                           Q2.5 > 0  ~ 'blue',
                           Q97.5 >0 & Q2.5 < 0 ~ 'black'))
    
  plot <- npi_data %>% ggplot(aes(y=NPI,x=Estimate,xmin=Q2.5,xmax=Q97.5,col=col)) +
    geom_pointinterval() + 
    geom_vline(aes(xintercept=0),linetype='dashed') +
    scale_color_manual(values = c('red'='red','blue'='blue','black'='black')) +
    theme_light() + theme(legend.position = 'none') 
  
  return(plot)
}

# Plot of vaccination effect sizes on transmission
plot_vaccination <- function(model)
{
  model_coefs             <- as_tibble(fixef(model))
  model_coefs$covariate   <- rownames(fixef(model))
  
  vacc_data                <- model_coefs %>% filter(str_detect(covariate,'vacc')) %>%
    mutate(age       = case_when(covariate == 'age_group0M10yrs' ~ 'Age 0-10 years',
                                 covariate == 'age_group11M18yrs' ~ 'Age 11-18 years',
                                 covariate == 'age_group19M39yrs' ~ 'Age 19-39 years',
                                 covariate == 'age_group40M59yrs' ~ 'Age 40-59 years',
                                 TRUE ~ 'Age over 60 years')) %>%
    mutate(col = case_when(Q97.5 < 0 ~ 'red',
                           Q2.5 > 0  ~ 'blue',
                           Q97.5 >0 & Q2.5 < 0 ~ 'black'))
  
  plot <- vacc_data %>% ggplot(aes(y=covariate,x=Estimate,xmin=Q2.5,xmax=Q97.5,col=col)) +
    geom_pointinterval() + 
    geom_vline(aes(xintercept=0),linetype='dashed') +
    scale_color_manual(values = c('red'='red','blue'='blue','black'='black')) +
    theme_light() + theme(legend.position = 'none') 
  
  return(plot)
}

# Plot effect size of personal characteristics (e.g. age)
plot_individual_attributes <- function(model)
{
  model_coefs             <- as_tibble(fixef(model))
  model_coefs$covariate   <- rownames(fixef(model))
  
  age_data                <- model_coefs %>% filter(str_detect(covariate,'age')) %>%
    mutate(age       = case_when(covariate == 'age_group0M10yrs' ~ 'Age 0-10 years',
                                 covariate == 'age_group11M18yrs' ~ 'Age 11-18 years',
                                 covariate == 'age_group19M39yrs' ~ 'Age 19-39 years',
                                 covariate == 'age_group40M59yrs' ~ 'Age 40-59 years',
                                 TRUE ~ 'Age over 60 years')) %>%
    mutate(col = case_when(Q97.5 < 0 ~ 'red',
                           Q2.5 > 0  ~ 'blue',
                           Q97.5 >0 & Q2.5 < 0 ~ 'black'))
  
  plot <- age_data %>% ggplot(aes(y=age,x=Estimate,xmin=Q2.5,xmax=Q97.5,col=col)) +
    geom_pointinterval() + 
    geom_vline(aes(xintercept=0),linetype='dashed') +
    scale_color_manual(values = c('red'='red','blue'='blue','black'='black')) +
    theme_light() + theme(legend.position = 'none') 
  
  return(plot)
}

