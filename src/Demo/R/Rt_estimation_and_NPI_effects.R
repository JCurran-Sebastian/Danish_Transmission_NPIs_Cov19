# Rt estimation and NPI effect plots
# Corresponds to Figure 4 (reproduction numbers and overdispersion) and
# Figure 5 (NPI effectiveness) in Curran-Sebastian et al.
# Set working directory to R/ before running.

library(tidyverse)
library(tidybayes)
library(rstan)
library(patchwork)
library(zoo)
library(MASS)
library(lme4)
library(jtools)


source('./stan_utility_fns.R')
source('./plotting_fns.R')
source('./partial_sequencing_adjustment.R')
options(mc.cores=parallel::detectCores())

# Zero-inflated negative binomial Stan model. ZI component accounts for
# individuals with zero secondary infections for structural reasons (e.g.
# very short follow-up) rather than low R. Fitted weekly to estimate
# time-varying R and overdispersion k (Figure 4A,B).
model_zi <- stan_model('./zero_inflated_neg_bin.stan')

# Functions ----------------------------------------------------------------------

# Wrapper around rstan::sampling that returns an empty tibble on failure.
# Used when looping over many trees/weeks where occasional divergences are
# expected and should not abort the loop.
run_nb_inference <- function(stan_model,data_list,isoweek_in,file,
                             k_results_random_tmp)
{
  tryCatch(
    expr = {
      fit       <- sampling(stan_model,data=data_list,iter=2000,warmup = 1000,
                            chains = 4, seed=123, refresh = 0)

      fit_summary           <- as_tibble(summary(fit)$summary)
      fit_summary$parameter <- rownames(summary(fit)$summary)
      fit_summary$isoweek   <- isoweek_in
      fit_summary$random_id <- str_remove(str_remove(file,'nodelist_'),'.csv')
      fit_summary$model     <- stan_model@model_name
    },
    error = function(e){
      fit_summary           <- k_results_random_tmp
      fit_summary$isoweek   <- isoweek_in
      fit_summary$random_id <- str_remove(str_remove(file,'nodelist_'),'.csv')
      fit_summary$model     <- stan_model@model_name
    },
    warning = function(w){
      fit_summary           <- k_results_random_tmp
      fit_summary$isoweek   <- isoweek_in
      fit_summary$random_id <- str_remove(str_remove(file,'nodelist_'),'.csv')
      fit_summary$model     <- stan_model@model_name
    },
    finally = {
      if(dim(fit_summary)[1]==0)
      {
        fit_summary           <- k_results_random_tmp
        fit_summary$isoweek   <- isoweek_in
        fit_summary$random_id <- str_remove(str_remove(file,'nodelist_'),'.csv')
        fit_summary$model     <- stan_model@model_name
      }

      return(fit_summary)
    }
  )
}

# Load node-level attributes from the maximum-likelihood (ML) transmission tree.
# out_degree = number of secondary infections caused by each individual (individual Rc).
# Week index spans 2020-2021 continuously: isoweek + (year-2020)*52.
node_attrs <- read_csv('../data/node_attrs.csv') %>%
  dplyr::select(-c('...1')) %>%
  mutate(isoweek = isoweek(SampleDate) + (year(SampleDate)-2020)*52 )

# Diagnostic: out-degree distribution by sample date
node_attrs %>% group_by(SampleDate) %>%
  ggplot(aes(x=as.factor(SampleDate),y=out_degree)) +
  geom_boxplot() + theme_light()

# Summarise date range within each ISO week for downstream joins
node_isoweek_dates <- node_attrs %>% group_by(isoweek) %>% summarise(date_min=min(SampleDate),
                                                                     date_max=max(SampleDate))
unique_weeks  <- sort(node_attrs$isoweek |> unique())

# Aggregate Rt and infection ascertainment rate (IAR) estimates from the
# national case-count model (Supplementary Section 4.1), used to compare
# individual-level Rc with aggregate Rt (Figure 4A).
rt_iar_file <- read_csv('../data/rt_iar.csv')


# Section 1: R and k estimates from the ML tree ----------------------------------
# Fit the ZINB model for all weeks pooled ('all') and for each ISO week
# separately. Weekly fitting captures temporal variation driven by variant
# succession, vaccination rollout and NPIs.
k_results    <- as_tibble(NULL)
for(stan_model in c(model_zi))
  for(isoweek_in in c('all',unique_weeks) )
  {
    if(isoweek_in == 'all')
    {
      Rc        <- node_attrs %>% pull(out_degree)
    } else {
      Rc        <- node_attrs %>% filter(isoweek==isoweek_in) %>%
        pull(out_degree)
    }

    data_list <- list(N=length(Rc),y=Rc)
    fit       <- sampling(stan_model,data=data_list,iter=2000,warmup = 1000,
                          chains = 4, seed=123)

    fit_summary           <- as_tibble(summary(fit)$summary)
    fit_summary$parameter <- rownames(summary(fit)$summary)
    fit_summary$isoweek   <- isoweek_in
    fit_summary$random_id <- 'ML'   # ML = maximum-likelihood infector tree
    fit_summary$model     <- stan_model@model_name
    k_results             <- bind_rows(k_results,fit_summary)
  }

write_csv(k_results,'stan_output/Rc_param_estimates.csv')


# Reshape to long format: one row per (isoweek, parameter) with mean and 95% CrI
data <- k_results %>% filter(random_id=='ML' & isoweek!='all') %>%
  filter(parameter!='lp__') %>%
  dplyr::select(parameter,mean,`2.5%`,`97.5%`,isoweek) %>%
  pivot_wider(names_from = parameter, values_from = c(mean,`2.5%`,`97.5%`)) %>%
  pivot_longer(
    -isoweek,
    names_to = c("quantile", "parameter"),
    names_pattern = "(.*)_(.*)",
    values_to = "value"
  ) %>%
  filter(parameter %in% c('k','R','theta')) %>%
  pivot_wider(names_from = quantile,values_from = value) %>%
  left_join(node_isoweek_dates%>% mutate(isoweek=as.character(isoweek)),
            by=c('isoweek')) %>%
  #left_join(npi_data,by=c('date_min'='Date')) %>%
  left_join(rt_iar_file,by=c('date_max'='date'))




# Reload from disk so this plotting section can be run without re-fitting
k_results <- read_csv('stan_output/Rc_param_estimates.csv')

# Plot weekly Rc from the ML tree (Figure 4A, individual-level track)
r_plot_data <- k_results %>% filter(random_id=='ML' & isoweek!='all') %>%
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
  left_join(node_isoweek_dates%>%
              mutate(isoweek=as.character(isoweek)),
            by=c('isoweek')) #%>%


plot <- r_plot_data %>% group_by(parameter) %>% mutate(date_mid=date_min + (date_max-date_min)/2) %>%
  ggplot(aes(x=date_mid,y=mean,col=parameter,group=parameter)) + geom_line(lwd=1) +
  geom_ribbon(aes(x=date_mid,ymin=`2.5%`, ymax=`97.5%`, fill = parameter,
                  color = parameter), alpha=0.2, linetype = "dotted") +
  scale_colour_manual("",breaks = c("R",'Rc aggregate', 'R sim data'),
                      values=c('red4','lightblue3','darkgreen')) +
  scale_fill_manual("",breaks = c("R",'Rc aggregate', 'R sim data'),
                    values=c('red4','lightblue3','darkgreen')) +
  ylim(c(0,2.5)) + theme_light() + ylab('') + xlab('') +
  ggtitle('Reproduction Number Estimates')



plot


# Section 2: R and k across 100 sampled transmission trees ----------------------
# Repeats the ZINB model over each prioritised-settings tree in trees_output/.
# Averaging across trees accounts for uncertainty in infector assignment
# (Methods: Sampling and Analysis of Transmission Trees).
# Pre-computed results are in stan_output/ -- this loop takes ~20 minutes.
folder <- '../trees_output/'
files  <- list.files(folder)

k_results_random <- as_tibble(NULL)
for(file in files)
{
  node_attrs_tmp <- read_csv(paste0(folder,file), quote="'") %>%
    dplyr::select(-c('...1')) %>%
    mutate(isoweek = isoweek(SampleDate) + (year(SampleDate)-2020)*52 )

  node_isoweek_dates_tmp <- node_attrs_tmp %>% group_by(isoweek) %>%
    summarise(date_min=min(SampleDate),
              date_max=max(SampleDate))
  unique_weeks <- sort(node_attrs_tmp$isoweek |> unique())

  for(stan_model in c(model_zi) )
    for(isoweek_in in c('all',unique_weeks) )
    {
      if(isoweek_in == 'all')
      {
        Rc        <- node_attrs_tmp %>% pull(out_degree)
      } else {
        Rc        <- node_attrs_tmp %>% filter(isoweek==isoweek_in) %>%
          pull(out_degree)
      }

      data_list <- list(N=length(Rc),y=Rc)
      fit       <- sampling(stan_model,data=data_list,iter=2000,
                            warmup = 1000,chains = 4, seed=123)

      fit_summary           <- as_tibble(summary(fit)$summary)
      fit_summary$parameter <- rownames(summary(fit)$summary)
      fit_summary$isoweek   <- isoweek_in
      fit_summary$random_id <- str_remove(str_remove(file,'nodelist_'),'.csv')
      fit_summary$model     <- stan_model@model_name
      k_results_random      <- bind_rows(k_results_random,fit_summary)
    }
}
write_csv(k_results_random,
          'stan_output/Rc_param_estimates_combined_prioritised_settings.csv')




# Reload from disk so plotting sections can run independently
k_results            <- read_csv('stan_output/Rc_param_estimates.csv')
k_results_random     <- read_csv('stan_output/Rc_param_estimates_combined_prioritised_settings.csv') %>% unique()

# Plot overdispersion parameter k across sampled trees (Figure 4B).
# Increasing k toward 1 indicates more homogeneous offspring distributions;
# values below 0.5 indicate substantial superspreading potential.
k_plot_data <- k_results_random %>% filter(isoweek!='all') %>% filter(parameter!='lp__') %>%
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
  left_join(node_isoweek_dates%>% mutate(isoweek=as.character(isoweek)),by=c('isoweek')) #%>%
  #left_join(Rc_agg,by=c('date_max'='date'))

plot_k <- k_plot_data  %>% mutate(date_mid=date_min + (date_max-date_min)/2) %>%
  ggplot(aes(x=date_mid,y=mean, color = parameter,group = random_id)) + geom_line() +
  geom_ribbon(aes(x=date_mid,ymin=`2.5%`, ymax=`97.5%`, fill = parameter,
                  color = parameter), alpha=0.2, linetype = "dotted") +
  scale_colour_manual("",breaks = c("k"), values=c('lightblue3')) +
  scale_fill_manual("",breaks = c("k"), values=c('lightblue3')) +
  theme_light() + theme(legend.position = 'none')
plot_k

# Plot R across sampled trees to show robustness to infector assignment
plot_r <- r_plot_data %>% group_by(parameter) %>%
  mutate(date_mid=date_min + (date_max-date_min)/2) %>%
  ggplot(aes(x=date_mid,y=mean,col=parameter,group=parameter)) + geom_line(lwd=1) +
  geom_ribbon(aes(x=date_mid,ymin=`2.5%`, ymax=`97.5%`, fill = parameter,
                  color = parameter), alpha=0.2, linetype = "dotted") +
  scale_colour_manual("",breaks = c("R",'Rc aggregate', 'R sim data'),
                      values=c('red4','lightblue3','darkgreen')) +
  scale_fill_manual("",breaks = c("R",'Rc aggregate', 'R sim data'),
                    values=c('red4','lightblue3','darkgreen')) +
  theme_light() + ylab('') + xlab('') +
  ggtitle('Reproduction Number Estimates')

plot_r

  # plot population level rather an output of model.
  k_results            <- read_csv('stan_output/Rc_param_estimates.csv')
  k_results_random     <- read_csv('stan_output/Rc_param_estimates_combined_prioritised_settings.csv') %>% unique()

  r_plot_data <-k_results_random %>% filter(isoweek!='all') %>% filter(parameter!='lp__') %>%
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
    left_join(node_isoweek_dates%>% mutate(isoweek=as.character(isoweek)),by=c('isoweek')) #%>%
    #left_join(Rc_agg,by=c('date_max'='date')) %>%
    #left_join(tmp)

  plot <- r_plot_data %>%
    mutate(date_mid=date_min + (date_max-date_min)/2) %>%
    ggplot(aes(x=date_mid,y=mean,col=random_id)) + geom_line() +
    #geom_line(aes(x=date_mid,y=Rc, col='green')) + ylim(c(0,2.25)) +
    theme_light() + theme(legend.position = 'none')
  plot <- r_plot_data %>% group_by(parameter) %>%
    mutate(date_mid=date_min + (date_max-date_min)/2)%>%
    ggplot(aes(x=date_mid,y=mean,col=parameter,group=random_id)) + geom_line(lwd=1) +
    geom_ribbon(aes(x=date_mid,ymin=`2.5%`, ymax=`97.5%`, fill = parameter,
                    color = parameter), alpha=0.2, linetype = "dotted") +
    #geom_line(aes(y=Rc,col='Rc aggregate'), linetype='dashed') +
    #geom_ribbon(aes(x=date_mid,ymin=Rc_lb,ymax=Rc_ub,col='Rc aggregate',fill='Rc aggregate'),alpha=0.2, linetype='dashed') +
    #geom_point(aes(x=date_mid,y=outdegree_average),shape=8) +
    scale_colour_manual("",breaks = c("R",'Rc aggregate', 'R sim data'), values=c('red4','lightblue3','darkgreen')) +
    scale_fill_manual("",breaks = c("R",'Rc aggregate', 'R sim data'), values=c('red4','lightblue3','darkgreen')) +
    theme_light() + ylab('') + xlab('') + ggtitle('Reproduction Number Estimates')
  plot


# Section 3: NPI and risk-factor effects ----------------------------------------
# Forest plots based on posterior summaries from the Zero-Inflated Beta-Binomial
# (ZIBB) model fit in Python (numpyro/JAX, Section 4 of demo_pipeline.ipynb).
# The model regresses individual out-degrees on NPI stringency indices from the
# Oxford Covid Government Response Tracker (C1, C2, C4, C6, C8, H6),
# vaccination status, age group, SARS-CoV-2 variant, and Danish region (random
# effect). Coefficients beta_i represent the change in secondary infections per
# unit change in NPI stringency. See paper Figure 5 and Methods.


  # Plotting helper functions ----------------------------------------------------

  # Builds a combined forest plot from the MCMC posterior summary CSV.
  # plot_type='reduced_main_text' shows NPI and risk-factor panels only (Figure 5).
  # include_C3 / exclude_H8 toggle whether 'Cancel public events' (C3) and
  # 'Protection of elderly' (H8) are displayed; both are excluded in the main
  # analysis as C3 was only active early in the study period.
  # return_components=TRUE returns a named list of ggplot objects so panels can
  # be assembled flexibly with patchwork.
  create_multi_model_plot <- function(res2,plot_type = NA, include_C3 = TRUE, exclude_H8 = FALSE, return_components=FALSE, nrow_in=1 )
  {
    res_coef <- res2 %>% filter(!str_starts(coef_name,'mu')) %>%
      separate_wider_delim(coef_name, '[', names = c('coef_name','coef_id'),too_few = "align_start") %>%
      mutate(coef_id = str_replace(coef_id,']',''),
             setting = get_setting(coef_name),
             NPI = get_NPI(coef_id,setting,coef_name, include_C3, exclude_H8 )) %>%
      filter( !(coef_name=='beta_age_group' |                # filter out variables which we have transformed
                  coef_name=='beta_vacc_status' |            # use _full instead.
                  coef_name=='beta_variant' |
                  str_starts(coef_name,'log_alpha')) ) %>%
      mutate(coef_id = ifelse(is.na(coef_id), 0, coef_id)) %>%
      mutate(group = get_group(coef_name,coef_id)) %>%
      mutate(model = str_wrap(model, width = 15))

    res_coef['Risk factor'] <- NA
    res_coef[str_detect(res_coef$coef_name,'age_group'),]$`Risk factor` <- 'Age'
    res_coef[str_detect(res_coef$coef_name,'vacc_status'),]$`Risk factor` <- 'Vaccination'
    res_coef[str_detect(res_coef$coef_name,'hurdle'),]$`Risk factor` <- 'Hurdle rate'
    res_coef[str_detect(res_coef$coef_name,'alpha'),]$`Risk factor` <- 'Overdispersion'
    res_coef[str_detect(res_coef$coef_name,'variant'),]$`Risk factor` <- 'Variant'
    res_coef[str_detect(res_coef$coef_name,'holiday'),]$`Risk factor` <- 'Holiday'

    res_coef[str_detect(res_coef$coef_name,'hurdle'),]$group <- unlist(lapply(res_coef[str_detect(res_coef$coef_name,'hurdle'),]$coef_name, FUN = function(x){str_split(x,'_')[[1]][2]}))
    res_coef[str_detect(res_coef$coef_name,'alpha'),]$group <- unlist(lapply(res_coef[str_detect(res_coef$coef_name,'alpha'),]$coef_name, FUN = function(x){str_split(x,'_')[[1]][2]}))
    res_coef[str_detect(res_coef$coef_name,'holiday'),]$group <- unlist(lapply(res_coef[str_detect(res_coef$coef_name,'holiday'),]$coef_name, FUN = function(x){str_split(x,'_')[[1]][2]}))

    res_coef <- res_coef %>% mutate(setting=case_when(setting=='address'~'Household',
                                                      setting=='other'~'Community',
                                                      setting=='family'~'Family (n/h)',
                                                      setting=='tt'~'Total Transmission',
                                                      TRUE~str_to_title(setting)),
                                    group=case_when(group=='address'~'Household',
                                                    group=='other'~'Community',
                                                    group=='family'~'Family (n/h)',
                                                    setting=='tt'~'Total Transmission',
                                                    TRUE~str_to_title(group)))

    # Transparency reduced for very wide CrIs to avoid visual clutter
    npi_plt <- res_coef %>% filter(!is.na(NPI)) %>%
      mutate(alpha_size=case_when(`hdi_97.5%`-`hdi_2.5%`>4 ~ 0.8,
                                  `hdi_97.5%`-`hdi_2.5%`>2 ~ 0.9, TRUE ~1)) %>%
      ggplot(aes(y=NPI,x=mean,xmin=`hdi_2.5%`,xmax=`hdi_97.5%`,
                 col=model,group=model,alpha=alpha_size)) +
      geom_pointinterval(position = position_dodge(width=0.2)) +
      geom_vline(aes(alpha = alpha_size, xintercept=0),linetype='dashed') +
      scale_color_lancet() + scale_alpha(guide = 'none') +
      facet_wrap(~setting,nrow=nrow_in,scales='free_x') + theme_bw() +
      xlab('Effectiveness &beta;<sub>i</sub> <br> (change in number of secondary infections) <br> for each NPI') + ylab('') +
      theme(axis.title.x = element_markdown(),
            plot.tag = element_text(face='bold'))

    design <- matrix(c(3,3,3,3,1,2,2,2),4,2)
    zi_overdispersion_plt <- res_coef %>% filter(is.na(NPI) & `Risk factor` %in% c('Hurdle rate','Overdispersion') ) %>%
      filter(mean<1000) %>%
      mutate(`Risk factor` = case_when(`Risk factor`=='Overdispersion'&group=='Household' ~ 'Overdispersion Beta Binomial',
                                       `Risk factor`=='Overdispersion'&group!='Household' ~ 'Overdispersion Negative Binomial',
                                       TRUE ~ `Risk factor`),
             `Risk factor` = str_replace(`Risk factor`,'Hurdle rate','Zero-inflation probability')) %>%
      group_by(`Risk factor`) %>%
      ggplot(aes(y=group,x=mean,xmin=`hdi_2.5%`,xmax=`hdi_97.5%`,
                 col=model, group=model)) +
      geom_pointinterval(aes(alpha = alpha_size),
                         position = position_dodge(width=0.4)) +
      geom_vline(aes(xintercept=1),linetype='dashed') +
      scale_color_lancet() +
      theme_light() +  ylab('') +
      facet_manual(~`Risk factor`, scales = 'free',design = design) + theme_bw() + xlab('') +
      theme(plot.tag = element_text(face='bold'))

    group_order <- c( NA,'Workplace', 'School', 'Household', 'Family (n/h)','Community',
                      "0-10yrs","11-18yrs","19-39yrs","40-59yrs",'>60yrs',
                      "Two Dose","One Dose",'Unvaccinated',
                      "Omicron",'Delta',"Eta","Alpha","Wildtype",
                      "Nordjydanmark", "Midjydanmark", "Syddanmark", "Hovedstaden", "Sjaelland")

    # Risk factor IRRs: exponentiate log-scale coefficients from the NB linear predictor
    risk_factor_data <- res_coef %>% filter(is.na(NPI) & !(`Risk factor` %in% c('Hurdle rate','Overdispersion', NA)) ) %>%
      mutate(mean      = exp(mean),
             `hdi_2.5%`  = exp(`hdi_2.5%`),
             `hdi_97.5%` = exp(`hdi_97.5%`)) %>%
      filter(!(`Risk factor`=='Variant' & mean  > 7) ) %>%  # filter out large value for omicron (as we have too little data for that case)
      filter(!(`Risk factor`=='Holiday' & `hdi_97.5%`  > 2.2 ) ) %>%
      filter(str_ends(coef_name,'_full')) %>%
      filter(!(str_starts(coef_name,'holiday')&coef_id==0)) %>%
      mutate(group=factor(group,levels=group_order))

    shaded_data <- risk_factor_data %>%
      filter(abs(mean-1)==0) %>%
      distinct(`Risk factor`, group)

    risk_factor_plt <- ggplot(risk_factor_data,aes(y=group,x=mean,xmin=`hdi_2.5%`,xmax=`hdi_97.5%`)) +
      geom_pointinterval(aes(col=model,group=model),position = position_dodge(width=0.2)) +
      geom_vline(aes(xintercept=1),linetype='dashed') +
      scale_color_lancet() +
      theme_light() + theme(legend.position = 'none') + ylab('') +
      facet_wrap(~`Risk factor`, scales = 'free',nrow=1) + theme_bw() + xlab('IRR') +
      theme(plot.tag = element_text(face='bold'))

    random_effects_plt <- res_coef %>% filter(str_starts(coef_name,'group_intercept')) %>%
      ggplot(aes(y=group,x=mean,xmin=`hdi_2.5%`,xmax=`hdi_97.5%`,col=model,group=model)) +
      geom_pointinterval(position = position_dodge(width=0.4)) +
      geom_vline(aes(xintercept=0),linetype='dashed') +
      scale_color_lancet() +
      theme_light() + theme(legend.position = 'none') +
      facet_wrap(~setting, nrow = 1, scales = 'free_x') + ylab('Regional Random Effects') + theme_bw() + xlab('') + ylab('') +
      theme(plot.tag = element_text(face='bold'))# + xlim(c(-1.25,1.25))

    if(return_components)
    {
      return(list(npi_plt               = npi_plt,
                  risk_factor_plt       = risk_factor_plt,
                  random_effects_plt    = random_effects_plt,
                  zi_overdispersion_plt = zi_overdispersion_plt))
    }

    if(is.na(plot_type)) {
      plt <- npi_plt / risk_factor_plt / random_effects_plt / zi_overdispersion_plt +
        plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

      return(plt)
    } else if( plot_type == 'reduced_main_text' ){
      plt <- npi_plt / risk_factor_plt +
        plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

      return(plt)
    }


  }

  # Extracts transmission setting from a coefficient name.
  # 'beta_school_...' -> 'school', 'group_intercept_address' -> 'address'.
  # Returns NA for individual-level covariates (age, vaccination, variant).
  get_setting <- function(x)
  {
    setting <- rep(NA, length(x))
    for(i in 1:length(x))
    {
      if((str_starts(x[i],'beta_')|str_starts(x[i],'mu_'))&!(str_detect(x[i],'age')|str_detect(x[i],'vacc')|str_detect(x[i],'variant')))

        setting[i] <- str_split(x[i],'_')[[1]][2]

      if((str_starts(x[i],'group_intercept')))

        setting[i] <- str_replace(x[i],'group_intercept_','')
    }

    return(setting)
  }

  # Creates a small annotation panel with a downward arrow and a transmission
  # percentage label, used to show what share of all transmissions occurred in
  # each setting (displayed between panels A and B in Figure 5).
  create_arrow_panel <- function(x)
  {
    arrow_panel_commuity <- ggplot() + annotate('segment',x=0.5,xend=0.5,y=0.5,yend=0,arrow=arrow(length=unit(0.2,'cm')), linewidth=1) +
      xlab('') + ylab('') +
      annotate('text', x=0.52,y=0.45,label= '% of transmission',size=3.5,fontface='italic') +
      annotate('text', x=0.52,y=0.25,label= x,size=3.5) +
      annotate('text', x=0.55,y=0.25,label= '',size=3.5) +
      theme_void() + theme(plot.tag = element_blank())
  }

  # Maps MCMC coefficient array indices to Oxford Covid Government Response
  # Tracker NPI labels. The index order matches the design matrix
  # X = [C1, C2, C4, C6, C8, H6] used in the ZIBB regression (Methods).
  # C3 (cancel public events) and H8 (protection of elderly) are included only
  # when the corresponding flags are set; both are excluded in the main analysis.
  get_NPI <- function(x,y,z, include_C3=TRUE, exclude_H8 = TRUE)
  {
    NPI_list <- list()

    if(include_C3 & !exclude_H8 )
    {
      NPI_list[['0']] <- 'School restrictions'
      NPI_list[['1']] <- 'Workplace restrictions'
      NPI_list[['2']] <- 'Cancel public events'
      NPI_list[['3']] <- 'Restriction on gathering size'
      NPI_list[['4']] <- 'Stay at home order guidance'
      NPI_list[['5']] <- 'Restrictions on international travel'
      NPI_list[['6']] <- 'Facial coverings'
      NPI_list[['7']] <- 'Protection of elderly people'
    } else if(include_C3 & exclude_H8 ){
      NPI_list[['0']] <- 'School restrictions'
      NPI_list[['1']] <- 'Workplace restrictions'
      NPI_list[['2']] <- 'Cancel public events'
      NPI_list[['3']] <- 'Restriction on gathering size'
      NPI_list[['4']] <- 'Stay at home order guidance'
      NPI_list[['5']] <- 'Restrictions on international travel'
      NPI_list[['6']] <- 'Facial coverings'
    }  else if(!include_C3 & exclude_H8 ){
      NPI_list[['0']] <- 'School restrictions'
      NPI_list[['1']] <- 'Workplace restrictions'
      NPI_list[['2']] <- 'Restriction on gathering size'
      NPI_list[['3']] <- 'Stay at home order guidance'
      NPI_list[['4']] <- 'Restrictions on international travel'
      NPI_list[['5']] <- 'Facial coverings'
    } else {
      NPI_list[['0']] <- 'School restrictions'
      NPI_list[['1']] <- 'Workplace restrictions'
      NPI_list[['2']] <- 'Restriction on gathering size'
      NPI_list[['3']] <- 'Stay at home order guidance'
      NPI_list[['4']] <- 'Restrictions on international travel'
      NPI_list[['5']] <- 'Facial coverings'
      NPI_list[['6']] <- 'Protection of elderly people'
    }


    NPIs <- rep(NA,length(x))

    for(i in 1:length(NPIs))
      if(!is.na(x[i])&!is.na(y[i])&!str_starts(z[i],'group_intercept'))
        NPIs[i] <- NPI_list[[x[i]]]

    return(NPIs)
  }

  # Maps coefficient array indices to human-readable labels for categorical
  # covariates. Reference categories: 60+ yrs (age), unvaccinated (vaccination),
  # Wildtype (variant), as specified in the ZIBB model priors (Methods).
  get_group <- function(name,id)
  {
    age_group <- vacc_status <- variant <- region <- list()

    age_group[['0']] <- '>60yrs'
    age_group[['1']] <- '40-59yrs'
    age_group[['2']] <- '19-39yrs'
    age_group[['3']] <- '11-18yrs'
    age_group[['4']] <- '0-10yrs'

    vacc_status[['0']] <- 'unvaccinated'
    vacc_status[['1']] <- 'One dose'
    vacc_status[['2']] <- 'Two dose'

    variant[['0']] <- 'Variant 1'

    region[['0']] <- 'Region 1'


    group <- rep(NA, length(name))

    for(i in 1:length(group))
    {
      if(str_detect(name[i],'age_group')){

        group[i] <- age_group[[id[i]]]}


      if(str_detect(name[i],'vacc_status'))
        group[i] <- vacc_status[[id[i]]]

      if(str_detect(name[i],'variant'))
        group[i] <- variant[[id[i]]]

      if(str_detect(name[i],'group_intercept'))
        group[i] <- region[[id[i]]]
    }

    return(group)
  }




  # Load NPI effect estimates --------------------------------------------------

  # Posterior summary CSV produced by demo_pipeline.ipynb Section 4 (ZIBB model).
  # 'total_transmission' files rename 'beta_school' -> 'beta_tt' so that the
  # total-transmission and setting-specific results share the same plotting code.
  models <- list(
    'Prioritised settings ZI-BB'='mcmc_object_zibb_total_transmission_re_summary.csv'
  )
  overall_res <- as_tibble(NULL)

  for(model in names(models))
  {
    tmp <- read_csv(models[[model]]) %>% rename(coef_name=`...1`)
    tmp$model <- model

    overall_res <- bind_rows(overall_res,tmp)

    if(str_detect(models[[model]],'total_transmission'))
      overall_res <- overall_res %>% mutate(coef_name = str_replace(coef_name,"beta_school","beta_tt"))
  }


  models_tt <- list(               # tt = Total Transmission
    'Prioritised settings ZI-BB'='mcmc_object_zibb_total_transmission_re_summary.csv'
  )

  overall_tt_res <- as_tibble(NULL)

  for(model in names(models_tt))
  {
    tmp <- read_csv(models_tt[[model]]) %>% rename(coef_name=`...1`)
    tmp$model <- model

    overall_tt_res <- bind_rows(overall_tt_res,tmp)

    if(str_detect(models[[model]],'total_transmission'))
      overall_tt_res <- overall_tt_res %>% mutate(coef_name = str_replace(coef_name,"beta_school","beta_tt"))
  }
  plt_settings <- create_multi_model_plot(overall_res,
                                          plot_type = 'reduced_main_text',
                                          include_C3 = FALSE,
                                          exclude_H8 = TRUE,
                                          return_components = TRUE,
                                          nrow_in = 1
                                          )
  plt_tt       <- create_multi_model_plot(overall_tt_res,
                                          plot_type = 'reduced_main_text',
                                          include_C3 = FALSE,
                                          exclude_H8 = TRUE,
                                          return_components = TRUE,
                                          nrow_in = 3)

  # patchwork layout: row 'a' = total-transmission NPI panel, row 'b' = setting-
  # specific NPI panel, row 'c' = risk factors. Rows 'd'-'h' = transmission-
  # proportion arrow annotations between panels A and B.
  design_main <- '
   aaaaa
   aaaaa
   aaaaa
   aaaaa
   defgh
   bbbbb
   bbbbb
   bbbbb
   bbbbb
   ccccc
   ccccc
   ccccc
   ccccc'

  # Compute the fraction of all transmissions occurring in each setting from
  # the ML tree, used to annotate the composite figure (Figure 5 between A and B).
  data_ps_ml_ms <- read_csv('data_ps_ml_ms.csv')

  ps_ml_ms_tmp <- c(sum(data_ps_ml_ms$out_degree_family,na.rm = TRUE),
                    sum(data_ps_ml_ms$out_degree_address,na.rm = TRUE),
                    sum(data_ps_ml_ms$out_degree_school,na.rm = TRUE),
                    sum(data_ps_ml_ms$out_degree_workplace,na.rm = TRUE) ) /
    sum(data_ps_ml_ms$out_degree,na.rm = TRUE) * 100
  ps_ml_ms <- paste0(round(c(100-sum(ps_ml_ms_tmp),ps_ml_ms_tmp),1),'%')
  arrow_panel_commuity  <- create_arrow_panel(ps_ml_ms[1])
  arrow_panel_family_nh <- create_arrow_panel(ps_ml_ms[2])
  arrow_panel_household <- create_arrow_panel(ps_ml_ms[3])
  arrow_panel_school    <- create_arrow_panel(ps_ml_ms[4])
  arrow_panel_workplace <- create_arrow_panel(ps_ml_ms[5])

  plt <- plt_tt$npi_plt + # Some NPIs missing
    plt_settings$npi_plt + # Aren't these the same??
    plt_settings$risk_factor_plt + # These have overlapping age groups as NA for some reason...
    arrow_panel_commuity + arrow_panel_family_nh + arrow_panel_household +
    arrow_panel_school + arrow_panel_workplace +
    plot_annotation(tag_levels = 'A') +
    plot_layout(guides = 'collect',design = design_main) &
    theme(legend.position = 'bottom')
  plt
