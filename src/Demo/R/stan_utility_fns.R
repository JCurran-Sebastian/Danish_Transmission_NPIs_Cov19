## stan utility functions

loo.stanfit <-
  function(x,
           pars = "log_lik",
           ...,
           save_psis = FALSE,
           cores = getOption("mc.cores", 1)) {
    stopifnot(length(pars) == 1L)
    LLarray <- loo::extract_log_lik(stanfit = x,
                                    parameter_name = pars,
                                    merge_chains = FALSE)
    r_eff <- loo::relative_eff(x = exp(LLarray), cores = cores)
    loo::loo.array(LLarray,
                   r_eff = r_eff,
                   cores = cores,
                   save_psis = save_psis)
  }

run_stan_for_R_and_k <- function(reg_data, stan_model, out_degree_variable = 'out_degree')
{
  k_results    <- as_tibble(NULL)
  unique_weeks <- unique(reg_data$isoweek)
  
  for(isoweek_in in c('all', as.vector(unique_weeks)) )
  {
    if(isoweek_in == 'all')
    {
      Rc        <- reg_data %>% pull(!!out_degree_variable)
    } else {
      Rc        <- reg_data %>% filter(isoweek==as.double(isoweek_in)) %>% pull(!!out_degree_variable)
    }
    
    data_list <- list(N=length(Rc),y=Rc)
    fit       <- sampling(stan_model,data=data_list,iter=2000,warmup = 1000,chains = 4, seed=123)
    
    fit_summary           <- as_tibble(summary(fit)$summary)
    fit_summary$parameter <- rownames(summary(fit)$summary)
    fit_summary$isoweek   <- isoweek_in
    fit_summary$random_id <- 'ML'
    fit_summary$model     <- stan_model@model_name
    k_results             <- bind_rows(k_results,fit_summary)
  }  
  
  return(k_results)
}