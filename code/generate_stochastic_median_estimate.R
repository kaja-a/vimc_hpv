# Generate median estimates from stochastic estimates and diagnostic plots

# ------------------------------------------------------------------------------
library (data.table)
library (tictoc)
library (patchwork)
library (ggpubr)
library (countrycode)


scenarios <- c ("hpv-routine-default",
                "hpv-routine-bestcase",
                "hpv-campaign-default",
                "hpv-campaign-bestcase", 
                "hpv-no-vaccination")

# loop through different scenarios
for (scenario in scenarios) {
  
  # read central estimate
  setwd ("D:/GitHub/vimc_hpv/output/16Feb2020")
  
  # set central estimate file
  central_file <- switch (
    scenario, 
    "hpv-routine-default"   = "central-burden-vaccination_201910gavi-4_hpv-routine-default.csv", 
    "hpv-routine-bestcase"  = "central-burden-vaccination_201910gavi-4_hpv-routine-bestcase.csv", 
    "hpv-campaign-default"  = "central-burden-vaccination_201910gavi-4_hpv-campaign-default.csv", 
    "hpv-campaign-bestcase" = "central-burden-vaccination_201910gavi-4_hpv-campaign-bestcase.csv", 
    "hpv-no-vaccination"    = "central-burden-novaccination_201910gavi-4_hpv-no-vaccination.csv"
    )
  
  # set stochastic estimate file
  stochastic_file <- switch (
    scenario, 
    "hpv-routine-default"   = "stochastic-burden-vaccination_201910gavi-4_hpv-routine-default.csv", 
    "hpv-routine-bestcase"  = "stochastic-burden-vaccination_201910gavi-4_hpv-routine-bestcase.csv", 
    "hpv-campaign-default"  = "stochastic-burden-vaccination_201910gavi-4_hpv-campaign-default.csv", 
    "hpv-campaign-bestcase" = "stochastic-burden-vaccination_201910gavi-4_hpv-campaign-bestcase.csv", 
    "hpv-no-vaccination"    = "stochastic-burden-novaccination_201910gavi-4_hpv-no-vaccination.csv"
  )
                          
  # central_file <- "central-burden-vaccination_201910gavi-4_hpv-routine-default.csv"
  
  # read central estimate
  central <- fread (central_file)
  # scenario <- "hpv-routine-default"
  
    # read stochastic estimate
  setwd ("D:/GitHub/vimc_hpv/output_psa/16Feb2020")
  
  psa_runs      <- 200   # number of runs for psa
  countries_n   <- 112   # number of countries (should be detected from central estimate file)
  # countries_n <- 1
  chunk <- 112 / 8 # number of countries in a chunk (14)
  # chunk <- 2 # DEBUG
  
  tic ()
  # read stochastic estimates file
  # stochastic_file <- "stochastic-burden-vaccination_201910gavi-4_hpv-routine-default.csv"
  header_row <- fread (stochastic_file, 
                       nrows = 0, 
                       header = T)
  toc ()
  
  rows_per_country_run <- 7401 # this should be estimated from the central file
  rows_per_country <- rows_per_country_run * psa_runs
  
  pdf (paste0 ("plots/plots_", scenario, ".pdf"))
  
  dat_median_full <- data.table ("disease" = character(), 
                                 "year" = numeric(), 
                                 "age" = numeric(), 
                                 "country" = character(), 
                                 "country_name" = character(), 
                                 "cases.median" = numeric (), 
                                 "cases.ci_low" = numeric (),   
                                 "cases.ci_high" = numeric (), 
                                 "deaths.median" = numeric (), 
                                 "deaths.ci_low" = numeric (), 
                                 "deaths.ci_high" = numeric (), 
                                 "dalys.median" = numeric (),
                                 "dalys.ci_low" = numeric (),
                                 "dalys.ci_high" = numeric (), 
                                 "rows" = numeric (), 
                                 "cohort_size" = numeric (), 
                                 "cases" = numeric (), 
                                 "dalys" = numeric (), 
                                 "deaths" = numeric ()
  )
  
  
  
  
  for (country_chunk in 0:7)  {  
  # for (country_chunk in 0:0)  { # DEBUG
    
    print (country_chunk)
    
    tic ()
    
    chunk_dat <- fread (stochastic_file, 
                        skip = 1 + ( (country_chunk) * rows_per_country * chunk), 
                        nrows = rows_per_country * chunk, 
                        header = F)
    
    colnames (chunk_dat) <- colnames (header_row)
    
    print (unique (chunk_dat [, country]) )
    
    for (country_code in unique (chunk_dat [, country])) {
      
      print (country_code)
      
      dat_central <- central [country == country_code]
      
      dat <- chunk_dat [country == country_code]
      
      dat_median <- dat [, list (cases.median   = median   (cases), 
                                 cases.ci_low   = quantile (cases, 0.025), 
                                 cases.ci_high  = quantile (cases, 0.975), 
                                 deaths.median  = median   (deaths), 
                                 deaths.ci_low  = quantile (deaths, 0.025), 
                                 deaths.ci_high = quantile (deaths, 0.975), 
                                 dalys.median   = median   (dalys), 
                                 dalys.ci_low   = quantile (dalys, 0.025), 
                                 dalys.ci_high  = quantile (dalys, 0.975) 
      ), 
      by = list (disease, year, age, country, country_name)]
      
      dat_median [, rows := as.numeric (rownames(dat_median))]
      
      dat_median <- dat_median [dat_central, on = .(disease = disease,
                                                    year = year, 
                                                    age = age, 
                                                    country = country,
                                                    country_name = country_name)]
      
      # combine country specific data to full data table
      dat_median_full <- rbindlist (list (dat_median_full, dat_median), 
                                    use.names = TRUE)
      
      # cases
      plot_cases <- ggplot (dat_median, aes (x = rows)) + 
        geom_point (aes (y = cases.median / cases.median, colour = "median / median"), size = 1, alpha = 0.5) + 
        geom_point (aes (y = cases.ci_low / cases.median, colour = "95ci_low / median"), size = 1, alpha = 0.5) + 
        geom_point (aes (y = cases.ci_high / cases.median, colour = "95ci_high / median"), size = 1, alpha = 0.5) + 
        geom_point (aes (y = cases / cases.median, colour = "mean / median"), size = 1, alpha = 0.5) + 
        labs (
          x= "rows (cases)",
          y= "cases ratio", 
          title = "cases ratio ", 
          subtitle = "cases (low_ci, mean, median, high_ci) / cases (median)", 
          colour = "cases ratio") + 
        ylim (0, NA)
      
      # deaths
      plot_deaths <- ggplot (dat_median, aes (x = rows)) + 
        geom_point (aes (y=deaths.median / deaths.median, colour = "median / median"), size = 1, alpha = 1) + 
        geom_point (aes (y=deaths.ci_low / deaths.median, colour = "95ci_low / median"), size = 1, alpha = 1) + 
        geom_point (aes (y=deaths.ci_high / deaths.median, colour = "95ci_high / median"), size = 1, alpha = 1) + 
        geom_point (aes (y=deaths / deaths.median, colour = "mean / median"), size = 1, alpha = 1) + 
        labs (
          x= "rows (deaths)",
          y= "deaths ratio", 
          title = "deaths ratio ", 
          subtitle = "deaths (low_ci, mean, median, high_ci) / deaths (median)", 
          colour = "deaths ratio") + 
        ylim (0, NA)
      
      # dalys
      plot_dalys <- ggplot (dat_median, aes (x = rows)) + 
        geom_point (aes (y=dalys.median / dalys.median, colour = "median / median"), size = 1, alpha = 1) + 
        geom_point (aes (y=dalys.ci_low / dalys.median, colour = "95ci_low / median"), size = 1, alpha = 1)  + 
        geom_point (aes (y=dalys.ci_high / dalys.median, colour = "95ci_high / median"), size = 1, alpha = 1) + 
        geom_point (aes (y=dalys / dalys.median, colour = "mean / median"), size = 1, alpha = 1) + 
        labs (
          x= "rows (dalys)",
          y= "dalys ratio", 
          title = "dalys ratio ", 
          subtitle = "dalys (low_ci, mean, median, high_ci) / dalys (median)", 
          colour = "dalys ratio") + 
        ylim (0, NA)
      
      # combine plots
      combined_plot <- list (plot_cases, plot_deaths, plot_dalys)
      
      # arrange plots in a single page
      q <- ggarrange (plotlist = combined_plot, nrow = 3, ncol = 1)
      
      print (
        annotate_figure (q, 
                         top = text_grob (paste0 (country_code, 
                                                  " / ", countrycode (country_code, 
                                                                      "iso3c", 
                                                                      "country.name"), 
                                                  " / ", scenario), 
                                          color = "black" 
                                          # , size = 9
                         )))
      
      
      
    }
    
    toc ()
  }
  
  # keep requiste columns and save disease burden (median estimates) in vimc format
  dat_median_vimc <- dat_median_full [, .(disease, year, age, country, country_name, 
                                          cohort_size, 
                                          cases = cases.median, 
                                          deaths = deaths.median, 
                                          dalys = dalys.median ) ]
  
  central_median_file <- paste0 ("median/", 
                                 substr (central_file, 
                                         start = 1, 
                                         stop = str_length (central_file) - str_length (".csv")), 
                                 "_median.csv")
  
  fwrite (x = dat_median_vimc, 
          file = central_median_file)
  
  
  
  dev.off ()
  
}

