# Filename: vimc_hpv_run.R

# run program from source file location/folder

# VIMC HPV runs 2020
# Central runs, stochastic runs, and diagnostic plots

#-------------------------------------------------------------------------------
# load libraries 
library (data.table)
library (stringr)
library (tictoc)
library (foreach)
library (doParallel)
library (prime)

rm (list = ls ())  # clear workspace
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# functions
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Generate vaccine impact estimates (central run)
EstimateVaccineImpactVimcCentral <- function (vaccine_coverage_file,
                                              disease_burden_template_file,
                                              disease_burden_no_vaccination_file, 
                                              disease_burden_vaccination_file, 
                                              disease_burden_results_file, 
                                              campaign_vaccination, 
                                              routine_vaccination) {

  # read files -- vaccination coverage
  vimc_coverage <- fread (vaccine_coverage_file)
  
  # read file -- central disease burden template
  vimc_template <- fread (disease_burden_template_file)
  vimc_template <- vimc_template [country == "CHN"] # DEGBUG -- remove this line later
  
  # register batch data for vimc runs
  RegisterBatchDataVimc (vimc_coverage             = vimc_coverage, 
                         vimc_template             = vimc_template, 
                         use_campaigns             = campaign_vaccination, 
                         use_routine               = routine_vaccination,
                         restrict_to_coverage_data = FALSE, 
                         force                     = TRUE, 
                         psa                       = 0)
  
  # log file to keep track of simulation run
  # log_file <- paste0 (paste0 ("log/prime_log.log", 
  #                             disease_burden_results_file)) 
  log_file <- "log/prime_log.log"
    
  # start of parallelisation
  cl <- makeCluster (detectCores())   # registering number of cores
  registerDoParallel (cl)             # start of parallelisation
  
  # simulation runs through the batch of cohorts
  results <- BatchRun(countries                       = -1,
                      coverage                        = -1,
                      agevac                          = -1,
                      agecohort                       = -1,
                      sens                            = -1,
                      year_born                       = -1,
                      year_vac                        = -1,
                      runs                            = 1,
                      vaccine_efficacy_beforesexdebut = 1,
                      vaccine_efficacy_aftersexdebut  = 0,
                      log                             = log_file,
                      by_calendaryear                 = TRUE,
                      use_proportions                 = TRUE,
                      analyseCosts                    = FALSE,
                      psa                             = 0,
                      psa_vals                        = ".data.batch.psa",
                      unwpp_mortality                 = TRUE,
                      disability.weights              = "gbd_2017",
                      canc.inc                        = "2018",
                      vaccine                         = "4vHPV"
  )
  
  # end of parallelisation
  stopCluster (cl)    
  
  # convert results to vimc format
  convert_results <- OutputVimc (DT            = results, 
                                 calendar_year = TRUE, 
                                 vimc_template = vimc_template)
  
  # save full results
  fwrite (x    = results, 
          file = disease_burden_results_file)
  
  # Saving output for no vaccination scenario (vimc format)
  no_vaccination <- convert_results [scenario == "pre-vaccination"]
  no_vaccination <- no_vaccination  [, colnames (vimc_template), with=F]
  no_vaccination <- no_vaccination  [!is.na(deaths)]
  fwrite (x    = no_vaccination, 
          file = disease_burden_no_vaccination_file)
  
  # Saving output for vaccination scenario (vimc format)
  vaccination <- convert_results [scenario == "post-vaccination"]
  vaccination <- vaccination     [, colnames (vimc_template), with=F]
  vaccination <- vaccination     [!is.na(deaths)]
  fwrite (x    = vaccination, 
          file = disease_burden_vaccination_file)
  
  return ()

} # end of function -- EstimateVaccineImpactVimcCentral
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# CreatePsaData function will be inserted here

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Generate vaccine impact estimates (sensitivity analysis/stochastic runs)
EstimateVaccineImpactVimcStochastic <- function (diseaseBurdenCentralFile,
                                                 resultsFile,
                                                 psaDataFile,
                                                 countryCode,
                                                 diseaseBurdenStochasticFile,
                                                 fileNumber) {

  # 112 csv files with stochastic estimates; 1 file per country


} # end of function -- EstimateVaccineImpactVimcStochastic
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Diagnostic plots of vaccine impact estimates
# (central run) and (sensitivity analysis/stochastic runs)
PlotVaccineImpactVimc <- function (diseaseBurdenCentralFile,
                                   countryCode,
                                   diseaseBurdenStochasticFile,
                                   fileNumber) {

  # 1 pdf file per scenario containing all 112 countries

} # end of function -- PlotVaccineImpactVimc
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# main program - start 
#-------------------------------------------------------------------------------

tic ()  # start timer

# save current folder name and move to base folder
current_folder <- getwd ()
current_folder <- setwd ("D:/GitHub/vimc_hpv/code")  # DEBUG -- remove this line later
setwd ("../")
#-------------------------------------------------------------------------------

# montangu touchstone
touchstone <- "201910gavi-4"

#-------------------------------------------------------------------------------
# Generate vaccine impact estimates (central run)

# central burden template file
central_burden_template_file <- paste0 ("input/central-burden-template.", 
                                        touchstone, 
                                        ".HPV_LSHTM-Jit_standard.csv")

# scenarios
scenarios <- c ("hpv-routine-default",
                "hpv-routine-bestcase",
                "hpv-campaign-default",
                "hpv-campaign-bestcase")
# scenarios <- c ("hpv-routine-default")  # DEGBUG -- remove this line later

# loop through different scenarios
for (scenario in scenarios) {
  
  # vaccination coverage filename
  # for routine vaccination: filename should contain "routine" 
  # for campaign vaccination: filename should contain "campaign"
  # filename should not contain both routine and campaign
  vaccine_coverage_file <- paste0 ("input/coverage_", 
                                   touchstone, "_", 
                                   scenario, 
                                   ".csv")
  
  # routine scenario indicates (routine + campaign) vaccination
  # campaign scenario indicates (campaign) vaccination
  campaign <- TRUE
  if (grepl (pattern = "routine", 
             x       = vaccine_coverage_file, 
             fixed   = TRUE)) {
    
    routine  <- TRUE
    
  } else if (grepl (pattern = "campaign", 
                    x       = vaccine_coverage_file, 
                    fixed   = TRUE)) {
    
    routine  <- FALSE
    
    # use routine vaccination coverage file (since is also includes campaign vaccination)
    # this is done this way since RegisterBatchDataVimc (prime) is implemented in this way
    vaccine_coverage_file <- str_replace (string      = vaccine_coverage_file, 
                                          pattern     = "campaign", 
                                          replacement = "routine")
  } else {
    
    stop (paste0 ("vaccination coverage file does not correspond to routine or campaign vaccination: ", 
                  vaccine_coverage_file))
  }
  
  # disease burden results files -- vimc format
  # includes cohort size, cases, deaths, dalys
  central_burden_no_vaccination_file <-  paste0 ("output/central-burden-novaccination_", 
                                                 touchstone, 
                                                 "_", scenario, ".csv")
  
  central_burden_vaccination_file <-  paste0 ("output/central-burden-vaccination_", 
                                              touchstone, 
                                              "_", scenario, ".csv")
  
  # disease burden results files -- internal format 
  # includes cohort size, cases, deaths, dalys, ylds, ylls
  # includes burden results for vaccination and no vaccination 
  central_burden_results_file <-  paste0 ("output/central-burden-results_", 
                                              touchstone, 
                                              "_", scenario, ".csv")
  
  print (paste0 ("burden template:", central_burden_template_file))
  print (paste0 ("scenario: ", scenario))
  print (paste0 ("coverage: ", vaccine_coverage_file))
  
  EstimateVaccineImpactVimcCentral (
    vaccine_coverage_file              = vaccine_coverage_file,
    disease_burden_template_file       = central_burden_template_file,
    disease_burden_no_vaccination_file = central_burden_no_vaccination_file,
    disease_burden_vaccination_file    = central_burden_vaccination_file,
    disease_burden_results_file        = central_burden_results_file,
    routine_vaccination                = routine,
    campaign_vaccination               = campaign
  )
  
}




# 
# # scenario: hpv-campaign-bestcase (campaign)
# EstimateVaccineImpactVimcCentral (
#   vaccine_coverage_file              = "input/coverage_201910gavi-4_hpv-campaign-bestcase.csv",
#   disease_burden_template_file       = central_burden_template_file,
#   disease_burden_no_vaccination_file = "output/central-burden-novaccination.201910gavi-4_hpv-campaign-bestcase.csv",
#   disease_burden_vaccination_file    = "output/central-burden-vaccination.201910gavi-4_hpv-campaign-bestcase.csv",
#   results_file                       = "output/central-burden-results.201910gavi-4_hpv-campaign-bestcase.csv",
#   routine                            = FALSE,
#   campaign                           = TRUE
# )
# 
# # scenario: hpv-routine-default (routine + campaign)
# EstimateVaccineImpactVimcCentral (
#   vaccine_coverage_file              = "input/coverage_201910gavi-4_hpv-routine-default.csv",
#   disease_burden_template_file       = central_burden_template_file,
#   disease_burden_no_vaccination_file = "output/central-burden-novaccination.201910gavi-4_hpv-routine-default.csv",
#   disease_burden_vaccination_file    = "output/central-burden-vaccination.201910gavi-4_hpv-routine-default.csv",
#   results_file                       = "output/central-burden-results.201910gavi-4_hpv-routine-default.csv",
#   routine                            = TRUE,
#   campaign                           = TRUE
# )
# 
# # scenario: hpv-campaign-default (campaign)
# EstimateVaccineImpactVimcCentral (
#   vaccine_coverage_file              = "input/coverage_201910gavi-4_hpv-campaign-default.csv",
#   disease_burden_template_file       = central_burden_template_file,
#   disease_burden_no_vaccination_file = "output/central-burden-novaccination.201910gavi-4_hpv-campaign-default.csv",
#   disease_burden_vaccination_file    = "output/central-burden-vaccination.201910gavi-4_hpv-campaign-default.csv",
#   results_file                       = "output/central-burden-results.201910gavi-4_hpv-campaign-default.csv",
#   routine                            = FALSE,
#   campaign                           = TRUE
# )

#-------------------------------------------------------------------------------
# return to source file location/folder
setwd (current_folder)

toc ()  # note current timer and compute elapsed time
#-------------------------------------------------------------------------------
# main program - end 
#-------------------------------------------------------------------------------







