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

# rm (list = ls ())  # clear workspace (DEBUG / uncomment this line later)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# functions
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#' Generate vaccine impact estimates (VIMC central run)
#'
#' Genrate vaccine impact estimates for VIMC central runs. The inputs are 
#'   vaccine coverage and disease burden template files and outputs are 
#'   disease burden estimates (pre-vaccination and post-vaccination). 
#'
#' Three disease burden estimates are generated. 
#'   (i)   disease burden estimates for no vaccination (vimc format)
#'   (ii)  disease burden estimates for vaccination (vimc format)
#'   (iii) disease burden estimates for vaccination (pre- and post-vaccination)
#'   and includes YLDs and YLLs
#'
#' @param vaccine_coverage_file csv file (input), vaccine coverage data 
#'   (vimc format)
#' @param disease_burden_template_file csv file (input), disease burden template 
#'   (vimc format)
#' @param disease_burden_no_vaccination_file csv file (output), disease burden estimates
#'   for no vaccination (vimc format) 
#' @param disease_burden_vaccination_file csv file (output), disease burden estimates
#'   for vaccination (vimc format) 
#' @param disease_burden_results_file csv file (output), disease burden estimates
#'   pre-vaccination and post-vaccination
#' @param campaign_vaccination logical, indicates campaign vaccination 
#' @param routine_vaccination logical, indicates routine vaccination
#'
#' @return 
#' @export
#'
#' @examples EstimateVaccineImpactVimcCentral (
#'   vaccine_coverage_file = "coverage_hpv-routine-default.csv",
#'   disease_burden_template_file       = "central-burden-template.csv",
#'   disease_burden_no_vaccination_file = "central_burden_no_vaccination.csv",
#'   disease_burden_vaccination_file    = "central_burden_vaccination.csv",
#'   disease_burden_results_file        = "central_burden_results.csv",
#'   routine_vaccination                = TRUE,
#'   campaign_vaccination               = TRUE)

EstimateVaccineImpactVimcCentral <- function (vaccine_coverage_file,
                                              disease_burden_template_file,
                                              disease_burden_no_vaccination_file, 
                                              disease_burden_vaccination_file, 
                                              disease_burden_results_file, 
                                              campaign_vaccination, 
                                              routine_vaccination, 
                                              vaccine = "4vHPV") {

  # read files -- vaccination coverage
  vimc_coverage <- fread (vaccine_coverage_file)
  
  # read file -- central disease burden template
  vimc_template <- fread (disease_burden_template_file)
  vimc_template <- vimc_template [country == "CHN"] # DEBUG -- comment this line later
  
  # register batch data for vimc runs
  RegisterBatchDataVimc (vimc_coverage             = vimc_coverage, 
                         vimc_template             = vimc_template, 
                         use_campaigns             = campaign_vaccination, 
                         use_routine               = routine_vaccination,
                         restrict_to_coverage_data = FALSE, 
                         force                     = TRUE, 
                         psa                       = 0)
  
  # log file to keep track of simulation run
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
                      vaccine                         = vaccine
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
  no_vaccination <- no_vaccination  [, colnames (vimc_template), with=FALSE]
  no_vaccination <- no_vaccination  [!is.na(deaths)]
  fwrite (x    = no_vaccination, 
          file = disease_burden_no_vaccination_file)
  
  # Saving output for vaccination scenario (vimc format)
  vaccination <- convert_results [scenario == "post-vaccination"]
  vaccination <- vaccination     [, colnames (vimc_template), with=FALSE]
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
#' Generate vaccine impact estimates (VIMC stochastic runs)
#'
#' Generate vaccine impact estimates for VIMC stochastic runs/sensitivity 
#'   analysis. The inputs are central disease burden estimates, input 
#'   parameter distributions (latin hyper sampling), runs forsensitivity analysis,
#'   and filename for stochastic burden estimates. The outputs are stochastic 
#'   disease burden estimates (full results file plus a file per country). 
#'   
#'   Stochastic disease burden estimates are generated.
#'     (i)  full results files
#'            1 full results file for pre-vaccination (optional)
#'            1 full results file for post-vaccination
#'     (ii) 1 file per country for all runs
#'            1 file per country for all runs -- pre-vaccination (optional)
#'            1 file per country for all runs -- post-vaccination
#' 
#' @param disease_burden_template_file csv file (required), 
#'   central disease burden template (vimc format); add column with run_id to get
#'   stochastic disease burden template
#' @param centralBurdenResultsFile csv file (required), central disease burden 
#'   estimates pre and post vaccination
#' @param psaData data table (required), latin hyper cube sample of input parameters 
#' @param diseaseBurdenStochasticFolder character string (required),
#'   stochastic disease burden estimates folder (output)
#' @param diseaseBurdenStochasticFile character string (required), 
#'   stochastic disease burden estimates file(s) (output)
#' @param psa_runs integer (required), simulation runs for sensitivity analysis
#' @param countryCodes list (optional), If country codes are provided, 
#'   stochastic burden estimates are generated for these countries. 
#'   If set to -1, then stochastic burden estimates are generated for the 
#'   countries included in the central burden estimates. 
#' @param vaccination_scenario logical (required), generate stochastic burden 
#'   estimates for (vaccination) or (no vaccination) scenario
#' 
#' @return 
#' @export
#'  
EstimateVaccineImpactVimcStochastic <- function (disease_burden_template_file, 
                                                 centralBurdenResultsFile,
                                                 psaData,
                                                 diseaseBurdenStochasticFolder,
                                                 diseaseBurdenStochasticFile, 
                                                 psa_runs, 
                                                 countryCodes = -1, 
                                                 vaccination_scenario) {
  
  # read file -- central disease burden template
  vimc_template <- fread (disease_burden_template_file)
  
  # read in central burden results
  central_burden <- fread (centralBurdenResultsFile)
  
  # initialise empty table to save stochastic results in vimc format
  header <- vimc_template [0, ]                 # empty table with requisite 
                                                #   columns
  header [, run_id := numeric ()]               # add column for run id
  setcolorder (header, c("disease", "run_id"))  # reorder columns
  
  # save header to file for stochastic estimates of disease burden
  stochastic_file <- paste0 (diseaseBurdenStochasticFolder, 
                             diseaseBurdenStochasticFile, 
                             ".csv")
  fwrite (x      = header, 
          file   = stochastic_file, 
          append = FALSE)
  
  
  # extract burden estimates for pre- or post-vaccination
  if (vaccination_scenario) {
    
    # vaccination scenario
    central_burden <- central_burden [scenario == "post-vaccination"]
    
  } else {
    
    # no vaccination scenario
    central_burden <- central_burden [scenario == "pre-vaccination"]
  }
  
  # get country codes
  if (countryCodes == -1) {
    countryCodes <- unique (central_burden [, country])
  } 
  
  # generate stochastic burden estimates for each country
  for (country_code in countryCodes) {
    
    # get disease burden template for current country (of this loop)
    vimc_template_country <- vimc_template [country == country_code]
    
    # get central burden estimates for current country 
    central_burden_country <- central_burden [country == country_code]
    
    # get latin hyper cube sample of input parameters for curent country
    psadat_country <- psaData [country == country_code]
    
    # file to save stochastic burden estimates of current country 
    stochastic_file_country <- paste0 (diseaseBurdenStochasticFolder, 
                                       diseaseBurdenStochasticFile, 
                                       "_", country_code, 
                                       ".csv")
    
    # initialise header for stochastic burden estimates
    stochastic_burden_country <- header
    
    # generate post-hoc stochastic estimates for each run
    for (run_number in 1:psa_runs) {
      
      # initialise burden estimate to central burden estimates
      burden <- central_burden_country
      
      # ------------------------------------------------------------------------
      # probabilistic sensitivity analysis to generate stochastic estimates 
      
      # cases -- incidence
      burden [, inc.cecx := (inc.cecx * psadat_country [run_id == run_number, 
                                                        incidence_ratio]
                                      * psadat_country [run_id == run_number, 
                                                        hpv_distribution_ratio])]
      
      # deaths -- mortality 
      burden [, mort.cecx := (mort.cecx * psadat_country [run_id == run_number, 
                                                          mortality_ratio]
                                        * psadat_country [run_id == run_number, 
                                                          hpv_distribution_ratio])]
      # prevalence
      burden [, prev.cecx := (prev.cecx * psadat_country [run_id == run_number, 
                                                          prevalence_ratio]
                                        * psadat_country [run_id == run_number, 
                                                          hpv_distribution_ratio])]
      
      # YLL -- premature mortality
      burden [, lifey := (lifey * psadat_country [run_id == run_number, 
                                                  mortality_ratio]
                                * psadat_country [run_id == run_number, 
                                                  hpv_distribution_ratio])]
      
      # ------------------------------------------------------------------------
      # estimation of YLD -- disability
      
      # disability weights for different phases of cervical cancer
      # (diagnosis & therapy, controlled, metastatic, terminal)
      dw = list (diag       = psadat_country [run_id == run_number, 
                                              dw_diagnosis], 
                 control    = psadat_country [run_id == run_number, 
                                              dw_control], 
                 metastatic = psadat_country [run_id == run_number, 
                                              dw_metastatic], 
                 terminal   = psadat_country [run_id == run_number, 
                                              dw_terminal] 
                 )
      
      # duration of different phases of cervical cancer
      # (diagnosis & therapy, controlled, metastatic, terminal) -- unit in years
      disability.weights <- "gbd_2017"
      
      cecx_duration = list (diag       = data.disability_weights [Source == disability.weights &
                                                                    Sequela=="diagnosis",
                                                                  Duration],
                            metastatic = data.disability_weights [Source == disability.weights &
                                                                    Sequela=="metastatic",
                                                                  Duration],
                            terminal   = data.disability_weights [Source == disability.weights &
                                                                    Sequela=="terminal",
                                                                  Duration])
      # duration of controlled phases is based on remainder of 
      # time after attributing to other phases
      
      # YLD -- disability
      # combine yld contribution from (incidence, prevalence and mortality) cases
      burden [, disability := (inc.cecx  * dw$diag * cecx_duration$diag) + 
                              (prev.cecx * dw$control) +
                              (mort.cecx * ( (dw$metastatic * cecx_duration$metastatic) +
                                             (dw$terminal   * cecx_duration$terminal) ) ) ]

      # ------------------------------------------------------------------------
      
      # ------------------------------------------------------------------------
      
      # convert results to vimc format
      convert_results <- OutputVimc (DT            = burden, 
                                     calendar_year = TRUE, 
                                     vimc_template = vimc_template_country)
      
      # drop scenario (post-vaccination or pre-vaccination) column
      convert_results [, scenario := NULL]
      
      # add column for run id
      convert_results [, run_id := run_number]
      
      # add stochastic burden estimate adjusted by sensitivity analysis to the
      # stochastic burden table of all runs
      stochastic_burden_country <- rbind (stochastic_burden_country, 
                                          convert_results,
                                          use.names = TRUE)
      
    } # end of -- for (run_id in 1:psa_runs)
    
    # save stochastic burden estimates of current country
    # fwrite (x      = stochastic_burden_country, 
    #         file   = stochastic_file_country, 
    #         append = FALSE)
    
    # save stochastic burden estimates of current country to larger file with 
    # stochastic burden estimates of all countries
    fwrite (x      = stochastic_burden_country, 
            file   = stochastic_file, 
            append = TRUE)
    
  } # end of -- for (country_code in countryCodes)

  return ()
  
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
print (Sys.time ())  # current time

# save current folder name and move to base folder
current_folder <- getwd ()
# current_folder <- setwd ("D:/GitHub/vimc_hpv/code")  # DEBUG -- remove this line later
# current_folder <- setwd ("F:/201910gavi/hpv/code")   # DEBUG -- remove this line later
setwd ("../")
#-------------------------------------------------------------------------------

# montangu touchstone
touchstone <- "201910gavi-4"

# central burden template file
central_burden_template_file <- paste0 ("input/central-burden-template.", 
                                        touchstone, 
                                        ".HPV_LSHTM-Jit_standard.csv")

#-------------------------------------------------------------------------------
# initialisation for probabilistic sensitivity analysis
run_central   <- TRUE  # logical, run/not run central analysis
run_lhs       <- TRUE  # generate latin hyper sample of input parameters for psa
run_psa       <- TRUE  # logical, run/not run PSA for vaccination scenarios
run_psa_novac <- FALSE  # logical, run/not run PSA for no vaccination scenario
psa_runs      <- 200   # number of runs for psa
seed_state    <- 1
vaccine       <- "4vHPV"

# files to save input parameter distributions foFALSEr probabilistic sensitivity analysis
psadat_file      <- "output_psa/psadat.csv"      
psadat_vimc_file <- paste0 ("output_psa/stochastic_parameters_vimc_", 
                            touchstone, ".csv")  # vimc format


# generate input parameter distributions for probabilistic sensitivity analysis
if (run_lhs) {
  
  # Read in file with a column "country" containing iso3 country codes
  country_table <- fread (central_burden_template_file)
  
  # get unique list of country codes (iso3)
  country_codes <- unique (country_table [, country])
  
  # create psa data for probabilistic sensitivity analysis
  psadat_list <- CreatePsaData (country_codes    = country_codes,
                                vaccine          = vaccine,
                                psa_runs         = psa_runs,
                                seed_state       = seed_state, 
                                psadat_file      = psadat_file,
                                psadat_vimc_file = psadat_vimc_file)
}
#-------------------------------------------------------------------------------

# scenarios
scenarios <- c ("hpv-routine-default",
                "hpv-routine-bestcase",
                "hpv-campaign-default",
                "hpv-campaign-bestcase")
scenarios <- c ("hpv-routine-default")  # DEBUG -- comment this line later

# loop through different scenarios
for (scenario in scenarios) {
  
  #-----------------------------------------------------------------------------
  # Generate vaccine impact estimates (central run)
  tic ()  # start timer
  
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
  
  print (paste0 ("burden template: ", central_burden_template_file))
  print (paste0 ("scenario: ", scenario))
  print (paste0 ("coverage: ", vaccine_coverage_file))
  
  # Generate vaccine impact estimates (central run)
  # for vaccination scenario
  if (run_central) {
    
    EstimateVaccineImpactVimcCentral (
      vaccine_coverage_file              = vaccine_coverage_file,
      disease_burden_template_file       = central_burden_template_file,
      disease_burden_no_vaccination_file = central_burden_no_vaccination_file,
      disease_burden_vaccination_file    = central_burden_vaccination_file,
      disease_burden_results_file        = central_burden_results_file,
      routine_vaccination                = routine,
      campaign_vaccination               = campaign, 
      vaccine                            = vaccine
    )
  }
  
  toc ()  # note current timer and compute elapsed time
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Generate vaccine impact estimates (stochastic run - vaccination scenario)
  tic ()  # start timer
  
  # probabilistic sensitivity analysis for vaccination scenarios
  if (run_psa) {
    
    # directory for stochastic burden files
    stochastic_burden_dir <-  paste0 ("output_psa/",
                                      scenario, "/")

    
    # filename for stochastic burden estimates (without .csv at the end)
    stochastic_burden_vaccination_file <- paste0 ("stochastic-burden-vaccination_", 
                                                  touchstone, "_", 
                                                  scenario)
    
    # Generate vaccine impact estimates (sensitivity analysis/stochastic runs)
    # for vaccination scenario
    EstimateVaccineImpactVimcStochastic (
      disease_burden_template_file  = central_burden_template_file,
      centralBurdenResultsFile      = central_burden_results_file,
      psaData                       = psadat_list [["psadat"]],
      diseaseBurdenStochasticFolder = stochastic_burden_dir,
      diseaseBurdenStochasticFile   = stochastic_burden_vaccination_file, 
      psa_runs                      = psa_runs, 
      countryCodes                  = -1,
      vaccination_scenario          = TRUE
      ) 
  }
  
  toc ()
  #-----------------------------------------------------------------------------
  
} # end of loop -- for (scenario in scenarios)


#-----------------------------------------------------------------------------
# Generate vaccine impact estimates (stochastic run - vaccination scenario)
tic ()

# probabilistic sensitivity analysis for no vaccination scenario
if (run_psa_novac) {

  # directory for stochastic burden files
  stochastic_burden_dir <-  paste0 ("output_psa/",
                                    "hpv-no-vaccination", "/")
  
  # filename for stochastic burden estimates (without .csv at the end)
  stochastic_burden_novaccination_file <- paste0 ("stochastic-burden-novaccination_", 
                                                  touchstone, "_", 
                                                  "hpv-no-vaccination")
  
  # Generate vaccine impact estimates (sensitivity analysis/stochastic runs)
  # for no vaccination scenario
  EstimateVaccineImpactVimcStochastic (
    disease_burden_template_file  = central_burden_template_file,
    centralBurdenResultsFile      = central_burden_results_file,
    psaData                       = psadat_list [["psadat"]],
    diseaseBurdenStochasticFolder = stochastic_burden_dir,
    diseaseBurdenStochasticFile   = stochastic_burden_novaccination_file, 
    psa_runs                      = psa_runs, 
    countryCodes                  = -1,
    vaccination_scenario          = FALSE
  ) 
}

toc ()
#-----------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# return to source file location/folder
setwd (current_folder)

print (Sys.time ())  # current time
#-------------------------------------------------------------------------------
# main program - end 
#-------------------------------------------------------------------------------







