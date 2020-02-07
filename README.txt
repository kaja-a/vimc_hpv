VIMC HPV runs (Jan 2020)

Folders:
- input - files downloaded from montagu

            - central-burden-template.201910gavi-4.HPV_LSHTM-Jit_standard.csv
            - stochastic-burden-template.201910gavi-4.HPV_LSHTM-Jit_standard.csv

            - stochastic_template_params.csv

            - coverage_201910gavi-4_hpv-no-vaccination.csv

            - coverage_201910gavi-4_hpv-routine-default.csv
            - coverage_201910gavi-4_hpv-campaign-default.csv

            - coverage_201910gavi-4_hpv-routine-bestcase.csv  
            - coverage_201910gavi-4_hpv-campaign-bestcase.csv

              (routine refers to routine + campaign)

- output - files to upload to montagu (central runs)

            - psadat_vimc.csv   (stochastic parameters values)

- output_psa - files to upload to dropbox (VIMC)

- log - log files to keep track of simulation runs


----------------------------------------------------------------------
Scenarios (5): 

- Best case - Campaign              (ID: hpv-campaign-bestcase)

- Best case - Routine and Campaign  (ID: hpv-routine-bestcase)

- Default - Campaign                (ID: hpv-campaign-default)

- Default - Routine and Campaign    (ID: hpv-routine-default)

- No vaccination                    (ID: hpv-no-vaccination)

----------------------------------------------------------------------
Download:

- Vaccine coverage data (5)

- Central disease burden template (1)
    | disease |          year | age | country | country_name | cohort_size | cases | dalys | deaths |

- Stochastic disease burden template (1)
    | disease | run_id | year | age | country | country_name | cohort_size | cases | dalys | deaths |

- Stochastic parameters template (1)
    | run_id | <param_1> | <param_2> |

----------------------------------------------------------------------
Upload: 

- Central burden estimates (5)

- Stochastic burden estimates (5)

- Stochastic parameters (1)

----------------------------------------------------------------------

HPV:LSHTM-Jit:standard

For scenarios:
Best case - Campaign (hpv-campaign-bestcase)
Best case - Routine and Campaign (hpv-routine-bestcase)
Default - Campaign (hpv-campaign-default)
Default - Routine and Campaign (hpv-routine-default)
No vaccination (hpv-no-vaccination)

Expecting data on:
cohort_size
cases
dalys
deaths

For all combinations of:
112 countries: view list
101 years: 2000 - 2100
92 ages: 9 - 100
Not including cohorts born before 1900 or after 2030