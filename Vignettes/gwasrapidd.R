library(tidyverse)
library(gwasrapidd)

get_variants(variant_id = 'rs1050501') -> vobj

my_studies <- get_studies(study_id = 'GCST000858', variant_id = 'rs12752552')

get_studies(variant_id = 'rs12752552') -> vobj
