## effect of number of species

# vector with possible ASVs
available_ASVs <- all_years_rarefaction %>% 
  ungroup() %>% 
  select(Sequence) %>% 
  distinct() %>%
  pull(Sequence)

#
all_years_rarefaction %>% 
  filter(Sequence %in% sample(available_ASVs, 100)) %>% 
  define_rb()
