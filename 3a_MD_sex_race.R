library(readxl);library(tidyverse);library(stringr)

# AMR-AFR
#Read in sex and race information

# Clinical MOD records  #from Rebecca
MOD_clinical <- read_excel("MOD_subject_sex_race_info_7.7.23.xlsx") %>%   
  mutate(sex=ifelse(Sex == 'Female', 2, 1)) 

# PRISM info from Rebecca
PRISM <- read.csv("PRISM to MOD-PRISM id mapping_from MS 3.17.18 w MS sasmple ID info added.csv") %>% 
  rename(PRISM.ID=1,
         mom=3,
         baby=4,
         Sex=6) 

# GOH info from Emma
GOH1 <- read.csv("MOD_OBER_Plate7_et.csv")[1:21, ] %>% 
  mutate(ID=paste0(Participant.ID, "_p7")) %>% 
  mutate(sex=ifelse(Sex == 'Female', 2, 1)) %>% 
  mutate(race=ifelse(str_starts(Race, "Black"), 'Black', 
                     ifelse(str_starts(Race, "White"), 'White',
                            ifelse(Race == "AA", "Black", Race)))) %>% 
  select(7:9)

# GOH info from NATHAN
GOH2 <- read.table("ROBIs_Genotyped_052023_update.csv", header = TRUE) %>% 
  mutate(ID=paste0(ID, "_p14")) %>% 
  mutate(sex=ifelse(Sex == 'Female', 2, 1)) %>% 
  rename(race=3) %>% 
  select(1, 4, 3)

# NUAA info from Chris
NUAA <- read_csv("Inventory of Q576R Genotyped NU DNA to UC 8-9-21.csv", show_col_types = FALSE) %>%
  rename(uc_id=1) %>% 
  mutate(ID=paste0(uc_id, '_p9')) %>%
  mutate(sex=ifelse(Sex == 'Female', 2, 1)) %>% 
  mutate(race=ifelse(str_starts(Race, "Black"), 'Black', Race)) %>% 
  select(8:10)


#Add sex and race to samples with 5 digits ID

# Samples in AMR-AFR ----- 769
mod2123ID <- read.table("mod2123.list")   #

# Add sex to samples with *****MOD --5_digit ID
mod5dig <- mod2123ID %>% 
  rename(ID=V1) %>% 
  filter(str_detect(ID, '^[1-9]')) %>% 
  filter(nchar(ID) > 7) %>% 
  mutate(id = str_sub(ID, 1, 7)) %>% 
  left_join(MOD_clinical, by=c('id'='MOD_ID')) %>% 
  rename(race=Race) 

# Add sex and race info for 13807 appeared in the Alternative_ID column
mod5dig$sex[mod5dig$id == '13807-I'] <- 1    # 13807-I male
mod5dig$sex[mod5dig$id == '13807-M'] <- 2    # 13807-M female
mod5dig$race[str_starts(mod5dig$id, '13807-')] <- 'Black'    # 13807-M female

mod5d <- mod5dig %>% 
  select(1, 9, 7)

# Add sex and race to PRISM samples

mom <- subset(PRISM, mom!="") %>% 
  mutate(ID=paste0(mom, "-M")) %>% 
  mutate(sex=2,
         id=mom) %>% 
  mutate(race=ifelse(Race=='African-American', 'Black', Race)) %>% 
  select(9,7,8,10)


baby <- subset(PRISM, baby!="") %>% 
  mutate(ID=paste0(baby, "-I")) %>% 
  mutate(sex=ifelse(Sex == 'F', 2, 1),
         id=baby) %>% 
  mutate(race=ifelse(Race=='African-American', 'Black', Race)) %>% 
  select(9, 7, 8, 10)

write.csv(mom, 'EFZmom_sex_race_AFR.csv', row.names = FALSE, quote = FALSE)
write.csv(baby, 'EFZbaby_sex_race_AFR.csv', row.names = FALSE, quote = FALSE)

EFZ <- mod2123ID %>% 
  filter(str_starts(V1, 'EFZ')) %>% 
  mutate(ID=str_sub(V1, 1, 10)) %>% 
  left_join(rbind(mom, baby), by=c('ID'='id')) %>% 
  select(1, 4, 5) %>% 
  rename(ID=V1)

# Add sex and race to GOH

GOH <- mod2123ID %>% 
  inner_join(rbind(GOH1, GOH2),by=c('V1'='ID')) %>% 
  rename(ID=V1)


# Add sex and race to NAPs and dp/db samples   

#Checked all cell line and dTL samples female and Black like NAPs
#https://mnlab.uchicago.edu/mod/genotyping/E1-E10/full.html

mod_nuaa_prism_goh <- rbind(mod5d, NUAA, EFZ, GOH)

MOD_sex_race <- mod2123ID %>% 
  left_join(mod_nuaa_prism_goh, by=c('V1'='ID')) %>%
  mutate(Sex=ifelse(is.na(sex), 2, sex)) %>% 
  mutate(Race=ifelse(is.na(race), 'Black', race)) %>% 
  rename(ID=V1) %>% 
  select(1, 4,5)

write.csv(MOD_sex_race, 'MOD_sex_race_AFR.csv', row.names = FALSE, quote = FALSE)

