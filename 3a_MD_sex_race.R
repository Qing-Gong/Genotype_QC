library(readxl);library(tidyverse);library(stringr)

# AMR-AFR
#Read in sex and race information

# Clinical MD records 
MD_clinical <- read_excel("sex_race_info_23.xlsx") %>%   
  mutate(sex=ifelse(Sex == 'Female', 2, 1)) 

# PM info 
PM <- read.csv("PM.csv") %>% 
  rename(PM.ID=1,
         mom=3,
         baby=4,
         Sex=6) 

# GH info
GH1 <- read.csv("Plt7.csv")[1:21, ] %>% 
  mutate(ID=paste0(Participant.ID, "_p7")) %>% 
  mutate(sex=ifelse(Sex == 'Female', 2, 1)) %>% 
  mutate(race=ifelse(str_starts(Race, "Black"), 'Black', 
                     ifelse(str_starts(Race, "White"), 'White',
                            ifelse(Race == "AA", "Black", Race)))) %>% 
  select(7:9)

# GH info
GH2 <- read.table("R_update.csv", header = TRUE) %>% 
  mutate(ID=paste0(ID, "_p14")) %>% 
  mutate(sex=ifelse(Sex == 'Female', 2, 1)) %>% 
  rename(race=3) %>% 
  select(1, 4, 3)

# NAA info
NAA <- read_csv("INU.csv", show_col_types = FALSE) %>%
  rename(uc_id=1) %>% 
  mutate(ID=paste0(uc_id, '_p9')) %>%
  mutate(sex=ifelse(Sex == 'Female', 2, 1)) %>% 
  mutate(race=ifelse(str_starts(Race, "Black"), 'Black', Race)) %>% 
  select(8:10)


#Add sex and race to samples with 5 digits ID

# Samples in AMR-AFR ----- 769
MD2123ID <- read.table("MD2123.list")   #

# Add sex to samples with *****MD --5_digit ID
md5dig <- MD2123ID %>% 
  rename(ID=V1) %>% 
  filter(str_detect(ID, '^[1-9]')) %>% 
  filter(nchar(ID) > 7) %>% 
  mutate(id = str_sub(ID, 1, 7)) %>% 
  left_join(MD_clinical, by=c('id'='MD_ID')) %>% 
  rename(race=Race) 

# Add sex and race info for XXXXX appeared in the Alternative_ID column
md5dig$sex[md5dig$id == 'XXXXX-I'] <- 1    # XXXXX-I male
md5dig$sex[md5dig$id == 'XXXXX-M'] <- 2    # XXXXX-M female
md5dig$race[str_starts(md5dig$id, 'XXXXX-')] <- 'Black'    # XXXXX-M female

md5d <- md5dig %>% 
  select(1, 9, 7)

# Add sex and race to PM samples

mom <- subset(PM, mom!="") %>% 
  mutate(ID=paste0(mom, "-M")) %>% 
  mutate(sex=2,
         id=mom) %>% 
  mutate(race=ifelse(Race=='African-American', 'Black', Race)) %>% 
  select(9,7,8,10)


baby <- subset(PM, baby!="") %>% 
  mutate(ID=paste0(baby, "-I")) %>% 
  mutate(sex=ifelse(Sex == 'F', 2, 1),
         id=baby) %>% 
  mutate(race=ifelse(Race=='African-American', 'Black', Race)) %>% 
  select(9, 7, 8, 10)

write.csv(mom, 'EZmom_sex_race_AFR.csv', row.names = FALSE, quote = FALSE)
write.csv(baby, 'EZbaby_sex_race_AFR.csv', row.names = FALSE, quote = FALSE)

EZ <- MD2123ID %>% 
  filter(str_starts(V1, 'EZ')) %>% 
  mutate(ID=str_sub(V1, 1, 10)) %>% 
  left_join(rbind(mom, baby), by=c('ID'='id')) %>% 
  select(1, 4, 5) %>% 
  rename(ID=V1)

# Add sex and race to GOH

GH <- MD2123ID %>% 
  inner_join(rbind(GH1, GH2),by=c('V1'='ID')) %>% 
  rename(ID=V1)


# Add sex and race to dp/db samples   

#Checked all cell line and XX samples female and Black

mnpg <- rbind(md5d, NAA, EZ, GH)

MD_sex_race <- MD2123ID %>% 
  left_join(mnpg, by=c('V1'='ID')) %>%
  mutate(Sex=ifelse(is.na(sex), 2, sex)) %>% 
  mutate(Race=ifelse(is.na(race), 'Black', race)) %>% 
  rename(ID=V1) %>% 
  select(1, 4,5)

write.csv(MD_sex_race, 'MD_sex_race_AFR.csv', row.names = FALSE, quote = FALSE)

