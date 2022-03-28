#### R script to process PlateReader data ####
# by Charlotte Kunze
# 02.03.22


# packages
library(tidyverse)
library(readxl)
library(here)

#----------------------------------------------------------------------------------------------------------------#
#Einlesen des plate-Files, welches Informationen darueber enthaelt welches well im plate mit welchem Planktotron belegt ist/welches treatment wo ist
mono  <- read_excel("~/Desktop/Exp22/Pipettierschema.xlsx", 
                                      sheet = "MONO") %>%
  mutate(combination = paste('mono'))
duo1  <- read_excel("~/Desktop/Exp22/Pipettierschema.xlsx", 
                    sheet = "DUO1")%>%
  mutate(combination = paste('duo1'))
duo2  <- read_excel("~/Desktop/Exp22/Pipettierschema.xlsx", 
                    sheet = "DUO2")%>%
  mutate(combination = paste('duo2'))
quattro  <- read_excel("~/Desktop/Exp22/Pipettierschema.xlsx", 
                    sheet = "QUATTRO")%>%
  mutate(combination = paste('quattro'))

plastestr <- rbind(mono, duo1, duo2, quattro)
#----------------------------------------------------------------------------------------------------------------#
setwd("~/Desktop/Exp22/Data/txtData")

#Arbeitsverzeichnis setzen (die folgende Zeile setzt das Arbeitsverzeichnis an den Ort des R-Skriptes)
directory <- paste(getwd(),sep = '') # hier sind meine Daten

#read file (from [here][1])
file.list<- list.files(full.names = TRUE, pattern = "*.txt")


#Rausziehen der Dateipfade:
files <- dir(directory, recursive=TRUE, full.names=TRUE, pattern="\\.txt$")

#Erstellen eines leeren "Blanko-Objektes", welches spaeter alle Daten speichern soll
expData <- NULL

#----------------------------------------------------------------------------------------------------------------#

for (i in 1:length(files)){
  data <- read_table(files[i])  %>%
    slice(-c(1:37)) %>%
    slice(-c(7:9, 16:18, 25:27,34:36,43:45))%>%
    rename(dummy = 'Software Version	3.03.14') %>%
    separate(dummy, into = c('row', 'V1','V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'method'), sep = '\t') %>%
    mutate(fileName = file.list[i] ) %>%
    mutate(dummy = str_remove(fileName, '.txt')) %>%
    mutate(dummy = str_remove(dummy, './')) %>%
    mutate(V1 = as.numeric(gsub(',', '.', V1)),
           V2 = as.numeric(gsub(',', '.', V2)),
           V3 = as.numeric(gsub(',', '.', V3)),
           V4 = as.numeric(gsub(',', '.', V4)),
           V5 = as.numeric(gsub(',', '.', V5)),
           V6 = as.numeric(gsub(',', '.', V6)),
           V7 = as.numeric(gsub(',', '.', V7)),
           V8 = as.numeric(gsub(',', '.', V8))) %>%
    separate(dummy, into = c('treat', 'combination','sampling'), '_') %>%
    gather(key = column, value = value, c(2:9)) %>%
    mutate(sampling = str_remove(sampling,'s'),
           sampling = as.numeric(sampling)) %>%
    mutate(method = str_remove(method,'Read 1:'),
          method = str_remove(method,'Read 2:'))
    
  expData <- rbind(data, expData)
}


#adjust treatment levels
expData$treat[expData$treat == 'Inc']<-'inc'
expData$treat[expData$treat == 'Fluc']<-'fluc'
expData$treat[expData$treat == 'CS']<-'cs'

unique(expData$treat)



#zusammenfuegen mit platten info
expData_noMIX <- left_join(expData, plastestr, by = c('column', 'row', 'combination'))
is.na(expData_noMIX$ID)
#write.csv2(expData_noMIX, file = 'expData_noMIX.csv')


#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
setwd("~/Desktop/Exp22/Data/Kalibrierung_abt6")

#Arbeitsverzeichnis setzen (die folgende Zeile setzt das Arbeitsverzeichnis an den Ort des R-Skriptes)
directory <- paste(getwd(),sep = '') # hier sind meine Daten

#read file (from [here][1])
file.list<- list.files(full.names = TRUE, pattern = "*.txt")


#Rausziehen der Dateipfade:
files <- dir(directory, recursive=TRUE, full.names=TRUE, pattern="\\.txt$")

#Erstellen eines leeren "Blanko-Objektes", welches spaeter alle Daten speichern soll
blankData <- NULL

#----------------------------------------------------------------------------------------------------------------#

for (i in 1:length(files)){
  blank <- read_table(files[i])  %>%
    slice(-c(1:37)) %>%
    slice(-c(7:9, 16:18, 25:27,34:36,43:45))%>%
    rename(dummy = 'Software Version	3.03.14') %>%
    separate(dummy, into = c('row', 'V1','V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'method'), sep = '\t') %>%
    mutate(fileName = file.list[i] ) %>%
    mutate(dummy = str_remove(fileName, '.txt')) %>%
    mutate(dummy = str_remove(dummy, './')) %>%
    mutate(V1 = as.numeric(gsub(',', '.', V1)),
           V2 = as.numeric(gsub(',', '.', V2)),
           V3 = as.numeric(gsub(',', '.', V3)),
           V4 = as.numeric(gsub(',', '.', V4)),
           V5 = as.numeric(gsub(',', '.', V5)),
           V6 = as.numeric(gsub(',', '.', V6)),
           V7 = as.numeric(gsub(',', '.', V7)),
           V8 = as.numeric(gsub(',', '.', V8))) %>%
    separate(dummy, into = c('treat', 'combination','X'), '_') %>%
    gather(key = column, value = value, c(2:9)) %>%
        mutate(method = str_remove(method,'Read 1:'),
           method = str_remove(method,'Read 2:')) %>%
    rename(blank = value) %>%
    select(-fileName, -X)
  
  blankData <- rbind(blank, blankData)
 
  }
unique(blankData$treat)

#adjust treatment levels
blankData$treat[blankData$treat == 'incfluc']<-'flucinc'
unique(blankData$treat)



#zusammenfuegen mit platten info
names(blankData)
expData_noMix_mitBlank <- left_join(blankData, expData_noMIX, by = c('column', 'row', 'combination', 'method', 'treat'))
is.na(expData_noMIX$ID)
write.csv2(expData_noMix_mitBlank, file = 'expData_noMix_mitBlank.csv')




#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#### MIX DATA ####

mix  <- read_excel("~/Desktop/Exp22/Pipettierschema.xlsx", 
                   sheet = "MIX")%>%
  mutate(combination = paste('mix'))


#----------------------------------------------------------------------------------------------------------------#
setwd("~/Desktop/Exp22/Data/mixData")

#Arbeitsverzeichnis setzen (die folgende Zeile setzt das Arbeitsverzeichnis an den Ort des R-Skriptes)
directory <- paste(getwd(),sep = '') # hier sind meine Daten

#read file (from [here][1])
file.list<- list.files(full.names = TRUE, pattern = "*.txt")


#Rausziehen der Dateipfade:
files <- dir(directory, recursive=TRUE, full.names=TRUE, pattern="\\.txt$")

#Erstellen eines leeren "Blanko-Objektes", welches spaeter alle Daten speichern soll
mixData <- NULL

#----------------------------------------------------------------------------------------------------------------#

for (i in 1:length(files)){
  data1 <- read_table(files[i])  %>%
    slice(-c(1:37)) %>%
    slice(-c(7:9, 16:18, 25:27,34:36,43:45))%>%
    rename(dummy = 'Software Version	3.03.14') %>%
    separate(dummy, into = c('row', 'V1','V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'method'), sep = '\t') %>%
    mutate(fileName = file.list[i] ) %>%
    mutate(dummy = str_remove(fileName, '.txt')) %>%
    mutate(dummy = str_remove(dummy, './')) %>%
    mutate(V1 = as.numeric(gsub(',', '.', V1)),
           V2 = as.numeric(gsub(',', '.', V2)),
           V3 = as.numeric(gsub(',', '.', V3)),
           V4 = as.numeric(gsub(',', '.', V4)),
           V5 = as.numeric(gsub(',', '.', V5)),
           V6 = as.numeric(gsub(',', '.', V6)),
           V7 = as.numeric(gsub(',', '.', V7)),
           V8 = as.numeric(gsub(',', '.', V8))) %>%
    separate(dummy, into = c('combination','sampling'), '_') %>%
    gather(key = column, value = value, c(2:9)) %>%
    mutate(sampling = str_remove(sampling,'s'),
           sampling = as.numeric(sampling)) %>%
    mutate(method = str_remove(method,'Read 1:'),
           method = str_remove(method,'Read 2:'))
  
  mixData <- rbind(data1, mixData)
}

unique(mixData$combination)



#zusammenfuegen mit platten info
expData_MIX <- left_join(mixData, mix, by = c('column', 'row', 'combination'))
is.na(expData_MIX$ID)


MIX_blank <- read_table('~/Desktop/Exp22/Data/mixBlank/MIX_blank.txt') %>%
  slice(-c(1:37)) %>%
  slice(-c(7:9, 16:18, 25:27,34:36,43:45))%>%
  rename(dummy = 'Software Version	3.03.14') %>%
  separate(dummy, into = c('row', 'V1','V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'method'), sep = '\t') %>%
  mutate(V1 = as.numeric(gsub(',', '.', V1)),
         V2 = as.numeric(gsub(',', '.', V2)),
         V3 = as.numeric(gsub(',', '.', V3)),
         V4 = as.numeric(gsub(',', '.', V4)),
         V5 = as.numeric(gsub(',', '.', V5)),
         V6 = as.numeric(gsub(',', '.', V6)),
         V7 = as.numeric(gsub(',', '.', V7)),
         V8 = as.numeric(gsub(',', '.', V8))) %>%
  gather(key = column, value = value, c(2:9)) %>%
    mutate(method = str_remove(method,'Read 1:'),
         method = str_remove(method,'Read 2:')) %>%
  rename(blank = value)

MIX_mitBlank <- left_join(expData_MIX, MIX_blank, by = c('method', 'row', 'column'))
#write.csv2(expData_MIX, file = 'MIX_mitBlank.csv')

#daten zusammenfuehren
All_Exp22 <- bind_rows(expData_noMix_mitBlank, MIX_mitBlank)%>%
  filter(ID != 'MQ')
unique(All_Exp22$combination)
write.csv(All_Exp22, file = here('~/Desktop/Exp22/output','allData_wBlankExp22.csv'))


