# --------------------------------------------------------------------------------- get data from GitHub
wd <- "https://github.com/Datactuariat/OpenSinistres/blob/master/Data/"

# --------------------------------------------------------------------------------- required packages
library(stringr)
library(data.table)
library(readr)
library(dplyr)
library(zoo)

# --------------------------------------------------------------------------------- Loading RAW data

# La base Etalab de données des accidents corporels de la circulation d'une année donnée, est répartie
# en 4 rubriques sous la forme pour chacune d'elles d'un fichier au format csv.
#    1. La rubrique CARACTERISTIQUES qui décrit les circonstances générales de l’accident
#    2. La rubrique LIEUX qui décrit le lieu principal de l’accident
#    3. La rubrique VEHICULES impliqués
#    4. La rubrique USAGERS impliqués

Caracteristiques <- list(
  C2009       = read_delim(paste0(wd,"caracteristiques_2009.csv?raw=true"),delim="\t",col_types="ciiiiiiiiiicc??i"),
  C2010       = read_csv(paste0(wd,"caracteristiques_2010.csv?raw=true"),col_types="ciiiiiiiiiicc??i"),
  C2011       = read_csv(paste0(wd,"caracteristiques_2011.csv?raw=true"),col_types="ciiiiiiiiiicc??i"),
  C2012       = read_csv(paste0(wd,"caracteristiques_2012.csv?raw=true"),col_types="ciiiiiiiiiicc??i"),
  C2013       = read_csv(paste0(wd,"caracteristiques_2013.csv?raw=true"),col_types="ciiiiiiiiiicc??i"),
  C2014       = read_csv(paste0(wd,"caracteristiques_2014.csv?raw=true"),col_types="ciiiiiiiiiicc??i")
)

Usagers <- list(
  U2009       = read_csv(paste0(wd,"usagers_2009.csv?raw=true"),col_types="ciiiiiiiiiic"),
  U2010       = read_csv(paste0(wd,"usagers_2010.csv?raw=true"),col_types="ciiiiiiiiiic"),
  U2011       = read_csv(paste0(wd,"usagers_2011.csv?raw=true"),col_types="ciiiiiiiiiic"),
  U2012       = read_csv(paste0(wd,"usagers_2012.csv?raw=true"),col_types="ciiiiiiiiiic"),
  U2013       = read_csv(paste0(wd,"usagers_2013.csv?raw=true"),col_types="ciiiiiiiiiic"),
  U2014       = read_csv(paste0(wd,"usagers_2014.csv?raw=true"),col_types="ciiiiiiiiiic")
)

Vehicules <- list(
  V2009 = read_csv(paste0(wd,"vehicules_2009.csv?raw=true"),col_types="ciiiiiiic"),
  V2010 = read_csv(paste0(wd,"vehicules_2010.csv?raw=true"),col_types="ciiiiiiic"),
  V2011 = read_csv(paste0(wd,"vehicules_2011.csv?raw=true"),col_types="ciiiiiiic"),
  V2012 = read_csv(paste0(wd,"vehicules_2012.csv?raw=true"),col_types="ciiiiiiic"),
  V2013 = read_csv(paste0(wd,"vehicules_2013.csv?raw=true"),col_types="ciiiiiiic"),
  V2014 = read_csv(paste0(wd,"vehicules_2014.csv?raw=true"),col_types="ciiiiiiic")
)

Lieux <- list(
  L2009 = read_csv(paste0(wd,"lieux_2009.csv?raw=true"),col_types="c?????????????????"),
  L2010 = read_csv(paste0(wd,"lieux_2010.csv?raw=true"),col_types="c?????????????????"),
  L2011 = read_csv(paste0(wd,"lieux_2011.csv?raw=true"),col_types="c?????????????????"),
  L2012 = read_csv(paste0(wd,"lieux_2012.csv?raw=true"),col_types="c?????????????????"),
  L2013 = read_csv(paste0(wd,"lieux_2013.csv?raw=true"),col_types="c?????????????????"),
  L2014 = read_csv(paste0(wd,"lieux_2014.csv?raw=true"),col_types="c?????????????????")
)

# --------------------------------------------------------------------------------- rbind

Caracteristiques <- rbindlist(Caracteristiques,fill=TRUE)
Caracteristiques[,date:=as.Date(paste(2000+an,mois,jour,sep="-"),"%Y-%m-%d")]
Usagers <- rbindlist(Usagers,fill=TRUE)
Vehicules <- rbindlist(Vehicules)
Lieux <- rbindlist(Lieux,fill=TRUE)

# -------------------
cat("test")
