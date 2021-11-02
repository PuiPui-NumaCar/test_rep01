###############################################################################
# Basic setup
###############################################################################

# Source setup file
source(file.path("/xxx/_setup.R"))
source(file.path("/xxx/funcs.R"))

# Set working directory
study="/xxx/"
setwd(study)

# Load required libraries
pacman::p_load(haven, data.table, dplyr, stringr, striprtf, arsenal)

###############################################################################
# Read in required datasets
###############################################################################

library(foreign)
adsl = data.table(read.xport(file.path(setup_list$analysis, "adsl.xpt")))
adae = data.table(read.xport(file.path(setup_list$analysis, "adae_j.xpt"))
)[SAFFL=="Y" & !is.na(AEDECOD) & TRTEMFL=="Y"]

###############################################################################
# Derivations
###############################################################################

adae2 <- merge(adsl, adae,)

# Required variables
vars = c("USUBJID","AESOC","AEDECOD","AESOCJ","AEDECODJ","TRTEMFL",
         "TRT01AN","TRT01A")

grd = copy(adae2[,..vars])

# Function to count per arm and period and transpose the counts
cnt = function(trt,col)
{
  # For all columns except total patients
  if(col!=2)
  {
    grd = grd[TRT01AN %in% trt]
    
    grd1 = unique(grd[,.(USUBJID,TRT01AN,TRT01A)])
    cnt1 = grd1[, .N,by = c("TRT01AN", "TRT01A")][
      ,`:=` (AEDECOD = "Number of subjects with at least one event",
             AESOC = "AAA",
             ORD=1)]
    
    grd2 = unique(grd[,.(USUBJID,TRT01AN,TRT01A,AESOC)])
    cnt2 = grd2[, .N,by = c("TRT01AN", "TRT01A", "AESOC")][,ORD:=2]
    
    grd3 = unique(grd[,.(USUBJID,TRT01AN,TRT01A,AESOC,AEDECOD)])
    cnt3 = grd3[, .N,by = c("TRT01AN", "TRT01A", "AESOC","AEDECOD")][,ORD:=2]
    
    # If at least 1 AE in any grade
    if(nrow(cnt1)>0 | nrow(cnt2)>0) 
    {
      combine = data.table(rbind(cnt1,cnt2,cnt3,fill=TRUE))
      combine = combine[,`:=` (SUBORD = ifelse(is.na(AEDECOD),1,2),
                               AEDECOD = ifelse(is.na(AEDECOD),AESOC,AEDECOD))][order(AESOC)]
      trans_dset = dcast(combine,
                         ORD+SUBORD+AESOC+AEDECOD ~ TRT01AN,
                         value.var = "N")
      colnames(trans_dset) = paste(as.character(col), colnames(trans_dset), sep="_")
    } else                #If 0 AE in all grades
    {
      cols = paste(col,c("ORD","SUBORD","AESOC","AEDECOD", trt), sep="_")
      trans_dset =  data.table(ord = numeric(), subord = numeric(), aesoc = character(), aedecod = character(), t1 = character())
      setnames(trans_dset, colnames(trans_dset), cols)
    }
    
  } else if(col==2)     #For all patients
  {
    grd = grd[!(TRT01AN==0)]
    grd1 = unique(grd[,.(USUBJID)])
    cnt1 = data.table(N=grd1[,.N])[,`:=` (AEDECOD = "Number of subjects with at least one event",
                                          AESOC = "AAA")][,ORD:=1]
    grd2 = unique(grd[,.(USUBJID,AESOC)])
    cnt2 = grd2[,.N, by = c("AESOC")][,ORD:=2]
    grd3 = unique(grd[,.(USUBJID,AESOC,AEDECOD)])
    cnt3 = grd3[,.N, by = c("AESOC","AEDECOD")][,ORD:=2]
    combine = data.table(rbind(cnt1,cnt2,cnt3,fill=TRUE))
    combine = combine[,`:=` (SUBORD = ifelse(is.na(AEDECOD),1,2),
                             AEDECOD = ifelse(is.na(AEDECOD),AESOC,AEDECOD))][order(AESOC)]
    trans_dset = combine
    colnames(trans_dset) = paste(as.character(col), colnames(trans_dset), sep="_")
  }
  
  return(trans_dset)
}

# Overall study period- Act High arm
o_acth = cnt(trt = 81, col = 1)
# Overall study period- Act low arm
o_actl = cnt(trt = 54, col = 1)
# Overall study period- placebo arm
o_plac = cnt(trt = 0, col = 1)
# Overall study period- all patients count
o_all = cnt(trt = c(54,81,0), col = 2)

View(o_acth)
View(o_actl)
View(o_plac)
View(o_all)

# Add by variables and merge datasets
o_act = as.data.table(o_act)
setkeyv(o_acth, c("1_AESOC", "1_AEDECOD", "1_ORD","1_SUBORD"))
setkeyv(o_actl, c("1_AESOC", "1_AEDECOD", "1_ORD","1_SUBORD"))
setkeyv(o_plac, c("1_AESOC", "1_AEDECOD", "1_ORD","1_SUBORD"))
setkeyv(o_all, c("2_AESOC", "2_AEDECOD", "2_ORD","2_SUBORD"))

final = copy(o_acth[o_actl[o_plac[o_all]]])
setnames(final, c("1_AESOC","1_AEDECOD","1_ORD","1_SUBORD"), c("AESOC","AEDECOD","ORD","SUBORD"))

# Calculate big N values
Ntrt = copy(adsl[SAFFL=="Y", .N, by="TRT01AN"])
Ntot = copy(adsl[SAFFL=="Y", .N])

# Calculate percentages
final = final[, `:=` (C1 = percent(`1_81`, Ntrt[TRT01AN==81,N], 1),
                      C2 = percent(`1_54`, Ntrt[TRT01AN==54,N], 1),
                      C3 = percent(`1_0`, Ntrt[TRT01AN==0,N], 1),
                      C4 = percent(`2_N`, Ntot, 1)
)][order(ORD,AESOC,SUBORD)]
View(final)


