library(data.table)
library(foreign)

adsl = data.table(read.xport("./adsl.xpt"))
adae = data.table(read.xport("./adae_j.xpt")
       )[SAFFL=="Y" & !is.na(AEDECOD) & TRTEMFL=="Y"]

###############################################################################
# Derivations
###############################################################################
# Required variables
vars = c("USUBJID", "TRT01AN", "TRT01A", "AESOC", "AESOCJ", "AEDECOD", "AEDECODJ", "TRTEMFL")
grd <- merge(adsl, adae)[, ..vars]

# Function to count per arm and period and transpose the counts
cnt = function(trt, col) {
  grd = grd[TRT01AN %in% trt]
  
  # For all columns except total patients
  grd1 = unique(grd[, .(USUBJID, TRT01AN, TRT01A)])
  cnt1 = grd1[, .N, by = c("TRT01AN", "TRT01A")
         ][, `:=` (AEDECOD = "Number of subjects with at least one event", AESOC = "AAA", ORD = 1)]
  
  grd2 = unique(grd[, .(USUBJID, TRT01AN, TRT01A, AESOC, AESOCJ)])
  cnt2 = grd2[, .N, by = c("TRT01AN", "TRT01A", "AESOC", "AESOCJ")][, ORD:=2]
  
  grd3 = unique(grd[,.(USUBJID, TRT01AN, TRT01A, AESOC, AESOCJ, AEDECOD, AEDECODJ)])
  cnt3 = grd3[, .N,by = c("TRT01AN", "TRT01A", "AESOC", "AESOCJ", "AEDECOD", "AEDECODJ")][, ORD:=2]
  
  combine = data.table(rbind(cnt1, cnt2, cnt3, fill = TRUE)
            )[, `:=` (SUBORD   = ifelse(is.na(AEDECOD), 1, 2),
                      AEDECOD  = ifelse(is.na(AEDECOD),  AESOC,  AEDECOD),
                      AEDECODJ = ifelse(is.na(AEDECODJ), AESOCJ, AEDECODJ))
            ][order(AESOC)]
  
  trans_dset = dcast(combine, ORD + SUBORD + AESOC + AESOCJ + AEDECOD + AEDECODJ ~ TRT01AN, value.var = "N")
  colnames(trans_dset) = paste(as.character(col), colnames(trans_dset), sep = "_")
  
  return(trans_dset)
}

o_all = cnt(trt = c(54, 81, 0), col = 1) # Overall study period- all arms