
# load testdata and space pre-allocation
dat <- read.csv("TESTDATA.csv")
n_test <- nrow(dat)
sas <- matrix(NA, n_test,2)
v1 <- matrix(NA, n_test,2)
# v2 <- matrix(NA, n_test,2)

# load sas result from .lst file
sas_lst <- readLines("bscCI_test.lst")
offset_bscline <- 20
counter <- 1
ii <- 1
while (counter < length(sas_lst)) {
  txt <- sas_lst[counter]
  if (grepl("I=",txt)) {
    txt_bsc <- sas_lst[counter+offset_bscline]
    txt_bsc <- sub("\\(Blyth-Still-Casella\\)  = \\(", "", txt_bsc)
    txt_bsc <- sub("\\)", "", txt_bsc)
    txt_bsc <- strsplit(txt_bsc, ",")[[1]]
    txt_bsc <- trimws(txt_bsc)
    sas[ii,] <- as.numeric(txt_bsc)
    ii <- ii+1
    counter <- counter+offset_bscline
  }
  else {
    counter <- counter+1
  }
}

# R results
library(fastCI)
system.time(
  for (ii in 1:n_test) {
    v1[ii,] <- bscCI(dat$n[ii], dat$x[ii], 0.95)
  }
)

# library(BlythStillCasellaCI)
# system.time(
#   for (ii in 1:n_test) {
#     v2[ii,] <- blyth.still.casella(dat$n[ii],dat$x[ii], digits = 4, 0.05)
#   }
# )

# put the results together
# out_compare <- data.frame(sas = sas, v1 = v1, v2 = v2)
out_compare <- data.frame(sas = sas, v1 = v1)
# View(out_compare)

# compare
# bscCI
max(out_compare$sas.1-out_compare$v1.1) # 4.842831e-05
max(out_compare$sas.2-out_compare$v1.2) # 4.641143e-05

# BlythStillCasellaCI
# max(out_compare$sas.1-out_compare$v2.1) # 0.0081
# max(out_compare$sas.2-out_compare$v2.2) # 0.0106

# --- test some util ----
system.time({
  for (ii in 1:1000) {
    bicoln_raw(2000,ii)
  }
})

system.time({
  for (ii in 1:1000) {
    bicoln_mem(2000,ii)
  }
})

system.time({
  bbb_pvalue(100,120,1e-6,50,50,0.01)
})

system.time({
  bbb_fast_pvalue(100,120,1e-6,50,50,0.01,0.95)
})

