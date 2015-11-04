# This script extracts some relevant quantities for the ipython notebook example.ipynb
# from the Natsal-3 dataset.

rm(list=ls())

library(foreign)
library(questionr)
library(survey)

natsal3 <- read.dta("~/Documents/data/natsal-3/UKDA-7799-stata11/stata11/eul_natsal_2010_for_archive.dta") # 15162 rows

# make a variable for whether sexually active (same-sex or opposite-sex)
natsal3$active <- factor(NA, levels=c("yes", "no"))
natsal3$active[natsal3$everhet == "Heterosexual sex since age 13"] <- "yes"
natsal3$active[natsal3$eversam == "Yes"] <- "yes"
natsal3$active[(natsal3$everhet == "No heterosexual sex") & (natsal3$eversam == "No")] <- "no"

# set up survey design
n3design <- svydesign( 
		id = ~psu_scrm , 
		strata = ~strata ,
		data = natsal3 ,		
		weight = ~total_wt 
	)

# point estimate for proportion sexually active
svyby(~active, by= ~rsex+agrp2, n3design, svymean, na.rm=TRUE)
svyby(~active, by= ~rsex+agrp, n3design, svymean, na.rm=TRUE)

# 95% confidence interval for proportion sexually active
confint(svyby(~active, by= ~rsex+agrp2, n3design, svymean, na.rm=TRUE))
confint(svyby(~active, by= ~rsex+agrp, n3design, svymean, na.rm=TRUE))