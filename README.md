# Analysing chlamydia surveillance data for England.

This repository contains software for analysing chlamydia surveillance data, applied to data from England in 2012. Code consists of R scrips, STAN model files and Jupyter notebooks in the IPython language.

For the reader's convenience, the Jupyter notebooks have been downloaded as both LaTeX/PDF and html files, and the code and Figures can be viewed using a PDF reader or web browser. To run or edit the notebooks and generate figures, the necessary software to run IPython and Jupyter notebooks must be installed (see http://ipython.org/install.html and http://jupyter.org/). The material

## Introducing the model

Files in the directory named **three_compartment_model** introduce the model and its properties.

## Example: chlamydia in England, 2012

Files in the directory named **england** illustrate the use of the model, with data from England in 2012.

As well as the IPython notebook, the directory contains an R script which is used to derive prior distributions for the proportions of men and women in different age groups who are sexually active, using data from the National Study of Sexual Attitudes and Lifestyles, Natsal-3.

The subdirectory named 'stan' contains STAN model files and R scripts to run them and obtain samples for natural clearance rates of chlamydia infection in men and women, based on the work of Price _et al._ presented in _Stat Med_ **32**:1547-1560 (2013).

## Local differences in chlamydia incidence, testing, diagnosis and prevalence

Files in the directory named **local_authorities.ipynb** use the model to investigate local differences in chlamydia testing and diagnosis rates, incidence and prevalence.

## Data

Natsal-3 data is available from the UK Data Archive:
http://www.data-archive.ac.uk/

Chlamydia testing data was originally downloaded from:
http://www.chlamydiascreening.nhs.uk/ps/data.asp
Checked 9 February 2016

English indices of deprivation 2010 downloaded from:
https://www.gov.uk/government/statistics/english-indices-of-deprivation-2010
Contains public sector information licensed under the Open Government Licence v3.0; http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/
Checked 9 February 2016

Shape files and coding for English local authorities were downloaded from:
https://geoportal.statistics.gov.uk/geoportal/catalog/content/filelist.page
Contains National Statistics data © Crown copyright and database right 2016
