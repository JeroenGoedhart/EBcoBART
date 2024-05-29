# EBcoBART
 EBcoBART R package
 R package to estimate prior variable weights (the probabilities that variables are selected in the splitting rules) for
 Bayesian Additive Regression Trees (BART). These prior variable weights are estimated using empirical Bayes and
 auxiliary information on the variables (termed co-data). For all details on the method see:

 https://arxiv.org/abs/2311.09997

 # Installation
#install.packages("remotes")
remotes::install_github("JeroenGoedhart/EBcoBART")
