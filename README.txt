Data and Code files to accomapny Gasch et al., Assessing impacts of crested wheatgrass and native
species establishment on soil characteristics in reclaimed land using Bayesian posterior predictive 
distributions, Land Degradation & Development.

Data File:
"CW_Data.csv"
Treatment: 
1 = native cool season grass seed mix; 
2 = crested wheatgrass seed mix; 
3 = undisturbed reference
Age: 
1 in Treatment 1 = 1.5 year old native cool season grass seed mix (not included in article);
2 in Treatment 1 = W14 (14 years old), native cool season grass seed mix;
3 in Treatment 1 = W26 (26 years old), native cool season grass seed mix;
1 in Treatment 2 = C11 (11 years old), crested wheatgrass seed mix;
2 in Treatment 2 = C16 (16 years old), crested wheatgrass seed mix;
3 in Treatment 2 = C29 (29 years old), crested wheatgrass seed mix
1 in Treatment 2 = UD, undisturbed reference area
Transect: 1-3 for each site (Treatment + Age combination)
Point: 1-4 for each transect
Depth: 1 = 0-5 cm; 2 = 5-15 cm soil depth
Macro: macroaggregate weight (g aggregate/g soil)
Micro: microaggregate weight (g aggregate/g soil)
CN: C:N ratio (unitless) of aggregates (macro + micro)
Total: total microbial abundance, determined by PLFA (nmol FA/g soil)
Eubac: bacteria abundance, determined by PLFA (nmol FA/g soil)
Actino: actinomycete abundance, determined by PLFA (nmol FA/g soil)
AMF: arbuscular mycorrhizal fungi abundance, determined by PLFA (nmol FA/g soil)
Fungi: fungi abundance, determined by PLFA (nmol FA/g soil)

Model Files:
Four model files are provided, each for a different scenario:
"cond_bugs_bvn_model.txt" for a bivariate model when one of the bivariate observations is missing, to be run by BUGS
"cond_bugs_bvn_model_ln.txt" for a log bivariate model when one of the bivariate observations is missing, to be run by BUGS
"cond_bvn_model.txt" for a bivariate model with complete observations, to be run by JAGS
"cond_bvn_model_ln.txt" for a log bivariate model with complete observations, to be run by JAGS

R Code File:
"cond_spec_mix_Total.R"
The script provided is for analysis of the Total microbial abundance, which requires use of each of the four model scenarios,
and therefore should provide an adaptable example of each scenario encountered within this data set. Script also includes
code for producing plots included in the article.


