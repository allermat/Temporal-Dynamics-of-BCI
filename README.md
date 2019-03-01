# Soruce code for research article: <br/>"To integrate or not to integrate:  Temporal dynamics of hierarchical Bayesian Causal Inference" 
**By Máté Aller and Uta Noppeney**

Article under review in Plos Biology. Preprint can be found [here](https://www.biorxiv.org/content/10.1101/504118v2).

## Description
- Class bci contains the methods for fitting Bayesian Causal Inference models to behavioural and EEG data
- Class mvpa contains the methods for performing Support Vector Regression and related statistics on EEG data
- Class mvpares handles the mvpa result data and plotting as an mvpares object. 
- The rest of the functions are for running the preprocessing, mvpa, bci, statistics and figure drawing, and some auxiliary functions. 

## Dependencies
Code is written in MATLAB, tested on MATLAB R2016a. 

1. [FieldTrip](http://www.fieldtriptoolbox.org/)
2. [SPM12](https://www.fil.ion.ucl.ac.uk/spm/)
3. [libsvm](https://www.csie.ntu.edu.tw/~cjlin/libsvm/)
4. [CircStat](https://www.jstatsoft.org/article/view/v031i10)
5. [MathWorks File Exchange functions](https://uk.mathworks.com/matlabcentral/fileexchange/?s_tid=gn_mlc_fx)
   - barwitherr
   - fminsearchbnd
   - jheapcl
   - numSubplots
   - parfor_progress
   - redblue
   - shadedErrorBar
   - sort_nat
   - suplabel

## Contributors
The code was mainly written by Máté Aller. Ulrik Beierholm contributed the Bayesian Causal Inference fitting code, Tim Rohe contributed the code for log-likelihood ratio statistic on circular indices (circ_AWTest.m). 
