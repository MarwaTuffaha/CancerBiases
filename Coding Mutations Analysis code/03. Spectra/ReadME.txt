There are two types of spectra:
1- The one that comes from only counting the number of mutations for each mutation type (referred to as the counts-spectrum)
2- The one that comes from dividing the above counts by the number of all possible such mutation in the genome (referred to as the mutation-rate-spectrum)

Step 1: Find Cancer-Gene Spectra (CGS) and NonCancer-Gene Spectra (NCGS) in germline and PCAWG, hyper/nonhyper, syn/nonsyn
Matlab code: "CGS_SynNonsyn_MutRates" and "CGS_SynNonsyn_Counts"
Matlab Data: "CGS_NCGS_MutRates" and "CGS_NCGS_Counts"

Step 2: Find corresponding spectra in all genes. Also, find hyper/nonhyper spectra except each individual tissue samples
Matlab code: "Spectra_find_MutRates" and "Spectra_find_Counts"
Matlab Data: "Spectra_MutRates", "Spectra_ExceptTissue_MutRates" and "Spectra_Counts", "Spectra_ExceptTissue_Counts"

Step 3: Correlate CGS and all-gene spectra and plot
Matlab code: "Correlate_CDMS_Spectra_MutRates" and "Correlate_CDMS_Spectra_Counts"
Since there are generally high correlations between CDMS and all-gene spectra (except at very small number of mutations; not enough data), it seems like there's nothing special about cancer genes and we can focus on all-gene spectra for further analysis.

Step 4: All-gene spectra for Germline and tissues and plot
Matlab code: "Correlate_Tissues_Spectra_MutRates" and "Correlate_Tissues_Spectra_Counts"



