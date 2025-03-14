Step 1: ((Preparing data))

Three functions are used here: 
effect2num.m	- Converts mutation effect to a numerical value: SILENT->0, MISSENSE->1, NONSENSE->2
comp_3mer.m	- Finds the complement 3mer if the middle ntd is not a pyrimidine
let2num.m	- Converts ntd's into numbers: A->1, C->2, G->3, T->4

-Converting genome coding content data into numerical form to reduce data memory size
and extract mutations only in cancer genes and those only in passenger genes
Matlab code: "Export_Codon_Mutations_inNumbers"
Matlab data: "Codon_Content_Mutations" and "Codon_Content_Mutations_3merFixed" and "Codon_Content_CancerGenes_Mutations"

-Combining all PCAWG coding mutations in one file and save numerically
Fix the 3mer column to be similar to COSMIC (pyrimidine in the middle, otherwise find the complement)
Matlab code: "Combining_PCAWG_Coding_Mutations"
Matlab data: "PCAWG_Coding_Mutations" and "PCAWG_coding_counts"

-Combining the above two files in one file containing all information about each mutation
Matlab code: "Combining_PCAWG_Content"
Matlab data: "Labeled_Coding_PCAWG"
Remark: Some PCAWG mutations were not found in the codon content big file ( # = 2022 , fraction = 0.0048 ) - Deleted


Step 2: ((Analysis))

-Lowest level analysis - syn/nonsyn
Matlab code: "analysis1_SynNonsyn"
Matlab Data: "tissue_mutations_count" (how many mutations are there for each tissue type)

A-3mer level analysis
Count syn/nonsyn - Ts/Tv for the 96 3mer mutations in PCAWG
and more detailed classification according to each tissue type
Matlab code: "analysis2_3merLevel"
Matlab data: "3mers" and "PCAWG_3mer_synNonsyn" and "PCAWG_3mer_synNonsyn_tissues" and "Labeled_Coding_PCAWG_hyperNonhyper" and "PCAWG_TvTs_hyperNonhyper" and "PCAWG_3mer_synNonsyn_hyperNonhyper"

B-Count syn/nonsyn - Ts/Tv for the 96 3mer mutations in the coding genome
Also, counting mutations in each tissue type (all together - regardless of 3mers)
Matlab code: "Count_3mer_Codon_Mutations_SynNonsyn"
Matlab data: "codon_3mer_synNonsyn"

C-Combine the above two results to find mutation rates of syn/nonsyn in PCAWG data
Matlab code: "SubstRates_SynNonsyn_3mers_PCAWG"

 
-Repeat steps A,B,C but label mutations according to being in ((cancer genes)) or not too
Goal is to see if result change is we delete mutations occuring in cancer genes (positively selected)
Matlab code1: "analysis2_CancerGenes"
Matlab code2: "Count_3mer_Codon_CancerGenes_SynNonsyn"
Matlab code3: "SubstRates_SynNonsyn_cancerGenes_3mers_PCAWG"
Matlab data: "PCAWG_3mer_synNonsyn_Tissues_cancerGenes", "PCAWG_3mer_Nonsyn_Tissues_HyperNonhyper_cancerGenes" 
		"codon_3mer_CancerGenes_synNonsyn" and "Count_3mer_mutations_tissues"




-Repeat again but:
On average, about 5% of mutations are in cancer genes.
Choose 5% random mutations and check if you get the same effect; compare to the positive selection effect seen in cancer gene mutations

-To get over the fact that there aren't any examples of syn mutations for 4 3mer mutations in the human genome
Instead of using the syn mutation rate as the neutral, use all coding mutations rate in noncancer genes as neutral.

Reduce skin and colon HM sample mutations to avoid dominance of these two tissues in the pooled HM mutation spectrum
Matlab code1: "Reduce_skin_colon_HMmutations"
MATLAB data: "Reduced_SkinColonHM"
Matlab code2: "SubstRates_SynNS_canGenes_3mer_PCAWG_Method2_SkinColonReduced"
MATLAB code3: "analysis2_Random5Percent_SkinColonReduced"
MATLAB data: "PCAWG_3mer_synNonsyn_Random5Percent_SkinColonReduced"
MATLAB code4: "SubstRates_SynNS_Rand5Perc_3mers_PCAWG_Method2_SkinColonReduced"


-Back to lowest level analysis - Ts/Tv extracted from 3mer analysis results
Matlab code: "analysis3_TvTs"
Matlab data: "TsTv_coding_result"


Now we have all the counts we need to extract 3mer spectra and compare
