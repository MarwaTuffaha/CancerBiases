Two functions are used here: 
comp_3mer.m	- Finds the complement 3mer if the middle ntd is not a pyrimidine
let2num.m	- Converts ntd's into numbers: A->1, C->2, G->3, T->4


Step 1:
Convert germline mutation data to numerical values (in Germline Dataset folder)
These are whole-genome mutations
Matlab code: "Export_Germline_Mutations_inNumbers"
Matlab data: "Germline_Mutations_Numbers"


Step 2:
Combine the above file with the codon content numerical file
Germline mutations not found in codon content are most probably noncoding -> ignored
First, fix the 3mer column to be similar to COSMIC (pyrimidine in the middle, otherwise find the complement)
Then classify 3mer mutations to syn/nonsyn
Matlab code: "Combining_Germline_Content"
Matlab Data: "Labeled_coding_Germline"


Step 3:
3mer level analysis
Matlab code: "analysis_Germline"
Matlab Data: "Germline_synNonsyn" and "Germline_synNonsyn_CancerGenes"
