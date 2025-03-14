% MutationRates_SynNonsyn_3mers_PCAWG

mers=load('3mers.mat');
U=mers.U;
n=length(U);

codon_muts=load('Path_To\Data_Extracted_files\codon_3mer_synNonsyn');
C_syn=codon_muts.C_syn;
C_nonsyn=codon_muts.C_nonsyn;

PCAWG_muts=load('Peth_To\Data_Extracted_files\PCAWG_3mer_synNonsyn');
P_syn=PCAWG_muts.P_syn;
P_nonsyn=PCAWG_muts.P_nonsyn;

syn_mutRate=P_syn./C_syn;
nonsyn_mutRate=P_nonsyn./C_nonsyn;
exp_nonsyn=syn_mutRate.*C_nonsyn;
diff_nonsyn=P_nonsyn-exp_nonsyn; 
reldiff_nonsyn=diff_nonsyn./exp_nonsyn; % relative difference
% If difference is positive, the number of nonsyn mutations in PCAWG is lower than we expect 
% meaning there is possibly a suppression of nonsynonymous mutations in
% cancer

ind1=1:16; ind2=17:32; ind3=33:48; ind4=49:64; ind5=65:80; ind6=81:96;
positives=length(find(diff_nonsyn>0));
detailed_positives=[length(find(diff_nonsyn(ind1)>0)),length(find(diff_nonsyn(ind2)>0)),length(find(diff_nonsyn(ind3)>0)),length(find(diff_nonsyn(ind4)>0)),length(find(diff_nonsyn(ind5)>0)),length(find(diff_nonsyn(ind6)>0))];
frac_positives=detailed_positives/16;
sum_positives=[sum(reldiff_nonsyn(find(diff_nonsyn(ind1)>0))),sum(reldiff_nonsyn(16+find(diff_nonsyn(ind2)>0))),sum(reldiff_nonsyn(16*2+find(diff_nonsyn(ind3)>0))),sum(reldiff_nonsyn(16*3+find(diff_nonsyn(ind4)>0))),sum(reldiff_nonsyn(16*4+find(diff_nonsyn(ind5)>0))),sum(reldiff_nonsyn(16*5+find(diff_nonsyn(ind6)>0)))];
sum_negatives=[sum(reldiff_nonsyn(find(diff_nonsyn(ind1)<0))),sum(reldiff_nonsyn(16+find(diff_nonsyn(ind2)<0))),sum(reldiff_nonsyn(16*2+find(diff_nonsyn(ind3)<0))),sum(reldiff_nonsyn(16*3+find(diff_nonsyn(ind4)<0))),sum(reldiff_nonsyn(16*4+find(diff_nonsyn(ind5)<0))),sum(reldiff_nonsyn(16*5+find(diff_nonsyn(ind6)<0)))];

ts_reldiff_nonsyn=[reldiff_nonsyn(ind3);reldiff_nonsyn(ind5)];
tv_reldiff_nonsyn=[reldiff_nonsyn(ind1);reldiff_nonsyn(ind2);reldiff_nonsyn(ind4);reldiff_nonsyn(ind6)];

figure
plot(ind1,diff_nonsyn(ind1),'o',ind2,diff_nonsyn(ind2),'o',ind3,diff_nonsyn(ind3),'o',ind4,diff_nonsyn(ind4),'o',ind5,diff_nonsyn(ind5),'o',ind6,diff_nonsyn(ind6),'o','LineWidth',1)
yline(0)
title("Difference in # of nonsyn mutations (cancer-expected)")
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])


figure
plot(ind1,reldiff_nonsyn(ind1),'o',ind2,reldiff_nonsyn(ind2),'o',ind3,reldiff_nonsyn(ind3),'o',ind4,reldiff_nonsyn(ind4),'o',ind5,reldiff_nonsyn(ind5),'o',ind6,reldiff_nonsyn(ind6),'o','LineWidth',1)
yline(0)
title("Relative difference in number of nonsyn mutations (cancer-expected)/expected")
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])


%% Now give mutations from each tissue type same weight to avoid biases

PCAWG_muts_t=load('Path_To\Data_Extracted_files\PCAWG_3mer_synNonsyn_tissues.mat');
Pt_syn=PCAWG_muts_t.counts_syn;
Pt_nonsyn=PCAWG_muts_t.counts_nonsyn;
Pt_all=Pt_nonsyn+Pt_syn;

syn_mutRate_t=Pt_syn./C_syn;
mean_syn_mutRate_t=mean(syn_mutRate_t,2);
exp_nonsyn_t=mean_syn_mutRate_t.*C_nonsyn;
mean_Pt_nonsyn=mean(Pt_nonsyn,2);
diff_nonsyn_t=mean_Pt_nonsyn-exp_nonsyn_t; 

figure
plot(ind1,diff_nonsyn_t(ind1),'o',ind2,diff_nonsyn_t(ind2),'o',ind3,diff_nonsyn_t(ind3),'o',ind4,diff_nonsyn_t(ind4),'o',ind5,diff_nonsyn_t(ind5),'o',ind6,diff_nonsyn_t(ind6),'o','LineWidth',1)
yline(0)
title(["Difference in # of nonsyn mutations (cancer-expected)"; "normalized by tissue type"])
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])


% cannot use relative difference because of many zero counts for synonymous
% The following code gives a figure that shows only 3mers with no such
% problems
reldiff_nonsyn_t=diff_nonsyn_t./exp_nonsyn_t; % relative difference


figure
plot(ind1,reldiff_nonsyn_t(ind1),'o',ind2,reldiff_nonsyn_t(ind2),'o',ind3,reldiff_nonsyn_t(ind3),'o',ind4,reldiff_nonsyn_t(ind4),'o',ind5,reldiff_nonsyn_t(ind5),'o',ind6,reldiff_nonsyn_t(ind6),'o','LineWidth',1)
yline(0)
title(["Relative Difference in # of nonsyn mutations (cancer-expected)/expected"; "normalized by tissue type"])
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])

%% Instead, try removing the colon and skin (indeces=7,18)

Pt_syn(:,7)=[];
Pt_syn(:,18)=[];
Pt_nonsyn(:,7)=[];
Pt_nonsyn(:,18)=[];

Pt_syn_r=sum(Pt_syn,2);
Pt_nonsyn_r=sum(Pt_nonsyn,2);

syn_mutRate_r=Pt_syn_r./C_syn;
nonsyn_mutRate_r=Pt_nonsyn_r./C_nonsyn;
exp_nonsyn_r=syn_mutRate_r.*C_nonsyn;
diff_nonsyn_r=Pt_nonsyn_r-exp_nonsyn_r; 
reldiff_nonsyn_r=diff_nonsyn_r./exp_nonsyn_r; % relative difference

figure
plot(ind1,diff_nonsyn_r(ind1),'o',ind2,diff_nonsyn_r(ind2),'o',ind3,diff_nonsyn_r(ind3),'o',ind4,diff_nonsyn_r(ind4),'o',ind5,diff_nonsyn_r(ind5),'o',ind6,diff_nonsyn_r(ind6),'o','LineWidth',1)
yline(0)
title(["Difference in # of nonsyn mutations (cancer-expected)"; "skin and colon removed"])
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])


figure
plot(ind1,reldiff_nonsyn_r(ind1),'o',ind2,reldiff_nonsyn_r(ind2),'o',ind3,reldiff_nonsyn_r(ind3),'o',ind4,reldiff_nonsyn_r(ind4),'o',ind5,reldiff_nonsyn_r(ind5),'o',ind6,reldiff_nonsyn_r(ind6),'o','LineWidth',1)
yline(0)
title(["Relative difference in number of nonsyn mutations (cancer-expected)/expected";"skin and colon removed"])
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])










