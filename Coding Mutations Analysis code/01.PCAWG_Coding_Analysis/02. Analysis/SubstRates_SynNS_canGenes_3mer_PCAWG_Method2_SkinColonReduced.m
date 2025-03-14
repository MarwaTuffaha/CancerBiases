%To get over the fact that there aren't any examples of syn mutations for 4 3mer mutations in the human genome
%Instead of using the syn mutation rate as the neutral, use all coding mutations rate in noncancer genes as neutral.

mers=load('3mers.mat');
U=mers.U;
n=length(U);

codon_muts=load('Path_To\Data_Extracted_files\codon_3mer_CancerGenes_synNonsyn');
C_syn_c=codon_muts.C_syn_c;
C_nonsyn_c=codon_muts.C_nonsyn_c;
C_syn_nc=codon_muts.C_syn_nc;
C_nonsyn_nc=codon_muts.C_nonsyn_nc;
C_nc=C_syn_nc+C_nonsyn_nc;

load('Path_To\Data_Extracted_files\PCAWG_3mer_Nonsyn_Tissues_HyperNonhyper_cancerGenes.mat',"P_syn_c_hyper","P_syn_c_nonhyper","P_nonsyn_c_hyper","P_nonsyn_c_nonhyper","P_syn_nc_hyper","P_syn_nc_nonhyper","P_nonsyn_nc_hyper","P_nonsyn_nc_nonhyper")

P_syn_c_nonhyper=sum(P_syn_c_nonhyper,2);
P_nonsyn_c_nonhyper=sum(P_nonsyn_c_nonhyper,2);
P_syn_nc_nonhyper=sum(P_syn_nc_nonhyper,2);
P_nonsyn_nc_nonhyper=sum(P_nonsyn_nc_nonhyper,2);

% Replace skin and colon HM counts with their reduced versions
load('Path_To\Data_Extracted_files\Reduced_SkinColonHM.mat',"P_colon_HM_nonsyn_nc","P_colon_HM_syn_nc","P_colon_HM_nonsyn_c","P_colon_HM_syn_c","P_skin_HM_nonsyn_nc","P_skin_HM_syn_nc","P_skin_HM_nonsyn_c","P_skin_HM_syn_c")

P_syn_c_hyper(:,7)=P_colon_HM_syn_c;
P_syn_nc_hyper(:,7)=P_colon_HM_syn_nc;
P_nonsyn_c_hyper(:,7)=P_colon_HM_nonsyn_c;
P_nonsyn_nc_hyper(:,7)=P_colon_HM_nonsyn_nc;

P_syn_c_hyper(:,18)=P_skin_HM_syn_c;
P_syn_nc_hyper(:,18)=P_skin_HM_syn_nc;
P_nonsyn_c_hyper(:,18)=P_skin_HM_nonsyn_c;
P_nonsyn_nc_hyper(:,18)=P_skin_HM_nonsyn_nc;

% now continue as before

P_syn_c_hyper=sum(P_syn_c_hyper,2);
P_nonsyn_c_hyper=sum(P_nonsyn_c_hyper,2);
P_syn_nc_hyper=sum(P_syn_nc_hyper,2);
P_nonsyn_nc_hyper=sum(P_nonsyn_nc_hyper,2);

P_nc_hyper=P_syn_nc_hyper+P_nonsyn_nc_hyper;
P_nc_nonhyper=P_syn_nc_nonhyper+P_nonsyn_nc_nonhyper;

P_c_hyper=P_syn_c_hyper+P_nonsyn_c_hyper;
P_c_nonhyper=P_syn_c_nonhyper+P_nonsyn_c_nonhyper;

P_hyper=P_nc_hyper+P_c_hyper;
P_nonhyper=P_nc_nonhyper+P_c_nonhyper;

ind1=1:16; ind2=17:32; ind3=33:48; ind4=49:64; ind5=65:80; ind6=81:96;

% Find mutation rates for noncancer genes to be used as the neutral
% expectation

mutRate_nc_hyper=P_nc_hyper./C_nc;
mutRate_nc_nonhyper=P_nc_nonhyper./C_nc;

exp_syn_c_hyper=mutRate_nc_hyper.*C_syn_c;
exp_nonsyn_c_hyper=mutRate_nc_hyper.*C_nonsyn_c;

diff_syn_c_hyper=P_syn_c_hyper-exp_syn_c_hyper; 
diff_nonsyn_c_hyper=P_nonsyn_c_hyper-exp_nonsyn_c_hyper; 
reldiff_syn_c_hyper=diff_syn_c_hyper./exp_syn_c_hyper;
reldiff_nonsyn_c_hyper=diff_nonsyn_c_hyper./exp_nonsyn_c_hyper;

exp_syn_c_nonhyper=mutRate_nc_nonhyper.*C_syn_c;
exp_nonsyn_c_nonhyper=mutRate_nc_nonhyper.*C_nonsyn_c;

diff_syn_c_nonhyper=P_syn_c_nonhyper-exp_syn_c_nonhyper; 
diff_nonsyn_c_nonhyper=P_nonsyn_c_nonhyper-exp_nonsyn_c_nonhyper; 
reldiff_syn_c_nonhyper=diff_syn_c_nonhyper./exp_syn_c_nonhyper;
reldiff_nonsyn_c_nonhyper=diff_nonsyn_c_nonhyper./exp_nonsyn_c_nonhyper;


figure
plot(ind1,reldiff_syn_c_hyper(ind1),'o',ind2,reldiff_syn_c_hyper(ind2),'o',ind3,reldiff_syn_c_hyper(ind3),'o',ind4,reldiff_syn_c_hyper(ind4),'o',ind5,reldiff_syn_c_hyper(ind5),'o',ind6,reldiff_syn_c_hyper(ind6),'o','LineWidth',1)
yline(0)
title("Relative difference in number of Syn mutations in cancer genes in Hypermutators (cancer-expected)/expected")
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])


figure
plot(ind1,reldiff_nonsyn_c_hyper(ind1),'o',ind2,reldiff_nonsyn_c_hyper(ind2),'o',ind3,reldiff_nonsyn_c_hyper(ind3),'o',ind4,reldiff_nonsyn_c_hyper(ind4),'o',ind5,reldiff_nonsyn_c_hyper(ind5),'o',ind6,reldiff_nonsyn_c_hyper(ind6),'o','LineWidth',1)
yline(0)
title("Relative difference in number of nonsyn mutations in cancer genes in Hypermutators (cancer-expected)/expected")
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])



figure
plot(ind1,reldiff_syn_c_nonhyper(ind1),'o',ind2,reldiff_syn_c_nonhyper(ind2),'o',ind3,reldiff_syn_c_nonhyper(ind3),'o',ind4,reldiff_syn_c_nonhyper(ind4),'o',ind5,reldiff_syn_c_nonhyper(ind5),'o',ind6,reldiff_syn_c_nonhyper(ind6),'o','LineWidth',1)
yline(0)
title("Relative difference in number of Syn mutations in cancer genes in Nonhypermutators (cancer-expected)/expected")
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])


figure
plot(ind1,reldiff_nonsyn_c_nonhyper(ind1),'o',ind2,reldiff_nonsyn_c_nonhyper(ind2),'o',ind3,reldiff_nonsyn_c_nonhyper(ind3),'o',ind4,reldiff_nonsyn_c_nonhyper(ind4),'o',ind5,reldiff_nonsyn_c_nonhyper(ind5),'o',ind6,reldiff_nonsyn_c_nonhyper(ind6),'o','LineWidth',1)
yline(0)
title("Relative difference in number of nonsyn mutations in cancer genes in Nonhypermutators (cancer-expected)/expected")
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])



save("Path_To\Data_Extracted_files\RelDiff_CancerGenes_Method2_SkinColonHMReduced","reldiff_syn_c_hyper","reldiff_nonsyn_c_hyper","reldiff_syn_c_nonhyper","reldiff_nonsyn_c_nonhyper","diff_syn_c_hyper","diff_nonsyn_c_hyper","diff_syn_c_nonhyper","diff_nonsyn_c_nonhyper")









