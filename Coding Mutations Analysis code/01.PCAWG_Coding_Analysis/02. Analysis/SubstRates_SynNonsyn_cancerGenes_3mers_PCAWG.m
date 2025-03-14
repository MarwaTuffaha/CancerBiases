% Substitution Rates_SynNonsyn_3mers_PCAWG Cancer genes

mers=load('3mers.mat');
U=mers.U;
n=length(U);

codon_muts=load('Path_To\Data_Extracted_files\codon_3mer_CancerGenes_synNonsyn');
C_syn_c=codon_muts.C_syn_c;
C_nonsyn_c=codon_muts.C_nonsyn_c;
C_syn_nc=codon_muts.C_syn_nc;
C_nonsyn_nc=codon_muts.C_nonsyn_nc;

PCAWG_muts=load('Path_To\Data_Extracted_files\PCAWG_3mer_synNonsyn_Tissues_cancerGenes.mat');
P_syn_c=sum(PCAWG_muts.P_syn_c,2);
P_syn_nc=sum(PCAWG_muts.P_syn_nc,2);
P_nonsyn_c=sum(PCAWG_muts.P_nonsyn_c,2);
P_nonsyn_nc=sum(PCAWG_muts.P_nonsyn_nc,2);

PCAWG_muts=load('Path_To\Data_Extracted_files\PCAWG_3mer_Nonsyn_Tissues_HyperNonhyper_cancerGenes.mat');
P_syn_c_hyper=sum(PCAWG_muts.P_syn_c_hyper,2);
P_syn_c_nonhyper=sum(PCAWG_muts.P_syn_c_nonhyper,2);
P_nonsyn_c_hyper=sum(PCAWG_muts.P_nonsyn_c_hyper,2);
P_nonsyn_c_nonhyper=sum(PCAWG_muts.P_nonsyn_c_nonhyper,2);
P_syn_nc_hyper=sum(PCAWG_muts.P_syn_nc_hyper,2);
P_syn_nc_nonhyper=sum(PCAWG_muts.P_syn_nc_nonhyper,2);
P_nonsyn_nc_hyper=sum(PCAWG_muts.P_nonsyn_nc_hyper,2);
P_nonsyn_nc_nonhyper=sum(PCAWG_muts.P_nonsyn_nc_nonhyper,2);

% check mutation rates for noncancer genes

syn_mutRate_nc=P_syn_nc./C_syn_nc;
nonsyn_mutRate_nc=P_nonsyn_nc./C_nonsyn_nc;
exp_nonsyn_nc=syn_mutRate_nc.*C_nonsyn_nc;
diff_nonsyn_nc=P_nonsyn_nc-exp_nonsyn_nc; 
reldiff_nonsyn_nc=diff_nonsyn_nc./exp_nonsyn_nc; % relative difference
% If difference is positive, the number of nonsyn mutations in PCAWG is lower than we expect 
% meaning there is possibly a suppression of nonsynonymous mutations in
% cancer

ind1=1:16; ind2=17:32; ind3=33:48; ind4=49:64; ind5=65:80; ind6=81:96;
% figure
% plot(1:96,exp_nonsyn,'o',1:96,P_nonsyn,'o')

% figure
% plot(1:96,diff_nonsyn,'o')
% yline(0)
% 
% figure
% plot(1:96,reldiff_nonsyn,'o')
% yline(0)

figure
plot(ind1,diff_nonsyn_nc(ind1),'o',ind2,diff_nonsyn_nc(ind2),'o',ind3,diff_nonsyn_nc(ind3),'o',ind4,diff_nonsyn_nc(ind4),'o',ind5,diff_nonsyn_nc(ind5),'o',ind6,diff_nonsyn_nc(ind6),'o','LineWidth',1)
yline(0)
title("Difference in number of nonsyn mutations in noncancer genes (cancer-expected)/expected")
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])


figure
plot(ind1,reldiff_nonsyn_nc(ind1),'o',ind2,reldiff_nonsyn_nc(ind2),'o',ind3,reldiff_nonsyn_nc(ind3),'o',ind4,reldiff_nonsyn_nc(ind4),'o',ind5,reldiff_nonsyn_nc(ind5),'o',ind6,reldiff_nonsyn_nc(ind6),'o','LineWidth',1)
yline(0)
title("Relative difference in number of nonsyn mutations in noncancer genes (cancer-expected)/expected")
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])

% Repeat for cancer genes

syn_mutRate_c=P_syn_c./C_syn_c;
nonsyn_mutRate_c=P_nonsyn_c./C_nonsyn_c;
exp_nonsyn_c=syn_mutRate_c.*C_nonsyn_c;
diff_nonsyn_c=P_nonsyn_c-exp_nonsyn_c; 
reldiff_nonsyn_c=diff_nonsyn_c./exp_nonsyn_c; % relative difference

figure
plot(ind1,diff_nonsyn_c(ind1),'o',ind2,diff_nonsyn_c(ind2),'o',ind3,diff_nonsyn_c(ind3),'o',ind4,diff_nonsyn_c(ind4),'o',ind5,diff_nonsyn_c(ind5),'o',ind6,diff_nonsyn_c(ind6),'o','LineWidth',1)
yline(0)
title("Difference in number of nonsyn mutations in cancer genes (cancer-expected)/expected")
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])


figure
plot(ind1,reldiff_nonsyn_c(ind1),'o',ind2,reldiff_nonsyn_c(ind2),'o',ind3,reldiff_nonsyn_c(ind3),'o',ind4,reldiff_nonsyn_c(ind4),'o',ind5,reldiff_nonsyn_c(ind5),'o',ind6,reldiff_nonsyn_c(ind6),'o','LineWidth',1)
yline(0)
title("Relative difference in number of nonsyn mutations in cancer genes (cancer-expected)/expected")
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])

%%

figure

a=[reldiff_nonsyn_c,reldiff_nonsyn_nc];
a(84,:)=[];
a(52,:)=[];
a(20,:)=[];
a(4,:)=[];
R= corrcoef(a);
corrCoeff=R(1,2);
plotting(a,"Rel diff in cancer genes","Rel diff in noncancer genes")

%% Repeat for hyper/nonhyper
syn_mutRate_c_hyper=P_syn_c_hyper./C_syn_c;
exp_nonsyn_c_hyper=syn_mutRate_c_hyper.*C_nonsyn_c;
diff_nonsyn_c_hyper=P_nonsyn_c_hyper-exp_nonsyn_c_hyper; 
reldiff_nonsyn_c_hyper=diff_nonsyn_c_hyper./exp_nonsyn_c_hyper; % relative difference

syn_mutRate_nc_hyper=P_syn_nc_hyper./C_syn_nc;
exp_nonsyn_nc_hyper=syn_mutRate_nc_hyper.*C_nonsyn_nc;
diff_nonsyn_nc_hyper=P_nonsyn_nc_hyper-exp_nonsyn_nc_hyper; 
reldiff_nonsyn_nc_hyper=diff_nonsyn_nc_hyper./exp_nonsyn_nc_hyper; % relative difference

syn_mutRate_c_nonhyper=P_syn_c_nonhyper./C_syn_c;
exp_nonsyn_c_nonhyper=syn_mutRate_c_nonhyper.*C_nonsyn_c;
diff_nonsyn_c_nonhyper=P_nonsyn_c_nonhyper-exp_nonsyn_c_nonhyper; 
reldiff_nonsyn_c_nonhyper=diff_nonsyn_c_nonhyper./exp_nonsyn_c_nonhyper; % relative difference

syn_mutRate_nc_nonhyper=P_syn_nc_nonhyper./C_syn_nc;
exp_nonsyn_nc_nonhyper=syn_mutRate_nc_nonhyper.*C_nonsyn_nc;
diff_nonsyn_nc_nonhyper=P_nonsyn_nc_nonhyper-exp_nonsyn_nc_nonhyper; 
reldiff_nonsyn_nc_nonhyper=diff_nonsyn_nc_nonhyper./exp_nonsyn_nc_nonhyper; % relative difference

figure
plot(ind1,diff_nonsyn_c_hyper(ind1),'o',ind2,diff_nonsyn_c_hyper(ind2),'o',ind3,diff_nonsyn_c_hyper(ind3),'o',ind4,diff_nonsyn_c_hyper(ind4),'o',ind5,diff_nonsyn_c_hyper(ind5),'o',ind6,diff_nonsyn_c_hyper(ind6),'o','LineWidth',1)
yline(0)
title("Difference in number of nonsyn mutations in cancer genes | Hyper | (cancer-expected)/expected")
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])

figure
plot(ind1,diff_nonsyn_c_nonhyper(ind1),'o',ind2,diff_nonsyn_c_nonhyper(ind2),'o',ind3,diff_nonsyn_c_nonhyper(ind3),'o',ind4,diff_nonsyn_c_nonhyper(ind4),'o',ind5,diff_nonsyn_c_nonhyper(ind5),'o',ind6,diff_nonsyn_c_nonhyper(ind6),'o','LineWidth',1)
yline(0)
title("Difference in number of nonsyn mutations in cancer genes | Nonhyper | (cancer-expected)/expected")
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])


figure
plot(ind1,reldiff_nonsyn_c_hyper(ind1),'o',ind2,reldiff_nonsyn_c_hyper(ind2),'o',ind3,reldiff_nonsyn_c_hyper(ind3),'o',ind4,reldiff_nonsyn_c_hyper(ind4),'o',ind5,reldiff_nonsyn_c_hyper(ind5),'o',ind6,reldiff_nonsyn_c_hyper(ind6),'o','LineWidth',1)
yline(0)
title("Relative difference in number of nonsyn mutations in cancer genes | Hyper | (cancer-expected)/expected")
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])

figure
plot(ind1,reldiff_nonsyn_c_nonhyper(ind1),'o',ind2,reldiff_nonsyn_c_nonhyper(ind2),'o',ind3,reldiff_nonsyn_c_nonhyper(ind3),'o',ind4,reldiff_nonsyn_c_nonhyper(ind4),'o',ind5,reldiff_nonsyn_c_nonhyper(ind5),'o',ind6,reldiff_nonsyn_c_nonhyper(ind6),'o','LineWidth',1)
yline(0)
title("Relative difference in number of nonsyn mutations in cancer genes | Nonhyper | (cancer-expected)/expected")
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])

figure
plot(ind1,reldiff_nonsyn_nc_hyper(ind1),'o',ind2,reldiff_nonsyn_nc_hyper(ind2),'o',ind3,reldiff_nonsyn_nc_hyper(ind3),'o',ind4,reldiff_nonsyn_nc_hyper(ind4),'o',ind5,reldiff_nonsyn_nc_hyper(ind5),'o',ind6,reldiff_nonsyn_nc_hyper(ind6),'o','LineWidth',1)
yline(0)
title("Relative difference in number of nonsyn mutations in noncancer genes | Hyper | (cancer-expected)/expected")
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])

figure
plot(ind1,reldiff_nonsyn_nc_nonhyper(ind1),'o',ind2,reldiff_nonsyn_nc_nonhyper(ind2),'o',ind3,reldiff_nonsyn_nc_nonhyper(ind3),'o',ind4,reldiff_nonsyn_nc_nonhyper(ind4),'o',ind5,reldiff_nonsyn_nc_nonhyper(ind5),'o',ind6,reldiff_nonsyn_nc_nonhyper(ind6),'o','LineWidth',1)
yline(0)
title("Relative difference in number of nonsyn mutations in noncancer genes | Nonhyper | (cancer-expected)/expected")
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])


function [] = plotting(a,b,c)

    R= corrcoef(a);
    s1 = plot(a(:,1), a(:,2), 'b+');
    set(s1, 'MarkerSize', 8, 'LineWidth', 2);
    s1 = ancestor(s1, 'axes');
    s1.YAxis.Exponent = 0;
    %%% regression line
    hold on
    l = lsline ;
    set(l,'LineWidth', 2)
    %%% axis display 
    xlabel(b, 'FontSize', 10)
    ylabel(c, 'FontSize', 10)
    set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')
    %t=text(min(a(:,1)),max(a(:,2)),{strcat('CorCoef = ',string(R(1,2)))});
    %t.FontSize=15;
    title(strcat("CorCoef=",string(R(1,2))))
    pbaspect([1 1 1])
end   



