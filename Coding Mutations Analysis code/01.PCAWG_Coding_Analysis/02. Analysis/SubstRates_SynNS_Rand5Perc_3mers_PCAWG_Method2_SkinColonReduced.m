% Compare bootstrapped to cancer-gene mutations
% to show that cancer-gene mutations are special (excess of nonsyn
% mutations in NHM)

load('Path_to\Data_Extracted_files\PCAWG_3mer_Nonsyn_Tissues_HyperNonhyper_cancerGenes.mat',"P_syn_c_hyper","P_syn_c_nonhyper","P_nonsyn_c_hyper","P_nonsyn_c_nonhyper","P_syn_nc_hyper","P_syn_nc_nonhyper","P_nonsyn_nc_hyper","P_nonsyn_nc_nonhyper")

P_syn_c_nonhyper=sum(P_syn_c_nonhyper,2);
P_nonsyn_c_nonhyper=sum(P_nonsyn_c_nonhyper,2);
P_syn_nc_nonhyper=sum(P_syn_nc_nonhyper,2);
P_nonsyn_nc_nonhyper=sum(P_nonsyn_nc_nonhyper,2);

% Replace skin and colon HM counts with their reduced versions
load('Path_to\Data_Extracted_files\Reduced_SkinColonHM.mat',"P_colon_HM_nonsyn_nc","P_colon_HM_syn_nc","P_colon_HM_nonsyn_c","P_colon_HM_syn_c","P_skin_HM_nonsyn_nc","P_skin_HM_syn_nc","P_skin_HM_nonsyn_c","P_skin_HM_syn_c")

P_syn_c_hyper(:,7)=P_colon_HM_syn_c;
P_syn_nc_hyper(:,7)=P_colon_HM_syn_nc;
P_nonsyn_c_hyper(:,7)=P_colon_HM_nonsyn_c;
P_nonsyn_nc_hyper(:,7)=P_colon_HM_nonsyn_nc;
P_syn_c_hyper(:,18)=P_skin_HM_syn_c;
P_syn_nc_hyper(:,18)=P_skin_HM_syn_nc;
P_nonsyn_c_hyper(:,18)=P_skin_HM_nonsyn_c;
P_nonsyn_nc_hyper(:,18)=P_skin_HM_nonsyn_nc;

P_syn_c_hyper=sum(P_syn_c_hyper,2);
P_nonsyn_c_hyper=sum(P_nonsyn_c_hyper,2);
P_syn_nc_hyper=sum(P_syn_nc_hyper,2);
P_nonsyn_nc_hyper=sum(P_nonsyn_nc_hyper,2);


load('Path_to\Data_Extracted_files\PCAWG_3mer_synNonsyn_Random5Percent_SkinColonReduced','B_syn_hyper','B_nonsyn_hyper','B_syn_nonhyper','B_nonsyn_nonhyper')

% assuming the distributions come from normal distribution
% and using Bonferroni correction because we test 96 values
% We look for confidence interval with alpha=5%/96=0.052%
% The corresponding z-score is 3.47

z=3.47;

mean_syn_hyper=mean(B_syn_hyper,2); sd_syn_hyper=std(B_syn_hyper,0,2);
mean_nonsyn_hyper=mean(B_nonsyn_hyper,2); sd_nonsyn_hyper=std(B_nonsyn_hyper,0,2);
mean_syn_nonhyper=mean(B_syn_nonhyper,2); sd_syn_nonhyper=std(B_syn_nonhyper,0,2);
mean_nonsyn_nonhyper=mean(B_nonsyn_nonhyper,2); sd_nonsyn_nonhyper=std(B_nonsyn_nonhyper,0,2);

pos_syn_hyper=[]; inside_syn_hyper=[]; neg_syn_hyper=[];
pos_nonsyn_hyper=[]; inside_nonsyn_hyper=[]; neg_nonsyn_hyper=[];
pos_syn_nonhyper=[]; inside_syn_nonhyper=[]; neg_syn_nonhyper=[];
pos_nonsyn_nonhyper=[]; inside_nonsyn_nonhyper=[]; neg_nonsyn_nonhyper=[];


for i=1:96
    if P_syn_c_hyper(i)<mean_syn_hyper(i)-z*sd_syn_hyper(i)
        neg_syn_hyper=[neg_syn_hyper,i];
    elseif P_syn_c_hyper(i)>mean_syn_hyper(i)+z*sd_syn_hyper(i)
        pos_syn_hyper=[pos_syn_hyper,i];
    else
        inside_syn_hyper=[inside_syn_hyper,i];
    end

    if P_nonsyn_c_hyper(i)<mean_nonsyn_hyper(i)-z*sd_nonsyn_hyper(i)
        neg_nonsyn_hyper=[neg_nonsyn_hyper,i];
    elseif P_nonsyn_c_hyper(i)>mean_nonsyn_hyper(i)+z*sd_nonsyn_hyper(i)
        pos_nonsyn_hyper=[pos_nonsyn_hyper,i];
    else
        inside_nonsyn_hyper=[inside_nonsyn_hyper,i];
    end

    if P_syn_c_nonhyper(i)<mean_syn_nonhyper(i)-z*sd_syn_nonhyper(i)
        neg_syn_nonhyper=[neg_syn_nonhyper,i];
    elseif P_syn_c_nonhyper(i)>mean_syn_nonhyper(i)+z*sd_syn_nonhyper(i)
        pos_syn_nonhyper=[pos_syn_nonhyper,i];
    else
        inside_syn_nonhyper=[inside_syn_nonhyper,i];
    end

    if P_nonsyn_c_nonhyper(i)<mean_nonsyn_nonhyper(i)-z*sd_nonsyn_nonhyper(i)
        neg_nonsyn_nonhyper=[neg_nonsyn_nonhyper,i];
    elseif P_nonsyn_c_nonhyper(i)>mean_nonsyn_nonhyper(i)+z*sd_nonsyn_nonhyper(i)
        pos_nonsyn_nonhyper=[pos_nonsyn_nonhyper,i];
    else
        inside_nonsyn_nonhyper=[inside_nonsyn_nonhyper,i];
    end
end



ind1=1:16; ind2=17:32; ind3=33:48; ind4=49:64; ind5=65:80; ind6=81:96;

figure
plot(ind1,mean_syn_hyper(ind1),'o',ind2,mean_syn_hyper(ind2),'o',ind3,mean_syn_hyper(ind3),'o',ind4,mean_syn_hyper(ind4),'o',ind5,mean_syn_hyper(ind5),'o',ind6,mean_syn_hyper(ind6),'o','LineWidth',1)
yline(0)
hold on
er = errorbar(1:96,mean_syn_hyper,z*sd_syn_hyper,z*sd_syn_hyper);
er.LineStyle = 'none';
er.LineWidth = 0.75;
plot(ind1,P_syn_c_hyper(ind1),'x',ind2,P_syn_c_hyper(ind2),'x',ind3,P_syn_c_hyper(ind3),'x',ind4,P_syn_c_hyper(ind4),'x',ind5,P_syn_c_hyper(ind5),'x',ind6,P_syn_c_hyper(ind6),'x','LineWidth',1)
title({"Number of syn mutations in 5% randomly-chosen HprM mutations [50,000 replicates]","Assumed normal distribution and used z-score (99.95% CI)"})
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])

figure
plot(ind1,mean_nonsyn_hyper(ind1),'o',ind2,mean_nonsyn_hyper(ind2),'o',ind3,mean_nonsyn_hyper(ind3),'o',ind4,mean_nonsyn_hyper(ind4),'o',ind5,mean_nonsyn_hyper(ind5),'o',ind6,mean_nonsyn_hyper(ind6),'o','LineWidth',1)
yline(0)
hold on
er = errorbar(1:96,mean_nonsyn_hyper,z*sd_nonsyn_hyper,z*sd_nonsyn_hyper);
er.LineStyle = 'none';
er.LineWidth = 0.75;
plot(ind1,P_nonsyn_c_hyper(ind1),'x',ind2,P_nonsyn_c_hyper(ind2),'x',ind3,P_nonsyn_c_hyper(ind3),'x',ind4,P_nonsyn_c_hyper(ind4),'x',ind5,P_nonsyn_c_hyper(ind5),'x',ind6,P_nonsyn_c_hyper(ind6),'x','LineWidth',1)
title({"Number of Nonsyn mutations in 5% randomly-chosen HprM mutations [50,000 replicates]","Assumed normal distribution and used z-score (99.95% CI)"})
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])



figure
plot(ind1,mean_syn_nonhyper(ind1),'o',ind2,mean_syn_nonhyper(ind2),'o',ind3,mean_syn_nonhyper(ind3),'o',ind4,mean_syn_nonhyper(ind4),'o',ind5,mean_syn_nonhyper(ind5),'o',ind6,mean_syn_nonhyper(ind6),'o','LineWidth',1)
yline(0)
hold on
er = errorbar(1:96,mean_syn_nonhyper,z*sd_syn_nonhyper,z*sd_syn_nonhyper);
er.LineStyle = 'none';
er.LineWidth = 0.75;
plot(ind1,P_syn_c_nonhyper(ind1),'x',ind2,P_syn_c_nonhyper(ind2),'x',ind3,P_syn_c_nonhyper(ind3),'x',ind4,P_syn_c_nonhyper(ind4),'x',ind5,P_syn_c_nonhyper(ind5),'x',ind6,P_syn_c_nonhyper(ind6),'x','LineWidth',1)
title({"Number of syn mutations in 5% randomly-chosen NHprM mutations [50,000 replicates]","Assumed normal distribution and used z-score (99.95% CI)"})
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])

figure
plot(ind1,mean_nonsyn_nonhyper(ind1),'o',ind2,mean_nonsyn_nonhyper(ind2),'o',ind3,mean_nonsyn_nonhyper(ind3),'o',ind4,mean_nonsyn_nonhyper(ind4),'o',ind5,mean_nonsyn_nonhyper(ind5),'o',ind6,mean_nonsyn_nonhyper(ind6),'o','LineWidth',1)
yline(0)
hold on
er = errorbar(1:96,mean_nonsyn_nonhyper,z*sd_nonsyn_nonhyper,z*sd_nonsyn_nonhyper);
er.LineStyle = 'none';
er.LineWidth = 0.75;
plot(ind1,P_nonsyn_c_nonhyper(ind1),'x',ind2,P_nonsyn_c_nonhyper(ind2),'x',ind3,P_nonsyn_c_nonhyper(ind3),'x',ind4,P_nonsyn_c_nonhyper(ind4),'x',ind5,P_nonsyn_c_nonhyper(ind5),'x',ind6,P_nonsyn_c_nonhyper(ind6),'x','LineWidth',1)
title({"Number of Nonsyn mutations in 5% randomly-chosen NHprM mutations [50,000 replicates]","Assumed normal distribution and used z-score (99.95% CI)"})
legend("C>A  Tv","C>G  Tv","C>T  Ts","T>A  Tv","T>C  Ts","T>G  Tv")
xticks(1:96)
xticklabels(["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT","ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"])




