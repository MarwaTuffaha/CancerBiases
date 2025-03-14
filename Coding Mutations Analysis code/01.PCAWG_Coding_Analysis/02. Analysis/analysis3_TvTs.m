P=load("Path_To\Data_Extracted_files\PCAWG_TvTs_hyperNonhyper");

P_tv_nonsyn_hyper_tissues=P.P_tv_nonsyn_hyper_tissues;
P_ts_nonsyn_hyper_tissues=P.P_ts_nonsyn_hyper_tissues;
P_tv_nonsyn_nonhyper_tissues=P.P_tv_nonsyn_nonhyper_tissues;
P_ts_nonsyn_nonhyper_tissues=P.P_ts_nonsyn_nonhyper_tissues;
P_tv_syn_hyper_tissues=P.P_tv_syn_hyper_tissues;
P_ts_syn_hyper_tissues=P.P_ts_syn_hyper_tissues;
P_tv_syn_nonhyper_tissues=P.P_tv_syn_nonhyper_tissues;
P_ts_syn_nonhyper_tissues=P.P_ts_syn_nonhyper_tissues;

P_tv_nonsyn_tissues=P_tv_nonsyn_hyper_tissues+P_tv_nonsyn_nonhyper_tissues;
P_ts_nonsyn_tissues=P_ts_nonsyn_hyper_tissues+P_ts_nonsyn_nonhyper_tissues;
P_tv_syn_tissues=P_tv_syn_hyper_tissues+P_tv_syn_nonhyper_tissues;
P_ts_syn_tissues=P_ts_syn_hyper_tissues+P_ts_syn_nonhyper_tissues;

P_tv_hyper_tissues=P_tv_syn_hyper_tissues+P_tv_nonsyn_hyper_tissues;
P_ts_hyper_tissues=P_ts_syn_hyper_tissues+P_ts_nonsyn_hyper_tissues;
P_tv_nonhyper_tissues=P_tv_syn_nonhyper_tissues+P_tv_nonsyn_nonhyper_tissues;
P_ts_nonhyper_tissues=P_ts_syn_nonhyper_tissues+P_ts_nonsyn_nonhyper_tissues;

P_tv_all_tissues=P_tv_hyper_tissues+P_tv_nonhyper_tissues;
P_ts_all_tissues=P_ts_hyper_tissues+P_ts_nonhyper_tissues;

tissue_types=["Bladder","Blood","Bone","Brain","Breast","Cervix","Colorectal","Esophagus","Gall Bladder","Head & neck","Kidney","Liver","Lung","Mesenchymal","Ovary","Pancreas","Prostate","Skin","Stomach","Uterus"];

nonsyn=P_tv_nonsyn_tissues+P_ts_nonsyn_tissues;
syn=P_tv_syn_tissues+P_ts_syn_tissues;
syn_frac=syn./(nonsyn+syn);


X = categorical(tissue_types);
X = reordercats(X,tissue_types);

figure
bar(X,[nonsyn;syn],'stacked')
plots=get(gca, 'Children');
legend([plots(1, 1),plots(2, 1)], { 'Synonymous','Nonsynonymous'})
title("Number of mutations - coding PCAWG")


freq_germ_ult = load('D:\Marwa\1Western\PhD\Research\Thesis\Analyzing Data\Germline Dataset\Freq_germ_ult.mat');
freq_germ_ult = freq_germ_ult.Freq;

tv_hyper_tissues=P_tv_hyper_tissues./(P_tv_hyper_tissues+P_ts_hyper_tissues);
ts_hyper_tissues=P_ts_hyper_tissues./(P_tv_hyper_tissues+P_ts_hyper_tissues);
tv_nonhyper_tissues=P_tv_nonhyper_tissues./(P_tv_nonhyper_tissues+P_ts_nonhyper_tissues);
ts_nonhyper_tissues=P_ts_nonhyper_tissues./(P_tv_nonhyper_tissues+P_ts_nonhyper_tissues);

tv_all_tissues=P_tv_all_tissues./(P_tv_all_tissues+P_ts_all_tissues);
ts_all_tissues=P_ts_all_tissues./(P_tv_all_tissues+P_ts_all_tissues);


tv_nonsyn_hyper_tissues=P_tv_nonsyn_hyper_tissues./(P_tv_nonsyn_hyper_tissues+P_ts_nonsyn_hyper_tissues);
ts_nonsyn_hyper_tissues=P_ts_nonsyn_hyper_tissues./(P_tv_nonsyn_hyper_tissues+P_ts_nonsyn_hyper_tissues);
tv_nonsyn_nonhyper_tissues=P_tv_nonsyn_nonhyper_tissues./(P_tv_nonsyn_nonhyper_tissues+P_ts_nonsyn_nonhyper_tissues);
ts_nonsyn_nonhyper_tissues=P_ts_nonsyn_nonhyper_tissues./(P_tv_nonsyn_nonhyper_tissues+P_ts_nonsyn_nonhyper_tissues);
tv_syn_hyper_tissues=P_tv_syn_hyper_tissues./(P_tv_syn_hyper_tissues+P_ts_syn_hyper_tissues);
ts_syn_hyper_tissues=P_ts_syn_hyper_tissues./(P_tv_syn_hyper_tissues+P_ts_syn_hyper_tissues);
tv_syn_nonhyper_tissues=P_tv_syn_nonhyper_tissues./(P_tv_syn_nonhyper_tissues+P_ts_syn_nonhyper_tissues);
ts_syn_nonhyper_tissues=P_ts_syn_nonhyper_tissues./(P_tv_syn_nonhyper_tissues+P_ts_syn_nonhyper_tissues);

ts_syn_tissues=P_ts_syn_tissues./(P_tv_syn_tissues+P_ts_syn_tissues);
ts_nonsyn_tissues=P_ts_nonsyn_tissues./(P_tv_nonsyn_tissues+P_ts_nonsyn_tissues);


save("Path_To\Data_Extracted_files\TsTv_coding_result","ts_all_tissues","ts_hyper_tissues","ts_nonhyper_tissues","ts_nonsyn_tissues","ts_nonsyn_hyper_tissues","ts_nonsyn_nonhyper_tissues","ts_syn_tissues","ts_syn_hyper_tissues","ts_syn_nonhyper_tissues")

figure

text(1*ones(20,1),ts_all_tissues,tissue_types)
text(2*ones(20,1),ts_hyper_tissues,tissue_types)
text(3*ones(20,1),ts_nonhyper_tissues,tissue_types)

ylabel("Transitions Frequency")
xlim([ 0 4])
xticks([ 1 2 3])
xticklabels({ 'All', 'Hyper', 'Nonhyper'})
hw=yline(1-freq_germ_ult(1),'-m',"Germline",'LineWidth',2);
hw.FontSize = 12;
hy=yline(1/3,'-b',"Neutral",'LineWidth',2);
hy.FontSize = 12;
title("Coding Mutations Analysis - PCAWG")

ylim([0 1])



figure

text(1*ones(20,1),ts_nonsyn_tissues,tissue_types)
text(2*ones(20,1),ts_nonsyn_hyper_tissues,tissue_types)
text(3*ones(20,1),ts_nonsyn_nonhyper_tissues,tissue_types)

ylabel("Transitions Frequency")
xlim([ 0 4])
xticks([ 1 2 3])
xticklabels({ 'All', 'Hyper', 'Nonhyper'})
hw=yline(1-freq_germ_ult(1),'-m',"Germline",'LineWidth',2);
hw.FontSize = 12;
hy=yline(1/3,'-b',"Neutral",'LineWidth',2);
hy.FontSize = 12;
title("Nonsynonymous Mutations Analysis - PCAWG")

ylim([0 1])


figure

text(1*ones(20,1),ts_syn_tissues,tissue_types)
text(2*ones(20,1),ts_syn_hyper_tissues,tissue_types)
text(3*ones(20,1),ts_syn_nonhyper_tissues,tissue_types)

ylabel("Transitions Frequency")
xlim([ 0 4])
xticks([ 1 2 3])
xticklabels({ 'All', 'Hyper', 'Nonhyper'})
hw=yline(1-freq_germ_ult(1),'-m',"Germline",'LineWidth',2);
hw.FontSize = 12;
hy=yline(1/3,'-b',"Neutral",'LineWidth',2);
hy.FontSize = 12;
title("Synonymous Mutations Analysis - PCAWG")

ylim([0 1])




