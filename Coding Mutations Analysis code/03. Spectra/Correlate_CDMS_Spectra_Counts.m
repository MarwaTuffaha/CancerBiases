% Correlate CDMS and spectra

load("Path_To\Data_Extracted_files\CDMS_Counts","CDMS_all","CDMS_tissues","CDMS_hyper","CDMS_tissues_hyper","CDMS_nonhyper","CDMS_tissues_nonhyper","CDMS_Germline")
load("Path_To\Data_Extracted_files\Spectra_Counts","spectrum_all","spectrum_tissues","spectrum_hyper","spectrum_tissues_hyper","spectrum_nonhyper","spectrum_tissues_nonhyper","spectrum_Germline")

tissue_types=["Bladdar","Blood","Bone","Brain","Breast","Cervix","Colorectal","Esophagus","Gall Bladder","Head & neck","Kidney","Liver","Lung","Mesenchymal","Ovary","Pancreas","Prostate","Skin","Stomach","Uterus"];

counts=load("Path_To\Data_Extracted_files\Count_3mer_mutations_tissues.mat");
counts_all_nonsyn=counts.count_all_tissues_nonsyn;
counts_CancerGene_nonsyn=counts.count_CancerGene_tissues_nonsyn;

figure
a=[spectrum_Germline,CDMS_Germline];
plotting_s(a,"Germline","CDMS (Germline)")


corrCoeff=zeros(20,1);
figure
for i=1:20
    subplot(4,5,i)
    a=[spectrum_tissues(:,i),CDMS_tissues(:,i)];
    R= corr(a,'Type','Spearman');
    corrCoeff(i)=R(1,2);
    plotting_s(a,tissue_types(i),strcat("CDMS ",num2str(i)))
    set(gca,'xtick',[])
end
sgtitle("All samples - Counts spectra")

figure
a=[counts_CancerGene_nonsyn',corrCoeff];
plotting_p(a,"# Cancer-genes nonsyn mutations","Pearson CorrCoeff(CDMS&AGS)/tissue")


%
figure
subplot(1,3,1)
a=[spectrum_all,CDMS_all];
plotting_s(a,"All","CDMS (All)")

subplot(1,3,2)
a=[spectrum_hyper,CDMS_hyper];
plotting_s(a,"Hyper","CDMS (Hyper)")

subplot(1,3,3)
a=[spectrum_nonhyper,CDMS_nonhyper];
plotting_s(a,"Nonhyper","CDMS (Nonhyper)")
sgtitle("Counts spectra")



load('Path_To\Data_Extracted_files\codon_3mer_CancerGenes_synNonsyn.mat','C_nonsyn_c','C_syn_c','C_nonsyn_nc','C_syn_nc')


figure
subplot(2,1,1)
a=[C_syn_c+C_syn_nc,C_syn_c];
plotting_s(a,"AGS Syn (Unbiased)","CGS Syn (Unbiased)")

subplot(2,1,2)
a=[C_nonsyn_c+C_nonsyn_nc,C_nonsyn_c];
plotting_s(a,"AGS NS (Unbiased)","CGS NS (Unbiased)")



