% Correlate CDMS and spectra

load("Path_To\Data_Extracted_files\CDMS_MutRates","CDMS_all","CDMS_tissues","CDMS_hyper","CDMS_tissues_hyper","CDMS_nonhyper","CDMS_tissues_nonhyper","CDMS_Germline")
load("Path_To\Data_Extracted_files\Spectra_MutRates","spectrum_all","spectrum_tissues","spectrum_hyper","spectrum_tissues_hyper","spectrum_nonhyper","spectrum_tissues_nonhyper","spectrum_Germline")

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
    R= corrcoef(a);
    corrCoeff(i)=R(1,2);
    plotting_s(a,tissue_types(i),strcat("CDMS ",num2str(i)))
    set(gca,'xtick',[])
end
sgtitle("All samples - Mutation rate Spectra")

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
sgtitle("Mutation rate Spectra")







