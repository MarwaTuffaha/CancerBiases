% Correlate spectra of different tissues

load("Path_To\Data_Extracted_files\Spectra_MutRates.mat","spectrum_all","spectrum_tissues","spectrum_hyper","spectrum_tissues_hyper","spectrum_nonhyper","spectrum_tissues_nonhyper","spectrum_Germline")
load("Path_To\Data_Extracted_files\Spectra_ExceptTissue_MutRates.mat","spectrum_ExceptTissue","spectrum_hyper_ExceptTissue","spectrum_nonhyper_ExceptTissue")

tissue_types=["Brain","Blood","Bone","Brain","Breast","Cervix","Colorectal","Esophagus","Gall Bladder","Head & neck","Kidney","Liver","Lung","Mesenchymal","Ovary","Pancreas","Prostate","Skin","Stomach","Uterus"];

figure
a=[spectrum_hyper,spectrum_nonhyper];
plotting_s(a,"Hyper","Nonhyper")

%% Germline

figure
subplot(1,2,1)
a=[spectrum_hyper,spectrum_Germline];
plotting_s(a,"Hyper","Germline")
subplot(1,2,2)
a=[spectrum_nonhyper,spectrum_Germline];
plotting_s(a,"Nonhyper","Germline")
sgtitle("Mutation rate spectra")

corrCoeff_Germline=zeros(20,1);
figure
for i=1:20
    subplot(4,5,i)
    a=[spectrum_tissues(:,i),spectrum_Germline];
    R= corr(a,'Type','Spearman');
    corrCoeff_Germline(i)=R(1,2);
    plotting_s(a,tissue_types(i),"Germline")
    set(gca,'xtick',[])
end
sgtitle("All samples - Mutation rate Spectra")
close

corrCoeff_Germline_hyper=zeros(20,1);
figure
for i=1:20
    subplot(4,5,i)
    a=[spectrum_tissues_hyper(:,i),spectrum_Germline];
    R= corr(a,'Type','Spearman');
    corrCoeff_Germline_hyper(i)=R(1,2);
    plotting_s(a,tissue_types(i),"Germline")
    set(gca,'xtick',[])
end
sgtitle("Hyper-mutated samples - Mutation rate Spectra")
close

corrCoeff_Germline_nonhyper=zeros(20,1);
figure
for i=1:20
    subplot(4,5,i)
    a=[spectrum_tissues_nonhyper(:,i),spectrum_Germline];
    R= corr(a,'Type','Spearman');
    corrCoeff_Germline_nonhyper(i)=R(1,2);
    plotting_s(a,tissue_types(i),"Germline")
    set(gca,'xtick',[])
end
sgtitle("Nonhyper-mutated samples - Mutation rate Spectra")
close

figure
text(1*ones(20,1),corrCoeff_Germline,tissue_types)
text(2*ones(20,1),corrCoeff_Germline_hyper,tissue_types)
text(3*ones(20,1),corrCoeff_Germline_nonhyper,tissue_types)
ylabel("Spearman Correlation Coeff.")
xlim([0.5 3.5])
xticks([ 1 2 3])
xticklabels({ sprintf('All samples'), sprintf('Hyper-mutated Samples'), sprintf('Nonhyper-mutated samples')})
title("Correlating spectra of different tissues with the Germline Spectrum - Mutation rate Spectra")

%%

corrCoeff_all=zeros(20,1);
figure
for i=1:20
    subplot(4,5,i)
    a=[spectrum_tissues(:,i),spectrum_ExceptTissue(:,i)];
    R= corr(a,'Type','Spearman');
    corrCoeff_all(i)=R(1,2);
    plotting_s(a,tissue_types(i),"All Samples")
    set(gca,'xtick',[])
end
sgtitle("Mutation rate Spectra")
close


corrCoeff_nonhyper=zeros(20,1);
figure
for i=1:20
    subplot(4,5,i)
    a=[spectrum_tissues_nonhyper(:,i),spectrum_nonhyper_ExceptTissue(:,i)];
    R= corr(a,'Type','Spearman');
    corrCoeff_nonhyper(i)=R(1,2);
    plotting_s(a,tissue_types(i),"Nonhyper Samples")
    set(gca,'xtick',[])
end
sgtitle("Mutation rate Spectra")
close

corrCoeff_hyper=zeros(20,1);
figure
for i=1:20
    subplot(4,5,i)
    a=[spectrum_tissues_hyper(:,i),spectrum_hyper_ExceptTissue(:,i)];
    R= corr(a,'Type','Spearman');
    corrCoeff_hyper(i)=R(1,2);
    plotting_s(a,tissue_types(i),"Hyper Samples")
    set(gca,'xtick',[])
end
sgtitle("Mutation rate Spectra")
close


figure
text(1*ones(20,1),corrCoeff_all,tissue_types)
text(2*ones(20,1),corrCoeff_hyper,tissue_types)
text(3*ones(20,1),corrCoeff_nonhyper,tissue_types)
ylabel("Spearman Correlation Coeff.")
xlim([0.5 3.5])
xticks([ 1 2 3])
xticklabels({ sprintf('All samples'), sprintf('Hyper-mutated Samples'), sprintf('Nonhyper-mutated samples')})
title("Correlating spectra of different tissues with Other Tissues spectra - Mutation rate Spectra")






