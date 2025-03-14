 % Find 3mer spectra in all genes

P=load('Path_To\Data_Extracted_files\PCAWG_3mer_synNonsyn_Tissues.mat');
P_nonsyn_tissues=P.counts_nonsyn; % Nonsynonymous mutation count for 3mers in each tissue type
P_syn_tissues=P.counts_syn; % Synonymous mutation count for 3mers in each tissue type
P_tissues=P_syn_tissues+P_nonsyn_tissues;
P=sum(P_tissues,2);

tissue_types=["Brain","Blood","Bone","Brain","Breast","Cervix","Colorectal","Esophagus","Gall Bladder","Head & neck","Kidney","Liver","Lung","Mesenchymal","Ovary","Pancreas","Prostate","Skin","Stomach","Uterus"];

Ph=load('Path_To\Data_Extracted_files\PCAWG_3mer_synNonsyn_hyperNonhyper.mat');
P_syn_tissues_hyper=Ph.P_syn_hyper_tissues; 
P_syn_tissues_nonhyper=Ph.P_syn_nonhyper_tissues; 
P_nonsyn_tissues_hyper=Ph.P_nonsyn_hyper_tissues; 
P_nonsyn_tissues_nonhyper=Ph.P_nonsyn_nonhyper_tissues; 

G=load('Path_To\Data_Extracted_files\Germline_synNonsyn.mat');
G_syn=G.G_syn;
G_nonsyn=G.G_nonsyn;

%% 
P_syn=sum(P_syn_tissues,2);
spectrum_syn_all= P_syn;
spectrum_syn_tissues=P_syn_tissues;

P_nonsyn=sum(P_nonsyn_tissues,2);
spectrum_nonsyn_all= P_nonsyn;
spectrum_nonsyn_tissues=P_nonsyn_tissues;

P_syn_hyper=sum(P_syn_tissues_hyper,2);
P_syn_nonhyper=sum(P_syn_tissues_nonhyper,2);
spectrum_syn_hyper= P_syn_hyper;
spectrum_syn_nonhyper= P_syn_nonhyper;
spectrum_syn_tissues_hyper=P_syn_tissues_hyper;
spectrum_syn_tissues_nonhyper=P_syn_tissues_nonhyper;

P_nonsyn_hyper=sum(P_nonsyn_tissues_hyper,2);
P_nonsyn_nonhyper=sum(P_nonsyn_tissues_nonhyper,2);
spectrum_nonsyn_hyper= P_nonsyn_hyper;
spectrum_nonsyn_nonhyper= P_nonsyn_nonhyper;
spectrum_nonsyn_tissues_hyper=P_nonsyn_tissues_hyper;
spectrum_nonsyn_tissues_nonhyper=P_nonsyn_tissues_nonhyper;

P_hyper=P_syn_hyper+P_nonsyn_hyper;
P_nonhyper=P_syn_nonhyper+P_nonsyn_nonhyper;
P_tissues_hyper=P_syn_tissues_hyper+P_nonsyn_tissues_hyper;
P_tissues_nonhyper=P_syn_tissues_nonhyper+P_nonsyn_tissues_nonhyper;

spectrum_syn_Germline=G_syn;
spectrum_nonsyn_Germline=G_nonsyn;

% Normalize
spectrum_syn_all=spectrum_syn_all/sum(spectrum_syn_all);
spectrum_syn_hyper=spectrum_syn_hyper/sum(spectrum_syn_hyper);
spectrum_syn_nonhyper=spectrum_syn_nonhyper/sum(spectrum_syn_nonhyper);
spectrum_syn_Germline=spectrum_syn_Germline/sum(spectrum_syn_Germline);
for i=1:20
    spectrum_syn_tissues(:,i)=spectrum_syn_tissues(:,i)/sum(spectrum_syn_tissues(:,i));
    spectrum_syn_tissues_hyper(:,i)=spectrum_syn_tissues_hyper(:,i)/sum(spectrum_syn_tissues_hyper(:,i));
    spectrum_syn_tissues_nonhyper(:,i)=spectrum_syn_tissues_nonhyper(:,i)/sum(spectrum_syn_tissues_nonhyper(:,i));
end

spectrum_nonsyn_all=spectrum_nonsyn_all/sum(spectrum_nonsyn_all);
spectrum_nonsyn_hyper=spectrum_nonsyn_hyper/sum(spectrum_nonsyn_hyper);
spectrum_nonsyn_nonhyper=spectrum_nonsyn_nonhyper/sum(spectrum_nonsyn_nonhyper);
spectrum_nonsyn_Germline=spectrum_nonsyn_Germline/sum(spectrum_nonsyn_Germline);
for i=1:20
    spectrum_nonsyn_tissues(:,i)=spectrum_nonsyn_tissues(:,i)/sum(spectrum_nonsyn_tissues(:,i));
    spectrum_nonsyn_tissues_hyper(:,i)=spectrum_nonsyn_tissues_hyper(:,i)/sum(spectrum_nonsyn_tissues_hyper(:,i));
    spectrum_nonsyn_tissues_nonhyper(:,i)=spectrum_nonsyn_tissues_nonhyper(:,i)/sum(spectrum_nonsyn_tissues_nonhyper(:,i));
end

spectrum_tissues_hyper=P_syn_tissues_hyper+P_nonsyn_tissues_hyper;
spectrum_tissues_nonhyper=P_syn_tissues_nonhyper+P_nonsyn_tissues_nonhyper;
spectrum_tissues=spectrum_tissues_hyper+spectrum_tissues_nonhyper;
spectrum_hyper=P_syn_hyper+P_nonsyn_hyper;
spectrum_nonhyper=P_syn_nonhyper+P_nonsyn_nonhyper;
spectrum_all=spectrum_hyper+spectrum_nonhyper;
spectrum_Germline=G_syn+G_nonsyn;

spectrum_all=spectrum_all/sum(spectrum_all);
spectrum_hyper=spectrum_hyper/sum(spectrum_hyper);
spectrum_nonhyper=spectrum_nonhyper/sum(spectrum_nonhyper);
spectrum_Germline=spectrum_Germline/sum(spectrum_Germline);
for i=1:20
    spectrum_tissues(:,i)=spectrum_tissues(:,i)/sum(spectrum_tissues(:,i));
    spectrum_tissues_hyper(:,i)=spectrum_tissues_hyper(:,i)/sum(spectrum_tissues_hyper(:,i));
    spectrum_tissues_nonhyper(:,i)=spectrum_tissues_nonhyper(:,i)/sum(spectrum_tissues_nonhyper(:,i));
end


save("Path_To\Data_Extracted_files\Spectra_Counts","spectrum_nonsyn_all","spectrum_nonsyn_tissues","spectrum_nonsyn_hyper","spectrum_nonsyn_tissues_hyper","spectrum_nonsyn_nonhyper","spectrum_nonsyn_tissues_nonhyper","spectrum_nonsyn_Germline","spectrum_all","spectrum_tissues","spectrum_hyper","spectrum_tissues_hyper","spectrum_nonhyper","spectrum_tissues_nonhyper","spectrum_Germline","spectrum_syn_all","spectrum_syn_tissues","spectrum_syn_hyper","spectrum_syn_tissues_hyper","spectrum_syn_nonhyper","spectrum_syn_tissues_nonhyper","spectrum_syn_Germline")


%% Now find spectra of all, hyper, nonhyper samples except the indexed tissue
spectrum_nonhyper_ExceptTissue=zeros(96,20);
spectrum_hyper_ExceptTissue=zeros(96,20);
spectrum_ExceptTissue=zeros(96,20);

P_ExceptTissue=zeros(96,20);
P_nonhyper_ExceptTissue=zeros(96,20);
P_hyper_ExceptTissue=zeros(96,20);

for i=1:20
    P_ExceptTissue(:,i)=P-P_tissues(:,i);
    spectrum_ExceptTissue(:,i)=P_ExceptTissue(:,i);
    P_nonhyper_ExceptTissue(:,i)=P_nonhyper-P_tissues_nonhyper(:,i);
    spectrum_nonhyper_ExceptTissue(:,i)=P_nonhyper_ExceptTissue(:,i);
    P_hyper_ExceptTissue(:,i)=P_hyper-P_tissues_nonhyper(:,i);
    spectrum_hyper_ExceptTissue(:,i)=P_hyper_ExceptTissue(:,i);
end

% Normalize
for i=1:20
    spectrum_ExceptTissue(:,i)=spectrum_ExceptTissue(:,i)/sum(spectrum_ExceptTissue(:,i));
    spectrum_hyper_ExceptTissue(:,i)=spectrum_hyper_ExceptTissue(:,i)/sum(spectrum_hyper_ExceptTissue(:,i));
    spectrum_nonhyper_ExceptTissue(:,i)=spectrum_nonhyper_ExceptTissue(:,i)/sum(spectrum_nonhyper_ExceptTissue(:,i));
end

save("Path_To\Data_Extracted_files\Spectra_ExceptTissue_Counts","spectrum_ExceptTissue","spectrum_hyper_ExceptTissue","spectrum_nonhyper_ExceptTissue")
