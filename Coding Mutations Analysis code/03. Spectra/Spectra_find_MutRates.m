 % Find 3mer spectra in all genes

P=load('Path_To\Data_Extracted_files\PCAWG_3mer_synNonsyn_Tissues.mat');
P_nonsyn_tissues=P.counts_nonsyn; % Nonsynonymous mutation count for 3mers in each tissue type
P_syn_tissues=P.counts_syn; % Synonymous mutation count for 3mers in each tissue type

tissue_types=["Brain","Blood","Bone","Brain","Breast","Cervix","Colorectal","Esophagus","Gall Bladder","Head & neck","Kidney","Liver","Lung","Mesenchymal","Ovary","Pancreas","Prostate","Skin","Stomach","Uterus"];

C=load('Path_To\Data_Extracted_files\codon_3mer_synNonsyn.mat');
C_nonsyn=C.C_nonsyn; % Nonsynonymous mutation count in genome for 3mers
C_syn=C.C_syn;

Ph=load('Path_To\Data_Extracted_files\PCAWG_3mer_synNonsyn_hyperNonhyper.mat');
P_syn_tissues_hyper=Ph.P_syn_hyper_tissues; 
P_syn_tissues_nonhyper=Ph.P_syn_nonhyper_tissues;
P_nonsyn_tissues_hyper=Ph.P_nonsyn_hyper_tissues; 
P_nonsyn_tissues_nonhyper=Ph.P_nonsyn_nonhyper_tissues; 

G=load('Path_To\Data_Extracted_files\Germline_synNonsyn.mat');
G_syn=G.G_syn;
G_nonsyn=G.G_nonsyn;

noSyn_ind=[4 20 52 84];

%% 

P_syn=sum(P_syn_tissues,2);
spectrum_syn_all= P_syn./C_syn;
spectrum_syn_tissues=P_syn_tissues./C_syn;

P_syn_hyper=sum(P_syn_tissues_hyper,2);
P_syn_nonhyper=sum(P_syn_tissues_nonhyper,2);
spectrum_syn_hyper= P_syn_hyper./C_syn;
spectrum_syn_nonhyper= P_syn_nonhyper./C_syn;
spectrum_syn_tissues_hyper=P_syn_tissues_hyper./C_syn;
spectrum_syn_tissues_nonhyper=P_syn_tissues_nonhyper./C_syn;

spectrum_syn_Germline=G_syn./C_syn;


% remove nan values from syn spectra caused by lack of such syn mutations
for i=4:-1:1
    spectrum_syn_all(noSyn_ind(i))=[];
    spectrum_syn_hyper(noSyn_ind(i))=[];
    spectrum_syn_nonhyper(noSyn_ind(i))=[];
    spectrum_syn_Germline(noSyn_ind(i))=[];
    spectrum_syn_tissues(noSyn_ind(i),:)=[];
    spectrum_syn_tissues_hyper(noSyn_ind(i),:)=[];
    spectrum_syn_tissues_nonhyper(noSyn_ind(i),:)=[];
end

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

P_nonsyn=sum(P_nonsyn_tissues,2);
spectrum_nonsyn_all= P_nonsyn./C_nonsyn;
spectrum_nonsyn_tissues=P_nonsyn_tissues./C_nonsyn;

P_nonsyn_hyper=sum(P_nonsyn_tissues_hyper,2);
P_nonsyn_nonhyper=sum(P_nonsyn_tissues_nonhyper,2);
spectrum_nonsyn_hyper= P_nonsyn_hyper./C_nonsyn;
spectrum_nonsyn_nonhyper= P_nonsyn_nonhyper./C_nonsyn;
spectrum_nonsyn_tissues_hyper=P_nonsyn_tissues_hyper./C_nonsyn;
spectrum_nonsyn_tissues_nonhyper=P_nonsyn_tissues_nonhyper./C_nonsyn;

spectrum_nonsyn_Germline=G_nonsyn./C_nonsyn;

% Normalize
spectrum_nonsyn_all=spectrum_nonsyn_all/sum(spectrum_nonsyn_all);
spectrum_nonsyn_hyper=spectrum_nonsyn_hyper/sum(spectrum_nonsyn_hyper);
spectrum_nonsyn_nonhyper=spectrum_nonsyn_nonhyper/sum(spectrum_nonsyn_nonhyper);
spectrum_nonsyn_Germline=spectrum_nonsyn_Germline/sum(spectrum_nonsyn_Germline);
for i=1:20
    spectrum_nonsyn_tissues(:,i)=spectrum_nonsyn_tissues(:,i)/sum(spectrum_nonsyn_tissues(:,i));
    spectrum_nonsyn_tissues_hyper(:,i)=spectrum_nonsyn_tissues_hyper(:,i)/sum(spectrum_nonsyn_tissues_hyper(:,i));
    spectrum_nonsyn_tissues_nonhyper(:,i)=spectrum_nonsyn_tissues_nonhyper(:,i)/sum(spectrum_nonsyn_tissues_nonhyper(:,i));
end

C_coding=C_syn+C_nonsyn;
G_coding=G_syn+G_nonsyn;
P_coding_tissues=P_syn_tissues+P_nonsyn_tissues;
P_coding_tissues_hyper=P_syn_tissues_hyper+P_nonsyn_tissues_hyper;
P_coding_tissues_nonhyper=P_syn_tissues_nonhyper+P_nonsyn_tissues_nonhyper;

P_coding=P_syn+P_nonsyn;
spectrum_coding_all= P_coding./C_coding;
spectrum_coding_tissues=P_coding_tissues./C_coding;


P_coding_hyper=P_syn_hyper+P_nonsyn_hyper;
P_coding_nonhyper=P_syn_nonhyper+P_nonsyn_nonhyper;
spectrum_coding_hyper= P_coding_hyper./C_coding;
spectrum_coding_nonhyper= P_coding_nonhyper./C_coding;
spectrum_coding_tissues_hyper=P_coding_tissues_hyper./C_coding;
spectrum_coding_tissues_nonhyper=P_coding_tissues_nonhyper./C_coding;

spectrum_coding_Germline=G_coding./C_coding;

% Normalize
spectrum_coding_all=spectrum_coding_all/sum(spectrum_coding_all);
spectrum_coding_hyper=spectrum_coding_hyper/sum(spectrum_coding_hyper);
spectrum_coding_nonhyper=spectrum_coding_nonhyper/sum(spectrum_coding_nonhyper);
spectrum_coding_Germline=spectrum_coding_Germline/sum(spectrum_coding_Germline);
for i=1:20
    spectrum_coding_tissues(:,i)=spectrum_coding_tissues(:,i)/sum(spectrum_coding_tissues(:,i));
    spectrum_coding_tissues_hyper(:,i)=spectrum_coding_tissues_hyper(:,i)/sum(spectrum_coding_tissues_hyper(:,i));
    spectrum_coding_tissues_nonhyper(:,i)=spectrum_coding_tissues_nonhyper(:,i)/sum(spectrum_coding_tissues_nonhyper(:,i));
end



save("Path_To\Data_Extracted_files\Spectra_MutRates","spectrum_syn_all","spectrum_syn_tissues","spectrum_syn_hyper","spectrum_syn_tissues_hyper","spectrum_syn_nonhyper","spectrum_syn_tissues_nonhyper","spectrum_syn_Germline","spectrum_nonsyn_all","spectrum_nonsyn_tissues","spectrum_nonsyn_hyper","spectrum_nonsyn_tissues_hyper","spectrum_nonsyn_nonhyper","spectrum_nonsyn_tissues_nonhyper","spectrum_nonsyn_Germline","spectrum_coding_all","spectrum_coding_tissues","spectrum_coding_hyper","spectrum_coding_tissues_hyper","spectrum_coding_nonhyper","spectrum_coding_tissues_nonhyper","spectrum_coding_Germline")


%% Now find spectra of all, hyper, nonhyper samples except the indexed tissue

spectrum_syn_nonhyper_ExceptTissue=zeros(96,20);
spectrum_syn_hyper_ExceptTissue=zeros(96,20);
spectrum_syn_ExceptTissue=zeros(96,20);

P_syn_ExceptTissue=zeros(96,20);
P_syn_nonhyper_ExceptTissue=zeros(96,20);
P_syn_hyper_ExceptTissue=zeros(96,20);

for i=1:20
    P_syn_ExceptTissue(:,i)=P_syn-P_syn_tissues(:,i);
    spectrum_syn_ExceptTissue(:,i)=P_syn_ExceptTissue(:,i)./C_syn;
    P_syn_nonhyper_ExceptTissue(:,i)=P_syn_nonhyper-P_syn_tissues_nonhyper(:,i);
    spectrum_syn_nonhyper_ExceptTissue(:,i)=P_syn_nonhyper_ExceptTissue(:,i)./C_syn;
    P_syn_hyper_ExceptTissue(:,i)=P_syn_hyper-P_syn_tissues_nonhyper(:,i);
    spectrum_syn_hyper_ExceptTissue(:,i)=P_syn_hyper_ExceptTissue(:,i)./C_syn;
end

for i=4:-1:1
spectrum_syn_ExceptTissue(noSyn_ind(i),:)=[];
spectrum_syn_hyper_ExceptTissue(noSyn_ind(i),:)=[];
spectrum_syn_nonhyper_ExceptTissue(noSyn_ind(i),:)=[];
end

% Normalize
for i=1:20
    spectrum_syn_ExceptTissue(:,i)=spectrum_syn_ExceptTissue(:,i)/sum(spectrum_syn_ExceptTissue(:,i));
    spectrum_syn_hyper_ExceptTissue(:,i)=spectrum_syn_hyper_ExceptTissue(:,i)/sum(spectrum_syn_hyper_ExceptTissue(:,i));
    spectrum_syn_nonhyper_ExceptTissue(:,i)=spectrum_syn_nonhyper_ExceptTissue(:,i)/sum(spectrum_syn_nonhyper_ExceptTissue(:,i));
end

spectrum_nonsyn_nonhyper_ExceptTissue=zeros(96,20);
spectrum_nonsyn_hyper_ExceptTissue=zeros(96,20);
spectrum_nonsyn_ExceptTissue=zeros(96,20);

P_nonsyn_ExceptTissue=zeros(96,20);
P_nonsyn_nonhyper_ExceptTissue=zeros(96,20);
P_nonsyn_hyper_ExceptTissue=zeros(96,20);

for i=1:20
    P_nonsyn_ExceptTissue(:,i)=P_nonsyn-P_nonsyn_tissues(:,i);
    spectrum_nonsyn_ExceptTissue(:,i)=P_nonsyn_ExceptTissue(:,i)./C_nonsyn;
    P_nonsyn_nonhyper_ExceptTissue(:,i)=P_nonsyn_nonhyper-P_nonsyn_tissues_nonhyper(:,i);
    spectrum_nonsyn_nonhyper_ExceptTissue(:,i)=P_nonsyn_nonhyper_ExceptTissue(:,i)./C_nonsyn;
    P_nonsyn_hyper_ExceptTissue(:,i)=P_nonsyn_hyper-P_nonsyn_tissues_nonhyper(:,i);
    spectrum_nonsyn_hyper_ExceptTissue(:,i)=P_nonsyn_hyper_ExceptTissue(:,i)./C_nonsyn;
end

% Normalize
for i=1:20
    spectrum_nonsyn_ExceptTissue(:,i)=spectrum_nonsyn_ExceptTissue(:,i)/sum(spectrum_nonsyn_ExceptTissue(:,i));
    spectrum_nonsyn_hyper_ExceptTissue(:,i)=spectrum_nonsyn_hyper_ExceptTissue(:,i)/sum(spectrum_nonsyn_hyper_ExceptTissue(:,i));
    spectrum_nonsyn_nonhyper_ExceptTissue(:,i)=spectrum_nonsyn_nonhyper_ExceptTissue(:,i)/sum(spectrum_nonsyn_nonhyper_ExceptTissue(:,i));
end

spectrum_coding_nonhyper_ExceptTissue=zeros(96,20);
spectrum_coding_hyper_ExceptTissue=zeros(96,20);
spectrum_coding_ExceptTissue=zeros(96,20);

P_coding_ExceptTissue=P_syn_ExceptTissue+P_nonsyn_ExceptTissue;
P_coding_nonhyper_ExceptTissue=P_syn_nonhyper_ExceptTissue+P_nonsyn_nonhyper_ExceptTissue;
P_coding_hyper_ExceptTissue=P_syn_hyper_ExceptTissue+P_nonsyn_hyper_ExceptTissue;

for i=1:20
    spectrum_coding_ExceptTissue(:,i)=P_coding_ExceptTissue(:,i)./C_coding;
    spectrum_coding_nonhyper_ExceptTissue(:,i)=P_coding_nonhyper_ExceptTissue(:,i)./C_coding;
    spectrum_coding_hyper_ExceptTissue(:,i)=P_coding_hyper_ExceptTissue(:,i)./C_coding;
end

% Normalize
for i=1:20
    spectrum_coding_ExceptTissue(:,i)=spectrum_coding_ExceptTissue(:,i)/sum(spectrum_coding_ExceptTissue(:,i));
    spectrum_coding_hyper_ExceptTissue(:,i)=spectrum_coding_hyper_ExceptTissue(:,i)/sum(spectrum_coding_hyper_ExceptTissue(:,i));
    spectrum_coding_nonhyper_ExceptTissue(:,i)=spectrum_coding_nonhyper_ExceptTissue(:,i)/sum(spectrum_coding_nonhyper_ExceptTissue(:,i));
end

save("Path_To\Data_Extracted_files\Spectra_ExceptTissue_MutRates","spectrum_nonsyn_ExceptTissue","spectrum_nonsyn_hyper_ExceptTissue","spectrum_nonsyn_nonhyper_ExceptTissue","spectrum_syn_ExceptTissue","spectrum_syn_hyper_ExceptTissue","spectrum_syn_nonhyper_ExceptTissue","spectrum_coding_ExceptTissue","spectrum_coding_hyper_ExceptTissue","spectrum_coding_nonhyper_ExceptTissue")
