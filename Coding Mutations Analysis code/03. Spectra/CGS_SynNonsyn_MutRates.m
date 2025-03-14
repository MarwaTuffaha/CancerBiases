 % Find Cancer-Gene Spectrum (CGS) syn/nonsyn
 % and NonCancer-Gene Spectrum (NCGS) syn/nonsyn
 % Substitution rate spectra (normalized by genome content)

Pc=load('Path_To\Data_Extracted_files\PCAWG_3mer_synNonsyn_Tissues_cancerGenes.mat');
P_syn_tissues_c=Pc.P_syn_c; % Nonsynonymous mutation count in cancer genes for 3mers in each tissue type
P_nonsyn_tissues_c=Pc.P_nonsyn_c; % Nonsynonymous mutation count in cancer genes for 3mers in each tissue type
P_syn_tissues_nc=Pc.P_syn_nc; % Nonsynonymous mutation count in cancer genes for 3mers in each tissue type
P_nonsyn_tissues_nc=Pc.P_nonsyn_nc; % Nonsynonymous mutation count in cancer genes for 3mers in each tissue type

tissue_types=["Brain","Blood","Bone","Brain","Breast","Cervix","Colorectal","Esophagus","Gall Bladder","Head & neck","Kidney","Liver","Lung","Mesenchymal","Ovary","Pancreas","Prostate","Skin","Stomach","Uterus"];

Ph=load('Path_To\Data_Extracted_files\PCAWG_3mer_Nonsyn_Tissues_HyperNonhyper_cancerGenes.mat');
P_syn_tissues_hyper_c=Ph.P_syn_c_hyper; 
P_syn_tissues_nonhyper_c=Ph.P_syn_c_nonhyper; 
P_nonsyn_tissues_hyper_c=Ph.P_nonsyn_c_hyper; 
P_nonsyn_tissues_nonhyper_c=Ph.P_nonsyn_c_nonhyper; 
P_syn_tissues_hyper_nc=Ph.P_syn_nc_hyper; 
P_syn_tissues_nonhyper_nc=Ph.P_syn_nc_nonhyper; 
P_nonsyn_tissues_hyper_nc=Ph.P_nonsyn_nc_hyper; 
P_nonsyn_tissues_nonhyper_nc=Ph.P_nonsyn_nc_nonhyper; 

Gc=load('Path_To\Data_Extracted_files\Germline_synNonsyn_CancerGenes.mat');
G_syn_c=Gc.G_syn_c; G_nonsyn_c=Gc.G_nonsyn_c;
G_syn_nc=Gc.G_syn_nc; G_nonsyn_nc=Gc.G_nonsyn_nc;

%% 
P_syn_c=sum(P_syn_tissues_c,2);
CGS_syn_all= P_syn_c./C_syn_c;
CGS_syn_tissues=P_syn_tissues_c./C_syn_c;
P_syn_nc=sum(P_syn_tissues_nc,2);
NCGS_syn_all= P_syn_nc./C_syn_nc;
NCGS_syn_tissues=P_syn_tissues_nc./C_syn_nc;

P_syn_hyper_c=sum(P_syn_tissues_hyper_c,2);
P_syn_nonhyper_c=sum(P_syn_tissues_nonhyper_c,2);
CGS_syn_hyper= P_syn_hyper_c./C_syn_c;
CGS_syn_nonhyper= P_syn_nonhyper_c./C_syn_c;
CGS_syn_tissues_hyper=P_syn_tissues_hyper_c./C_syn_c;
CGS_syn_tissues_nonhyper=P_syn_tissues_nonhyper_c./C_syn_c;

P_syn_hyper_nc=sum(P_syn_tissues_hyper_nc,2);
P_syn_nonhyper_nc=sum(P_syn_tissues_nonhyper_nc,2);
NCGS_syn_hyper= P_syn_hyper_nc./C_syn_nc;
NCGS_syn_nonhyper= P_syn_nonhyper_nc./C_syn_nc;
NCGS_syn_tissues_hyper=P_syn_tissues_hyper_nc./C_syn_nc;
NCGS_syn_tissues_nonhyper=P_syn_tissues_nonhyper_nc./C_syn_nc;

CGS_syn_Germline=G_syn_c./C_syn_c;
NCGS_syn_Germline=G_syn_nc./C_syn_nc;

P_nonsyn_c=sum(P_nonsyn_tissues_c,2);
CGS_nonsyn_all= P_nonsyn_c./C_nonsyn_c;
CGS_nonsyn_tissues=P_nonsyn_tissues_c./C_nonsyn_c;

P_nonsyn_nc=sum(P_nonsyn_tissues_nc,2);
NCGS_nonsyn_all= P_nonsyn_nc./C_nonsyn_nc;
NCGS_nonsyn_tissues=P_nonsyn_tissues_nc./C_nonsyn_nc;

P_nonsyn_hyper_c=sum(P_nonsyn_tissues_hyper_c,2);
P_nonsyn_nonhyper_c=sum(P_nonsyn_tissues_nonhyper_c,2);
CGS_nonsyn_hyper= P_nonsyn_hyper_c./C_nonsyn_c;
CGS_nonsyn_nonhyper= P_nonsyn_nonhyper_c./C_nonsyn_c;
CGS_nonsyn_tissues_hyper=P_nonsyn_tissues_hyper_c./C_nonsyn_c;
CGS_nonsyn_tissues_nonhyper=P_nonsyn_tissues_nonhyper_c./C_nonsyn_c;

P_nonsyn_hyper_nc=sum(P_nonsyn_tissues_hyper_nc,2);
P_nonsyn_nonhyper_nc=sum(P_nonsyn_tissues_nonhyper_nc,2);
NCGS_nonsyn_hyper= P_nonsyn_hyper_nc./C_nonsyn_nc;
NCGS_nonsyn_nonhyper= P_nonsyn_nonhyper_nc./C_nonsyn_nc;
NCGS_nonsyn_tissues_hyper=P_nonsyn_tissues_hyper_nc./C_nonsyn_nc;
NCGS_nonsyn_tissues_nonhyper=P_nonsyn_tissues_nonhyper_nc./C_nonsyn_nc;

CGS_nonsyn_Germline=G_nonsyn_c./C_nonsyn_c;
NCGS_nonsyn_Germline=G_nonsyn_nc./C_nonsyn_nc;


CGS_syn_tissues(isnan(CGS_syn_tissues))=0;
CGS_syn_tissues_hyper(isnan(CGS_syn_tissues_hyper))=0;
CGS_syn_tissues_nonhyper(isnan(CGS_syn_tissues_nonhyper))=0;

CGS_syn_all(isnan(CGS_syn_all))=0;
CGS_syn_hyper(isnan(CGS_syn_hyper))=0;
CGS_syn_nonhyper(isnan(CGS_syn_nonhyper))=0;
CGS_syn_Germline(isnan(CGS_syn_Germline))=0;

NCGS_syn_tissues(isnan(NCGS_syn_tissues))=0;
NCGS_syn_tissues_hyper(isnan(NCGS_syn_tissues_hyper))=0;
NCGS_syn_tissues_nonhyper(isnan(NCGS_syn_tissues_nonhyper))=0;

NCGS_syn_all(isnan(NCGS_syn_all))=0;
NCGS_syn_hyper(isnan(NCGS_syn_hyper))=0;
NCGS_syn_nonhyper(isnan(NCGS_syn_nonhyper))=0;
NCGS_syn_Germline(isnan(NCGS_syn_Germline))=0;


%Normalize
CGS_syn_all=CGS_syn_all/sum(CGS_syn_all);
CGS_syn_hyper=CGS_syn_hyper/sum(CGS_syn_hyper);
CGS_syn_nonhyper=CGS_syn_nonhyper/sum(CGS_syn_nonhyper);
CGS_syn_Germline=CGS_syn_Germline/sum(CGS_syn_Germline);

NCGS_syn_all=NCGS_syn_all/sum(NCGS_syn_all);
NCGS_syn_hyper=NCGS_syn_hyper/sum(NCGS_syn_hyper);
NCGS_syn_nonhyper=NCGS_syn_nonhyper/sum(NCGS_syn_nonhyper);
NCGS_syn_Germline=NCGS_syn_Germline/sum(NCGS_syn_Germline);

CGS_nonsyn_all=CGS_nonsyn_all/sum(CGS_nonsyn_all);
CGS_nonsyn_hyper=CGS_nonsyn_hyper/sum(CGS_nonsyn_hyper);
CGS_nonsyn_nonhyper=CGS_nonsyn_nonhyper/sum(CGS_nonsyn_nonhyper);
CGS_nonsyn_Germline=CGS_nonsyn_Germline/sum(CGS_nonsyn_Germline);

NCGS_nonsyn_all=NCGS_nonsyn_all/sum(NCGS_nonsyn_all);
NCGS_nonsyn_hyper=NCGS_nonsyn_hyper/sum(NCGS_nonsyn_hyper);
NCGS_nonsyn_nonhyper=NCGS_nonsyn_nonhyper/sum(NCGS_nonsyn_nonhyper);
NCGS_nonsyn_Germline=NCGS_nonsyn_Germline/sum(NCGS_nonsyn_Germline);

for i=1:20
    CGS_syn_tissues(:,i)=CGS_syn_tissues(:,i)/sum(CGS_syn_tissues(:,i));
    CGS_syn_tissues_hyper(:,i)=CGS_syn_tissues_hyper(:,i)/sum(CGS_syn_tissues_hyper(:,i));
    CGS_syn_tissues_nonhyper(:,i)=CGS_syn_tissues_nonhyper(:,i)/sum(CGS_syn_tissues_nonhyper(:,i));    
    CGS_nonsyn_tissues(:,i)=CGS_nonsyn_tissues(:,i)/sum(CGS_nonsyn_tissues(:,i));
    CGS_nonsyn_tissues_hyper(:,i)=CGS_nonsyn_tissues_hyper(:,i)/sum(CGS_nonsyn_tissues_hyper(:,i));
    CGS_nonsyn_tissues_nonhyper(:,i)=CGS_nonsyn_tissues_nonhyper(:,i)/sum(CGS_nonsyn_tissues_nonhyper(:,i));
end

for i=1:20
    NCGS_syn_tissues(:,i)=NCGS_syn_tissues(:,i)/sum(NCGS_syn_tissues(:,i));
    NCGS_syn_tissues_hyper(:,i)=NCGS_syn_tissues_hyper(:,i)/sum(NCGS_syn_tissues_hyper(:,i));
    NCGS_syn_tissues_nonhyper(:,i)=NCGS_syn_tissues_nonhyper(:,i)/sum(NCGS_syn_tissues_nonhyper(:,i));    
    NCGS_nonsyn_tissues(:,i)=NCGS_nonsyn_tissues(:,i)/sum(NCGS_nonsyn_tissues(:,i));
    NCGS_nonsyn_tissues_hyper(:,i)=NCGS_nonsyn_tissues_hyper(:,i)/sum(NCGS_nonsyn_tissues_hyper(:,i));
    NCGS_nonsyn_tissues_nonhyper(:,i)=NCGS_nonsyn_tissues_nonhyper(:,i)/sum(NCGS_nonsyn_tissues_nonhyper(:,i));
end





save("Path_To\Data_Extracted_files\CGS_NCGS_MutRates","CGS_nonsyn_all","CGS_nonsyn_tissues","CGS_nonsyn_hyper","CGS_nonsyn_tissues_hyper","CGS_nonsyn_nonhyper","CGS_nonsyn_tissues_nonhyper","CGS_nonsyn_Germline","CGS_syn_all","CGS_syn_tissues","CGS_syn_hyper","CGS_syn_tissues_hyper","CGS_syn_nonhyper","CGS_syn_tissues_nonhyper","CGS_syn_Germline","NCGS_nonsyn_all","NCGS_nonsyn_tissues","NCGS_nonsyn_hyper","NCGS_nonsyn_tissues_hyper","NCGS_nonsyn_nonhyper","NCGS_nonsyn_tissues_nonhyper","NCGS_nonsyn_Germline","NCGS_syn_all","NCGS_syn_tissues","NCGS_syn_hyper","NCGS_syn_tissues_hyper","NCGS_syn_nonhyper","NCGS_syn_tissues_nonhyper","NCGS_syn_Germline")

