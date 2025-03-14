% Example: Find bias reversal measures on the Ts:Tv, 1mer, and 3mer level
% for PCAWG data

Pc=load('Path_To\Data_Extracted_files\PCAWG_3mer_Nonsyn_Tissues_HyperNonhyper_cancerGenes.mat');
P_nonhyper_tissues_c=Pc.P_nonsyn_c_nonhyper;
P_hyper_tissues_c=Pc.P_nonsyn_c_hyper;

tissue_types=["Bladder","Blood","Bone","Brain","Breast","Cervix","Colorectal","Esophagus","Gall Bladder","Head & neck","Kidney","Liver","Lung","Mesenchymal","Ovary","Pancreas","Prostate","Skin","Stomach","Uterus"];

C=load('Path_To\Data_Extracted_files\codon_3mer_CancerGenes_synNonsyn.mat');
C_c=C.C_nonsyn_c;

Gc=load('Path_To\Data_Extracted_files\Germline_synNonsyn_CancerGenes.mat');
G_c=Gc.G_nonsyn_c;


% 1-mer counts
ind_1mer=[1:16 ; 17:32 ; 33:48 ; 49:64 ; 65:80 ; 81:96];
onemer_nonhyper_tissues_count=zeros(6,20);
onemer_hyper_tissues_count=zeros(6,20);
onemer_germline_count=zeros(6,1);
onemer_genome_count=zeros(6,1);
for i=1:6
    onemer_germline_count(i)=sum(G_c(ind_1mer(i,:)));
    onemer_genome_count(i)=sum(C_c(ind_1mer(i,:)));
end
for j=1:20
    for i=1:6
        onemer_nonhyper_tissues_count(i,j)=sum(P_nonhyper_tissues_c(ind_1mer(i,:),j));
        onemer_hyper_tissues_count(i,j)=sum(P_hyper_tissues_c(ind_1mer(i,:),j));
    end
end

% Ts:Tv counts
Ts_nonhyper_tissues_count=onemer_nonhyper_tissues_count(3,:)+onemer_nonhyper_tissues_count(5,:);
Tv_nonhyper_tissues_count=onemer_nonhyper_tissues_count(1,:)+onemer_nonhyper_tissues_count(2,:)+onemer_nonhyper_tissues_count(4,:)+onemer_nonhyper_tissues_count(6,:);
Ts_hyper_tissues_count=onemer_hyper_tissues_count(3,:)+onemer_hyper_tissues_count(5,:);
Tv_hyper_tissues_count=onemer_hyper_tissues_count(1,:)+onemer_hyper_tissues_count(2,:)+onemer_hyper_tissues_count(4,:)+onemer_hyper_tissues_count(6,:);
Ts_germline_count=onemer_germline_count(3)+onemer_germline_count(5);
Tv_germline_count=onemer_germline_count(1)+onemer_germline_count(2)+onemer_germline_count(4)+onemer_germline_count(6);
Ts_genome_count=onemer_genome_count(3)+onemer_genome_count(5);
Tv_genome_count=onemer_genome_count(1)+onemer_genome_count(2)+onemer_genome_count(4)+onemer_genome_count(6);


% Hyper all together and NHM all together
P_nonhyper_c=sum(P_nonhyper_tissues_c,2);
P_hyper_c=sum(P_hyper_tissues_c,2);
onemer_nonhyper_count=sum(onemer_nonhyper_tissues_count,2);
onemer_hyper_count=sum(onemer_hyper_tissues_count,2);
Ts_nonhyper_count=sum(Ts_nonhyper_tissues_count,2);
Tv_nonhyper_count=sum(Tv_nonhyper_tissues_count,2);
Ts_hyper_count=sum(Ts_hyper_tissues_count,2);
Tv_hyper_count=sum(Tv_hyper_tissues_count,2);

spect_3mer_nonhyper_c=P_nonhyper_c/sum(P_nonhyper_c);
spect_3mer_hyper_c=P_hyper_c/sum(P_hyper_c);
spect_3mer_germline_c=G_c/sum(G_c);
spect_1mer_nonhyper_c=onemer_nonhyper_count/sum(onemer_nonhyper_count);
spect_1mer_hyper_c=onemer_hyper_count/sum(onemer_hyper_count);
spect_1mer_germline_c=onemer_germline_count/sum(onemer_germline_count);
spect_Ts_nonhyper_c=Ts_nonhyper_count/(Ts_nonhyper_count+Tv_nonhyper_count);
spect_Ts_hyper_c=Ts_hyper_count/(Ts_hyper_count+Tv_hyper_count);
spect_Ts_germline_c=Ts_germline_count/(Ts_germline_count+Tv_germline_count);

spect_3mer_genome_c=C_c/sum(C_c);
spect_1mer_genome_c=onemer_genome_count/sum(onemer_genome_count);
spect_Ts_genome_c=Ts_genome_count/(Ts_genome_count+Tv_genome_count);


% Find bias reversal measure for above spectra


rev_3mer_hyper=rev_counts(spect_3mer_hyper_c,spect_3mer_germline_c,spect_3mer_genome_c); 
rev_3mer_nonhyper=rev_counts(spect_3mer_nonhyper_c,spect_3mer_germline_c,spect_3mer_genome_c);

rev_1mer_hyper=rev_counts(spect_1mer_hyper_c,spect_1mer_germline_c,spect_1mer_genome_c); 
rev_1mer_nonhyper=rev_counts(spect_1mer_nonhyper_c,spect_1mer_germline_c,spect_1mer_genome_c);

rev_Ts_hyper=rev_counts([spect_Ts_hyper_c;1-spect_Ts_hyper_c],[spect_Ts_germline_c;1-spect_Ts_germline_c],[spect_Ts_genome_c;1-spect_Ts_genome_c]); 
rev_Ts_nonhyper=rev_counts([spect_Ts_nonhyper_c;1-spect_Ts_nonhyper_c],[spect_Ts_germline_c;1-spect_Ts_germline_c],[spect_Ts_genome_c;1-spect_Ts_genome_c]);


