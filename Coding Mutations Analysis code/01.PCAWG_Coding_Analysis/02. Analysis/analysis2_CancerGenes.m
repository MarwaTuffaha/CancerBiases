% analysis of the combined file on the 3mer level: Labeled coding PCAWG mutations

% Column1  chromosome
% Column2  position
% Column3  reference ntd
% Column4  alternative ntd
% Column5  reference 3-mer
% Column6  mutation effect
% Column7  degeneracy
% Column8  0:trv | 1:trs
% Column9  0:syn | 1:nonsyn
% Column10 tissue type
% Column11  Cancer gene? (yes 1, no 0)
% Column12 Donor ID

H=load('Path_To\Data_Extracted_files\Labeled_Coding_PCAWG.mat');
H=H.H;
h=length(H);

tissue_types=["Bladder","Blood","Bone","Brain","Breast","Cervix","Colorectal","Esophagus","Gall Bladder","Head & neck","Kidney","Liver","Lung","Mesenchymal","Ovary","Pancreas","Prostate","Skin","Stomach","Uterus"];
t=length(tissue_types);


%%

mers=load('3mers.mat');
U=mers.U;
n=length(U);

%%
% Counting all-gene and cancer-gene mutations in each tissue type
% (all together - regardless of 3mers)
count_all_tissues_syn=zeros(1,20);
count_all_tissues_nonsyn=zeros(1,20);
count_CancerGene_tissues_syn=zeros(1,20);
count_CancerGene_tissues_nonsyn=zeros(1,20);

for i=1:20
    count_all_tissues_syn(i)=sum(H(:,10) == i & H(:,9)==0);
    count_all_tissues_nonsyn(i)=sum(H(:,10) == i & H(:,9)==1);
    count_CancerGene_tissues_syn(i)=sum(H(:,10) == i & H(:,9)==0 & H(:,11)==1);
    count_CancerGene_tissues_nonsyn(i)=sum(H(:,10) == i & H(:,9)==1 & H(:,11)==1);
end

save("Path_To\Data_Extracted_files\Count_3mer_mutations_tissues","count_all_tissues_syn","count_all_tissues_nonsyn","count_CancerGene_tissues_syn","count_CancerGene_tissues_nonsyn")


%% 
% Counting ts/tv - syn/nonsyn mutations for each 3mer mutation
% and classifying syn/nonsyn mutations to being in cancer gene or not
P_syn_c=zeros(n,t); %syn cancer gene
P_syn_nc=zeros(n,t); %syn non-cancer gene
P_nonsyn_c=zeros(n,t); %nonsyn cancer gene
P_nonsyn_nc=zeros(n,t); %nonsyn non-cancer gene
P_tv=zeros(n,t);
P_ts=zeros(n,t);


for i=1:h
    for j=1:n
        if H(i,5)==U(j,1) && H(i,4)==U(j,2)
            if H(i,8)==0
                P_tv(j,H(i,10))=P_tv(j,H(i,10))+1;
            else
                P_ts(j,H(i,10))=P_ts(j,H(i,10))+1;
            end
            if H(i,9)==0 && H(i,11)==1
                P_syn_c(j,H(i,10))=P_syn_c(j,H(i,10))+1;
            elseif H(i,9)==0 && H(i,11)==0
                P_syn_nc(j,H(i,10))=P_syn_nc(j,H(i,10))+1;
            elseif H(i,9)==1 && H(i,11)==1
                P_nonsyn_c(j,H(i,10))=P_nonsyn_c(j,H(i,10))+1; 
            else
                P_nonsyn_nc(j,H(i,10))=P_nonsyn_nc(j,H(i,10))+1; 
            end            
        end
    end
end

save('Path_To\Data_Extracted_files\PCAWG_3mer_synNonsyn_Tissues_cancerGenes','P_syn_c','P_syn_nc','P_nonsyn_c','P_nonsyn_nc')

%% Now count cancer and noncancer gene synonymous and nonsynonymous mutations in hyper and nonhyper samples

load('Path_To\Data_Extracted_files\Labeled_Coding_PCAWG_hyperNonhyper','H_hyper','H_nonhyper')

P_syn_c_hyper=zeros(n,t);
P_syn_c_nonhyper=zeros(n,t);
P_nonsyn_c_hyper=zeros(n,t);
P_nonsyn_c_nonhyper=zeros(n,t);
P_syn_nc_hyper=zeros(n,t);
P_syn_nc_nonhyper=zeros(n,t);
P_nonsyn_nc_hyper=zeros(n,t);
P_nonsyn_nc_nonhyper=zeros(n,t);

h=length(H_hyper);

for i=1:h
    for j=1:n
        if H_hyper(i,5)==U(j,1) && H_hyper(i,4)==U(j,2)
            if H_hyper(i,9)==1 && H_hyper(i,11)==1
                P_nonsyn_c_hyper(j,H_hyper(i,10))=P_nonsyn_c_hyper(j,H_hyper(i,10))+1; 
            elseif H_hyper(i,9)==1 && H_hyper(i,11)==0
                P_nonsyn_nc_hyper(j,H_hyper(i,10))=P_nonsyn_nc_hyper(j,H_hyper(i,10))+1; 
            elseif H_hyper(i,9)==0 && H_hyper(i,11)==1
                P_syn_c_hyper(j,H_hyper(i,10))=P_syn_c_hyper(j,H_hyper(i,10))+1; 
            else
                P_syn_nc_hyper(j,H_hyper(i,10))=P_syn_nc_hyper(j,H_hyper(i,10))+1; 
            end            
        end
    end
end

h=length(H_nonhyper);

for i=1:h
    for j=1:n
        if H_nonhyper(i,5)==U(j,1) && H_nonhyper(i,4)==U(j,2)
            if H_nonhyper(i,9)==1 && H_nonhyper(i,11)==1
                P_nonsyn_c_nonhyper(j,H_nonhyper(i,10))=P_nonsyn_c_nonhyper(j,H_nonhyper(i,10))+1; 
            elseif H_nonhyper(i,9)==1 && H_nonhyper(i,11)==0
                P_nonsyn_nc_nonhyper(j,H_nonhyper(i,10))=P_nonsyn_nc_nonhyper(j,H_nonhyper(i,10))+1; 
            elseif H_nonhyper(i,9)==0 && H_nonhyper(i,11)==1
                P_syn_c_nonhyper(j,H_nonhyper(i,10))=P_syn_c_nonhyper(j,H_nonhyper(i,10))+1; 
            else
                P_syn_nc_nonhyper(j,H_nonhyper(i,10))=P_syn_nc_nonhyper(j,H_nonhyper(i,10))+1; 
            end            
        end
    end
end


H_hyper_skin=H_hyper(H_hyper(:,10)==18,:);
H_hyper_colon=H_hyper(H_hyper(:,10)==7,:);

P_syn_c_hyper_skin=zeros(n,t);
P_nonsyn_c_hyper_skin=zeros(n,t);
P_syn_nc_hyper_skin=zeros(n,t);
P_nonsyn_nc_hyper_skin=zeros(n,t);
P_syn_c_hyper_colon=zeros(n,t);
P_nonsyn_c_hyper_colon=zeros(n,t);
P_syn_nc_hyper_colon=zeros(n,t);
P_nonsyn_nc_hyper_colon=zeros(n,t);


h=length(H_hyper_skin);

for i=1:h
    for j=1:n
        if H_hyper_skin(i,5)==U(j,1) && H_hyper_skin(i,4)==U(j,2)
            if H_hyper_skin(i,9)==1 && H_hyper_skin(i,11)==1
                P_nonsyn_c_hyper_skin(j,H_hyper_skin(i,10))=P_nonsyn_c_hyper_skin(j,H_hyper_skin(i,10))+1; 
            elseif H_hyper_skin(i,9)==1 && H_hyper_skin(i,11)==0
                P_nonsyn_nc_hyper_skin(j,H_hyper_skin(i,10))=P_nonsyn_nc_hyper_skin(j,H_hyper_skin(i,10))+1; 
            elseif H_hyper_skin(i,9)==0 && H_hyper_skin(i,11)==1
                P_syn_c_hyper_skin(j,H_hyper_skin(i,10))=P_syn_c_hyper_skin(j,H_hyper_skin(i,10))+1; 
            else
                P_syn_nc_hyper_skin(j,H_hyper_skin(i,10))=P_syn_nc_hyper_skin(j,H_hyper_skin(i,10))+1; 
            end            
        end
    end
end

h=length(H_hyper_colon);

for i=1:h
    for j=1:n
        if H_hyper_colon(i,5)==U(j,1) && H_hyper_colon(i,4)==U(j,2)
            if H_hyper_colon(i,9)==1 && H_hyper_colon(i,11)==1
                P_nonsyn_c_hyper_colon(j,H_hyper_colon(i,10))=P_nonsyn_c_hyper_colon(j,H_hyper_colon(i,10))+1; 
            elseif H_hyper_colon(i,9)==1 && H_hyper_colon(i,11)==0
                P_nonsyn_nc_hyper_colon(j,H_hyper_colon(i,10))=P_nonsyn_nc_hyper_colon(j,H_hyper_colon(i,10))+1; 
            elseif H_hyper_colon(i,9)==0 && H_hyper_colon(i,11)==1
                P_syn_c_hyper_colon(j,H_hyper_colon(i,10))=P_syn_c_hyper_colon(j,H_hyper_colon(i,10))+1; 
            else
                P_syn_nc_hyper_colon(j,H_hyper_colon(i,10))=P_syn_nc_hyper_colon(j,H_hyper_colon(i,10))+1; 
            end            
        end
    end
end


save('Path_To\Data_Extracted_files\PCAWG_3mer_Nonsyn_Tissues_HyperNonhyper_cancerGenes','P_syn_c_hyper','P_syn_c_nonhyper','P_nonsyn_c_hyper','P_nonsyn_c_nonhyper','P_syn_nc_hyper','P_syn_nc_nonhyper','P_nonsyn_nc_hyper','P_nonsyn_nc_nonhyper','P_syn_c_hyper_skin','P_nonsyn_c_hyper_skin','P_syn_nc_hyper_skin','P_nonsyn_nc_hyper_skin','P_syn_c_hyper_colon','P_nonsyn_c_hyper_colon','P_syn_nc_hyper_colon','P_nonsyn_nc_hyper_colon')


