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

tissue_types=["Bladder","Blood","Bone","Brain","Breast","Cervix","Colorectal","Esophagus","Gall Bladder","Head & neck","Kidney","Liver","Lung","Mesenchymal","Ovary","Pancreas","Prostate","Skin","Stomach","Uterus"];
h=length(H);

mers=load('3mers.mat');
U=mers.U;
n=length(U);

% Counting ts/tv - syn/nonsyn mutations for each 3mer mutation
P_syn=zeros(n,1);
P_nonsyn=zeros(n,1);
P_tv=zeros(n,1);
P_ts=zeros(n,1);

for i=1:h
    for j=1:n
        if H(i,5)==U(j,1) && H(i,4)==U(j,2)
            if H(i,8)==0
                P_tv(j)=P_tv(j)+1;
            else
                P_ts(j)=P_ts(j)+1;
            end
            if H(i,9)==0
                P_syn(j)=P_syn(j)+1;
            else
                P_nonsyn(j)=P_nonsyn(j)+1;
            end            
        end
    end
end

save('Path_To\Data_Extracted_files\PCAWG_3mer_synNonsyn','P_syn','P_nonsyn')


%% Now let's weight by the number of mutations in each tissue type

%counts=load("tissue_mutations_count.mat");
%counts=counts.counts;

% First, count number of mutations for each 3mer in each tissue
t=length(tissue_types);

counts=zeros(n,t); % rows: 3mers, columns: tissues
counts_syn=zeros(n,t);
counts_nonsyn=zeros(n,t);

for i=1:h % loop on mutations
    for j=1:n % loop on 3mers
        if (H(i,5)==U(j,1) && H(i,4)==U(j,2)) % matching the 3mer mutation
            counts(j,H(i,10))=counts(j,H(i,10))+1;
            if H(i,9)==0
                counts_syn(j,H(i,10))=counts_syn(j,H(i,10))+1;
            else
                counts_nonsyn(j,H(i,10))=counts_nonsyn(j,H(i,10))+1;
            end
        end
    end
end

save("Path_To\Data_Extracted_files\PCAWG_3mer_synNonsyn_tissues","counts_syn","counts_nonsyn")


%% Count syn/nonsyn mutations in hyper/nonhypermutated samples

fileID = fopen("Path_to\PCAWG\hypermutated.txt"); % ID of hypermutated samples
C = textscan(fileID,'%s');
fclose(fileID);

h=cellfun(@height,C(1));

Donor_ID=zeros(h,1);
for i=1:h
    Donor=char(C{1}(i)); Donor(1:2)=[];
    Donor_ID(i)=str2double(Donor);
end

H_nonhyper=H;
H_hyper=[];

for i=1:h
    index=find(H_nonhyper(:,12)==Donor_ID(i));
    for j=1:length(index)
        H_hyper=[H_hyper;H_nonhyper(index(j),:)];
    end
    H_nonhyper(index,:)=[];
end

H_noHyperExceptSkin=H;

for i=1:h
    index=find(H_noHyperExceptSkin(:,12)==Donor_ID(i));
    if ~isempty(index)
        if H_noHyperExceptSkin(index(1),10)~=18
            H_noHyperExceptSkin(index,:)=[];
        end
    end
end

H_noSkin=H;
index=find(H_noSkin(:,10)==18);
H_noSkin(index,:)=[];


H_noHyperColon=H;
for i=1:h
    index=find(H_noHyperColon(:,12)==Donor_ID(i));
    if ~isempty(index)
        if H_noHyperColon(index(1),10)==7
            H_noHyperColon(index,:)=[];
        end
    end
end

save('Path_To\Data_Extracted_files\Labeled_Coding_PCAWG_hyperNonhyper','H_hyper','H_nonhyper')

%% Now count syn/nonsyn mutations in hyper / nonhyper

h=length(H_hyper);
P_syn_hyper=zeros(n,1);
P_nonsyn_hyper=zeros(n,1);
for i=1:h % loop on mutations
    for j=1:n % loop on 3mers
        if (H_hyper(i,5)==U(j,1) && H_hyper(i,4)==U(j,2)) % matching the 3mer mutation
            if H_hyper(i,9)==0
                P_syn_hyper(j)=P_syn_hyper(j)+1;
            else
                P_nonsyn_hyper(j)=P_nonsyn_hyper(j)+1;
            end
        end
    end
end

h=length(H_nonhyper);
P_syn_nonhyper=zeros(n,1);
P_nonsyn_nonhyper=zeros(n,1);
for i=1:h % loop on mutations
    for j=1:n % loop on 3mers
        if (H_nonhyper(i,5)==U(j,1) && H_nonhyper(i,4)==U(j,2)) % matching the 3mer mutation
            if H_nonhyper(i,9)==0
                P_syn_nonhyper(j)=P_syn_nonhyper(j)+1;
            else
                P_nonsyn_nonhyper(j)=P_nonsyn_nonhyper(j)+1;
            end
        end
    end
end

h=length(H_noSkin);
P_syn_noSkin=zeros(n,1);
P_nonsyn_noSkin=zeros(n,1);
for i=1:h % loop on mutations
    for j=1:n % loop on 3mers
        if (H_noSkin(i,5)==U(j,1) && H_noSkin(i,4)==U(j,2)) % matching the 3mer mutation
            if H_noSkin(i,9)==0
                P_syn_noSkin(j)=P_syn_noSkin(j)+1;
            else
                P_nonsyn_noSkin(j)=P_nonsyn_noSkin(j)+1;
            end
        end
    end
end

h=length(H_noHyperColon);
P_syn_noHyperColon=zeros(n,1);
P_nonsyn_noHyperColon=zeros(n,1);
for i=1:h % loop on mutations
    for j=1:n % loop on 3mers
        if (H_noHyperColon(i,5)==U(j,1) && H_noHyperColon(i,4)==U(j,2)) % matching the 3mer mutation
            if H_noHyperColon(i,9)==0
                P_syn_noHyperColon(j)=P_syn_noHyperColon(j)+1;
            else
                P_nonsyn_noHyperColon(j)=P_nonsyn_noHyperColon(j)+1;
            end
        end
    end
end

h=length(H_noHyperExceptSkin);
P_syn_noHyperExceptSkin=zeros(n,1);
P_nonsyn_noHyperExceptSkin=zeros(n,1);
for i=1:h % loop on mutations
    for j=1:n % loop on 3mers
        if (H_noHyperExceptSkin(i,5)==U(j,1) && H_noHyperExceptSkin(i,4)==U(j,2)) % matching the 3mer mutation
            if H_noHyperExceptSkin(i,9)==0
                P_syn_noHyperExceptSkin(j)=P_syn_noHyperExceptSkin(j)+1;
            else
                P_nonsyn_noHyperExceptSkin(j)=P_nonsyn_noHyperExceptSkin(j)+1;
            end
        end
    end
end

%%

h=length(H_nonhyper);
P_syn_nonhyper_tissues=zeros(n,20);
P_nonsyn_nonhyper_tissues=zeros(n,20);

P_ts_syn_nonhyper_tissues=zeros(n,20);
P_tv_syn_nonhyper_tissues=zeros(n,20);
P_ts_nonsyn_nonhyper_tissues=zeros(n,20);
P_tv_nonsyn_nonhyper_tissues=zeros(n,20);

for i=1:h % loop on mutations
    for j=1:n % loop on 3mers
        if (H_nonhyper(i,5)==U(j,1) && H_nonhyper(i,4)==U(j,2)) % matching the 3mer mutation
            if H_nonhyper(i,9)==0
                P_syn_nonhyper_tissues(j,H_nonhyper(i,10))=P_syn_nonhyper_tissues(j,H_nonhyper(i,10))+1;
                if H_nonhyper(i,8)==0
                    P_tv_syn_nonhyper_tissues(j,H_nonhyper(i,10))=P_tv_syn_nonhyper_tissues(j,H_nonhyper(i,10))+1;
                else
                    P_ts_syn_nonhyper_tissues(j,H_nonhyper(i,10))=P_ts_syn_nonhyper_tissues(j,H_nonhyper(i,10))+1;
                end
            else
                P_nonsyn_nonhyper_tissues(j,H_nonhyper(i,10))=P_nonsyn_nonhyper_tissues(j,H_nonhyper(i,10))+1;
                if H_nonhyper(i,8)==0
                    P_tv_nonsyn_nonhyper_tissues(j,H_nonhyper(i,10))=P_tv_nonsyn_nonhyper_tissues(j,H_nonhyper(i,10))+1;
                else
                    P_ts_nonsyn_nonhyper_tissues(j,H_nonhyper(i,10))=P_ts_nonsyn_nonhyper_tissues(j,H_nonhyper(i,10))+1;
                end
            end


        end
    end
end

h=length(H_hyper);
P_syn_hyper_tissues=zeros(n,20);
P_nonsyn_hyper_tissues=zeros(n,20);

P_ts_syn_hyper_tissues=zeros(n,20);
P_tv_syn_hyper_tissues=zeros(n,20);
P_ts_nonsyn_hyper_tissues=zeros(n,20);
P_tv_nonsyn_hyper_tissues=zeros(n,20);

for i=1:h % loop on mutations
    for j=1:n % loop on 3mers
        if (H_hyper(i,5)==U(j,1) && H_hyper(i,4)==U(j,2)) % matching the 3mer mutation
            if H_hyper(i,9)==0
                P_syn_hyper_tissues(j,H_hyper(i,10))=P_syn_hyper_tissues(j,H_hyper(i,10))+1;
                if H_hyper(i,8)==0
                    P_tv_syn_hyper_tissues(j,H_hyper(i,10))=P_tv_syn_hyper_tissues(j,H_hyper(i,10))+1;
                else
                    P_ts_syn_hyper_tissues(j,H_hyper(i,10))=P_ts_syn_hyper_tissues(j,H_hyper(i,10))+1;
                end                
            else
                P_nonsyn_hyper_tissues(j,H_hyper(i,10))=P_nonsyn_hyper_tissues(j,H_hyper(i,10))+1;
                if H_hyper(i,8)==0
                    P_tv_nonsyn_hyper_tissues(j,H_hyper(i,10))=P_tv_nonsyn_hyper_tissues(j,H_hyper(i,10))+1;
                else
                    P_ts_nonsyn_hyper_tissues(j,H_hyper(i,10))=P_ts_nonsyn_hyper_tissues(j,H_hyper(i,10))+1;
                end
            end
        end
    end
end

P_tv_nonsyn_nonhyper_tissues=sum(P_tv_nonsyn_nonhyper_tissues,1);
P_ts_nonsyn_nonhyper_tissues=sum(P_ts_nonsyn_nonhyper_tissues,1);
P_tv_nonsyn_hyper_tissues=sum(P_tv_nonsyn_hyper_tissues,1);
P_ts_nonsyn_hyper_tissues=sum(P_ts_nonsyn_hyper_tissues,1);


P_tv_syn_nonhyper_tissues=sum(P_tv_syn_nonhyper_tissues,1);
P_ts_syn_nonhyper_tissues=sum(P_ts_syn_nonhyper_tissues,1);
P_tv_syn_hyper_tissues=sum(P_tv_syn_hyper_tissues,1);
P_ts_syn_hyper_tissues=sum(P_ts_syn_hyper_tissues,1);

%%


save("Path_To\Data_Extracted_files\PCAWG_3mer_synNonsyn_hyperNonhyper","P_syn_hyper","P_nonsyn_hyper","P_syn_nonhyper","P_nonsyn_nonhyper","P_syn_noSkin","P_nonsyn_noSkin","P_syn_noHyperColon","P_nonsyn_noHyperColon", "P_syn_noHyperExceptSkin" , "P_nonsyn_noHyperExceptSkin","P_syn_nonhyper_tissues","P_nonsyn_nonhyper_tissues","P_syn_hyper_tissues","P_nonsyn_hyper_tissues")

save("Path_To\Data_Extracted_files\PCAWG_TvTs_hyperNonhyper","P_tv_nonsyn_hyper_tissues","P_ts_nonsyn_hyper_tissues","P_tv_nonsyn_nonhyper_tissues","P_ts_nonsyn_nonhyper_tissues","P_tv_syn_hyper_tissues","P_ts_syn_hyper_tissues","P_tv_syn_nonhyper_tissues","P_ts_syn_nonhyper_tissues")

