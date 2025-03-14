% Classifying Germline mutations using content information
% Combining the two files

% A-1   C-2   G-3   T-4

% Column1  chromosome
% Column2  position
% Column3  reference ntd
% Column4  alternative ntd
% Column5  reference 3-mer
% Column6  mutation effect
% Column7  degeneracy
% Column8  0:trv | 1:trs
% Column9  0:syn | 1:nonsyn
% Column10 Cancer gene? (yes 1, no 0)
C=load('Path_To\Codon_Content_Mutations');
C=C.R;


% Column1  chromosome
% Column2  position
% Column3  reference ntd
% Column4  alternative ntd
G=load('Path_To\Data_Extracted_files\Germline_Mutations_Numbers.mat');
G=G.G;

tissue_types=["Bladder","Blood","Bone","Brain","Breast","Cervix","Colorectal","Esophagus","Gall Bladder","Head & neck","Kidney","Liver","Lung","Mesenchymal","Ovary","Pancreas","Prostate","Skin","Stomach","Uterus"];


%%

% Check which germline mutations in G are found in the
% coding genome C, and storing their indeces in C
t=0; % counts rows in G that were found in C (those are coding)
J=zeros(1,length(G)); % stores their indeces in C
F=zeros(1,length(G)); % stores their indeces in G
Ref=zeros(1,length(G)); % tests if the reference ntd is the same

N=zeros(1,length(G)); % Storing indeces of mutations in G that were not found in C
N_i=0;

k=1; % index on C


for i=1:length(G)
    for j=k:length(C)
        if or(or(C(j,1)>G(i,1),j==length(C)),C(j,1)==G(i,1)&&C(j,2)>G(i,2)) % break loop if next chromosom or end of genome reached
            N_i=N_i+1;
            N(N_i)=i;
            break
        end

        if G(i,1)==C(j,1) && G(i,2)==C(j,2) && G(i,4)==C(j,4)
            t=t+1;
            J(t)=j;
            F(t)=i;
            if G(i,3)==C(j,3)
                r_t=1;
            else
                r_t=0;
            end
            Ref(t)=r_t;
            k=j;
            break
        end
    end
end


% These should be equal
sum(Ref)
t
% They are not!
% There are 12 mismatches in chromosome 1, probably bacause authors for the
% germline paper used an older version of the reference genome

%%
% We will remove noncoding mutations
% and mutations that did not match the reference ntd

Ref(t+1:end)=[];
F(Ref==0)=0;
F(F==0)=[];
f=length(F);
K=zeros(f,4); % new matrix
for i=1:f
    K(i,:)=G(F(i),:);
end
G=K;

%%

% Column1  chromosome
% Column2  position
% Column3  reference ntd
% Column4  alternative ntd
% Column5  reference 3-mer
% Column6  mutation effect
% Column7  0:trv | 1:trs
% Column8  0:syn | 1:nonsyn
% Column9 Cancer gene? (yes 1, no 0)  

H=zeros(length(G),9); % combined matrix
H(:,1:4)=G(:,1:4);

J(Ref==0)=0;
J(J==0)=[];

for m=1:length(J)
    H(m,5)=C(J(m),5);
    H(m,6)=C(J(m),6);
    H(m,7)=C(J(m),8);
    H(m,8)=C(J(m),9);
    H(m,9)=C(J(m),10);
end


% Fix the 3mer column
% Unifying the complement mutations by using function comp_3mer
h=length(H);
HC=zeros(h,9);
for i=1:h
    HC(i,1)=H(i,1);
    HC(i,2)=H(i,2);
    HC(i,6)=H(i,6);
    HC(i,7)=H(i,7);
    HC(i,8)=H(i,8);
    HC(i,9)=H(i,9);
    if or(H(i,3)==2,H(i,3)==4)
        HC(i,3)=H(i,3);
        HC(i,4)=H(i,4);
        HC(i,5)=H(i,5);
    else
        [comp_ref,comp_alt,comp_mer]=comp_3mer(H(i,3),H(i,4),H(i,5));
        HC(i,3)=comp_ref;
        HC(i,4)=comp_alt;
        HC(i,5)=comp_mer;
    end
end
H=HC;


%%
save('Path_To\Data_Extracted_files\Labeled_coding_Germline','H')


%% 
% Checking if C is ordered well

indeces=[];
for l=1:22
indeces=[indeces,find(C(:,1)==l,1)];
end

indeces=[indeces,length(C)];

TF=[];
for l=1:22
TF = [TF,issorted(C(indeces(l):indeces(l+1)-1,2))];
end

