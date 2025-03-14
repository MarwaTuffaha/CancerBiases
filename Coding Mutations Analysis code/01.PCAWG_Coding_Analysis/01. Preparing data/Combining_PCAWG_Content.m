% Classifying PCAWG coding mutations using content information
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
C=load('Path_To\Data_Extracted_files\Codon_Content_Mutations');
C=C.R;

% Column1  chromosome
% Column2  position
% Column3  reference ntd
% Column4  alternative ntd
% Column5  degeneracy
% Column6  tissue type (using index of tissue type as a proxy)
% Column7  Cancer gene? (yes 1, no 0)
% Column8  Donor ID
P=load('Path_To\Data_Extracted_files\PCAWG_Coding_Mutations');
P=P.P;

tissue_types=["Bladder","Blood","Bone","Brain","Breast","Cervix","Colorectal","Esophagus","Gall Bladder","Head & neck","Kidney","Liver","Lung","Mesenchymal","Ovary","Pancreas","Prostate","Skin","Stomach","Uterus"];

% Sanity check that all PCAWG coding mutations in P are found in the the
% genome C, and storing their indeces in C
t=0; % counts rows in P that were found in C
J=[]; % stores their indeces
Ref=[]; % tests is the reference ntd is the same
Deg=[]; % tests if the degeneracy of the mutation is the same

N=[]; % Storing indeces of mutations in P that were not found in C

k=1; % index on C

for i=1:length(P)
    for j=k:length(C)
        if P(i,1)==C(j,1) && P(i,2)==C(j,2) && P(i,4)==C(j,4)
            t=t+1;
            J=[J,j];
            if P(i,3)==C(j,3)
                r_t=1;
            else
                r_t=0;
            end
            Ref=[Ref,r_t];
            if P(i,5)==C(j,7)
                d_t=1;
            else
                d_t=0;
            end
            Deg=[Deg,d_t];
            k=j;
            break
        end
        if j==length(C)
            N=[N,i];
        end
    end
end

%%
% These should be equal to the length of P
sum(Deg)
sum(Ref)
t
%%
% There are 2022 mutations (0.0048) in P whose details were not found in C
% We will ignore them - Let's delete the corresponding rows
n=length(N);
for i=1:n
    P(N(n+1-i),:)=[];
end

%%

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
% Column11 Cancer gene? (yes 1, no 0)
% Column12 Donor ID

H=zeros(length(P),12); % combined matrix
H(:,1:4)=P(:,1:4);
H(:,7)=P(:,5);
H(:,10)=P(:,6);
H(:,11)=P(:,7);
H(:,12)=P(:,8);

for m=1:length(J)
    H(m,5)=C(J(m),5);
    H(m,6)=C(J(m),6);
    H(m,8)=C(J(m),8);
    H(m,9)=C(J(m),9);
end

%%

% Fix the 3mer column
% Unifying the complement mutations by using function comp_3mer
h=length(H);
HC=zeros(h,12);
for i=1:h
    HC(i,1)=H(i,1);
    HC(i,2)=H(i,2);
    HC(i,6)=H(i,6);
    HC(i,7)=H(i,7);
    HC(i,8)=H(i,8);
    HC(i,9)=H(i,9);
    HC(i,10)=H(i,10);
    HC(i,11)=H(i,11);
    HC(i,12)=H(i,12);
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
% Export this final version
save('Path_To\Data_Extracted_files\Labeled_Coding_PCAWG','H')


