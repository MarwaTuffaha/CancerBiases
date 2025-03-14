% Count syn/nonsyn - Ts/Tv for the 96 3mer mutations in the coding genome
% Also, counting mutations in each tissue type (all together - regardless of 3mers)


mers=load('3mers.mat');
U=mers.U;
n=length(U);

% Column1  chromosome
% Column2  position
% Column3  reference ntd
% Column4  alternative ntd
% Column5  reference 3-mer
% Column6  mutation effect
% Column7  degeneracy
% Column8  0:trv | 1:trs
% Column9  0:syn | 1:nonsyn
% Column10 0:noncancer gene | 1:cancer gene

codon_muts=load('Path_to\Data_Extracted_files\Codon_Content_Mutations');
R=codon_muts.R;
r=length(R);

% First, fix the 3mer column
% Unifying the complement mutations by using function comp_3mer
RC=zeros(r,9);
for i=1:r
    RC(i,1)=R(i,1);
    RC(i,2)=R(i,2);
    RC(i,6)=R(i,6);
    RC(i,7)=R(i,7);
    RC(i,8)=R(i,8);
    RC(i,9)=R(i,9);
    if or(R(i,3)==2,R(i,3)==4)
        RC(i,3)=R(i,3);
        RC(i,4)=R(i,4);
        RC(i,5)=R(i,5);
    else
        [comp_ref,comp_alt,comp_mer]=comp_3mer(R(i,3),R(i,4),R(i,5));
        RC(i,3)=comp_ref;
        RC(i,4)=comp_alt;
        RC(i,5)=comp_mer;
    end
end


% Counting ts/tv - syn/nonsyn mutations for each 3mer mutation
C_syn=zeros(n,1);
C_nonsyn=zeros(n,1);
C_tv=zeros(n,1);
C_ts=zeros(n,1);

for i=1:r
    for j=1:n
        if RC(i,5)==U(j,1) && RC(i,4)==U(j,2)
            if RC(i,8)==0
                C_tv(j)=C_tv(j)+1;
            else
                C_ts(j)=C_ts(j)+1;
            end
            if RC(i,9)==0
                C_syn(j)=C_syn(j)+1;
            else
                C_nonsyn(j)=C_nonsyn(j)+1;
            end            
        end
    end
end



save('Path_to\Data_Extracted_files\codon_3mer_synNonsyn','C_syn','C_nonsyn')



