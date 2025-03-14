% Reduce the number of mutations skin and colon hypermutated samples 
% to avoid dominance by these two tissues in the pooled HM spectrum
% We will choose only a number of mutations equal to the maximum number of
% HM mutations in other tissues = stomach


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

H=load('Path_To\Data_Extracted_files\Labeled_Coding_PCAWG_hyperNonhyper.mat');
H=H.H_hyper;
h=length(H);

H_skin=H(H(:,10)==18,:);
h_skin=length(H_skin);
H_colon=H(H(:,10)==7,:);
h_colon=length(H_colon);

H_stomach=H(H(:,10)==19,:);
s=length(H_stomach); % size of the sample

% random indeces
p_skin = randsample(h_skin, s);
p_colon = randsample(h_colon, s);

% reduced matrices
H_skin_r = H_skin(p_skin, :); 
H_colon_r = H_colon(p_colon, :);

H(H(:,10)==18,:)=[];
H(H(:,10)==7,:)=[];
H_hyper_reduced=[H;H_skin_r;H_colon_r];

save('Path_To\Data_Extracted_files\Labeled_Coding_PCAWG_HM_SkinColonReduced.mat',"H_hyper_reduced")


%% Now find spectra
mers=load('Path_To\Data_Extracted_files\3mers.mat');
U=mers.U;
n=length(U);

P_skin_HM_syn=zeros(n,1);
P_skin_HM_nonsyn=zeros(n,1);

P_colon_HM_syn=zeros(n,1);
P_colon_HM_nonsyn=zeros(n,1);

for i=1:s
    for j=1:n
        if H_skin_r(i,5)==U(j,1) && H_skin_r(i,4)==U(j,2)
            if H_skin_r(i,9)==0
                P_skin_HM_syn(j)=P_skin_HM_syn(j)+1;
            else
                P_skin_HM_nonsyn(j)=P_skin_HM_nonsyn(j)+1;
            end            
        end
        if H_colon_r(i,5)==U(j,1) && H_colon_r(i,4)==U(j,2)
            if H_colon_r(i,9)==0
                P_colon_HM_syn(j)=P_colon_HM_syn(j)+1;
            else
                P_colon_HM_nonsyn(j)=P_colon_HM_nonsyn(j)+1;
            end            
        end        
    end
end


% distinguish cancer-gene from noncancer-gene mutations
P_skin_HM_syn_c=zeros(n,1);
P_skin_HM_nonsyn_c=zeros(n,1);
P_skin_HM_syn_nc=zeros(n,1);
P_skin_HM_nonsyn_nc=zeros(n,1);
P_colon_HM_syn_c=zeros(n,1);
P_colon_HM_nonsyn_c=zeros(n,1);
P_colon_HM_syn_nc=zeros(n,1);
P_colon_HM_nonsyn_nc=zeros(n,1);

for i=1:s
    for j=1:n
        if H_skin_r(i,11)==1 % cancer gene
            if H_skin_r(i,5)==U(j,1) && H_skin_r(i,4)==U(j,2)
                if H_skin_r(i,9)==0
                    P_skin_HM_syn_c(j)=P_skin_HM_syn_c(j)+1;
                else
                    P_skin_HM_nonsyn_c(j)=P_skin_HM_nonsyn_c(j)+1;
                end            
            end
        else % noncancer gene
            if H_skin_r(i,5)==U(j,1) && H_skin_r(i,4)==U(j,2)
                if H_skin_r(i,9)==0
                    P_skin_HM_syn_nc(j)=P_skin_HM_syn_nc(j)+1;
                else
                    P_skin_HM_nonsyn_nc(j)=P_skin_HM_nonsyn_nc(j)+1;
                end            
            end            
        end
        if H_colon_r(i,11)==1 % cancer gene
            if H_colon_r(i,5)==U(j,1) && H_colon_r(i,4)==U(j,2)
                if H_colon_r(i,9)==0
                    P_colon_HM_syn_c(j)=P_colon_HM_syn_c(j)+1;
                else
                    P_colon_HM_nonsyn_c(j)=P_colon_HM_nonsyn_c(j)+1;
                end            
            end  
        else % noncancer gene
            if H_colon_r(i,5)==U(j,1) && H_colon_r(i,4)==U(j,2)
                if H_colon_r(i,9)==0
                    P_colon_HM_syn_nc(j)=P_colon_HM_syn_nc(j)+1;
                else
                    P_colon_HM_nonsyn_nc(j)=P_colon_HM_nonsyn_nc(j)+1;
                end            
            end  
        end
    end
end

save("Path_To\Data_Extracted_files\Reduced_SkinColonHM","P_colon_HM_nonsyn_nc","P_colon_HM_syn_nc","P_colon_HM_nonsyn_c","P_colon_HM_syn_c","P_skin_HM_nonsyn_nc","P_skin_HM_syn_nc","P_skin_HM_nonsyn_c","P_skin_HM_syn_c")
