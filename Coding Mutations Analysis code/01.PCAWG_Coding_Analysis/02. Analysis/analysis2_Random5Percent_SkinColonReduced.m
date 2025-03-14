% analysis of the combined file on the 3mer level: Labeled coding PCAWG mutations
% Choose random 5% of the mutations to compare with cancer gene mutation
% results

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


H1=load('Path_To\Data_Extracted_files\Labeled_Coding_PCAWG_hyperNonhyper.mat');
H_nonhyper=H1.H_nonhyper;
h_nonhyper=length(H_nonhyper);

load('Path_To\Data_Extracted_files\Labeled_Coding_PCAWG_HM_SkinColonReduced.mat',"H_hyper_reduced")
H_hyper=H_hyper_reduced;
h_hyper=length(H_hyper);

tissue_types=["Bladder","Blood","Bone","Brain","Breast","Cervix","Colorectal","Esophagus","Gall Bladder","Head & neck","Kidney","Liver","Lung","Mesenchymal","Ovary","Pancreas","Prostate","Skin","Stomach","Uterus"];



%%

mers=load('3mers.mat');
U=mers.U;
n=length(U);

%%

% Repeat the following x times
% Choose random 5% of the mutations and count syn/nonsyn
% Store the results in the following 2 variables to be exported in the end

% Do for all mutations together, and for hyper and nonhyper mutations separated

x=50000;

B_syn_hyper=zeros(n,x); B_nonsyn_hyper=zeros(n,x); 

for k=1:x
    % Choose random 5% of the mutations
    number= round(0.05*h_hyper);
    r=zeros(number,1); %indeces
    for i=1:number
        r(i) = randi([1 h_hyper],1);
    end
    
    % Delete other mutations
    HD=H_hyper;
    HD=HD(r,:);

    % Counting ts/tv - syn/nonsyn mutations for each 3mer mutation
    for i=1:number
        for j=1:n
            if HD(i,5)==U(j,1) && HD(i,4)==U(j,2)
                if HD(i,9)==0
                    B_syn_hyper(j,k)=B_syn_hyper(j,k)+1;
                else
                    B_nonsyn_hyper(j,k)=B_nonsyn_hyper(j,k)+1;
                end            
            end
        end
    end 
end

B_syn_nonhyper=zeros(n,x); B_nonsyn_nonhyper=zeros(n,x);

for k=1:x
    % Choose random 5% of the mutations
    number= round(0.05*h_nonhyper);
    r=zeros(number,1); %indeces
    for i=1:number
        r(i) = randi([1 h_nonhyper],1);
    end
    
    % Delete other mutations
    HD=H_nonhyper;
    HD=HD(r,:);

    % Counting ts/tv - syn/nonsyn mutations for each 3mer mutation
    for i=1:number
        for j=1:n
            if HD(i,5)==U(j,1) && HD(i,4)==U(j,2)
                if HD(i,9)==0
                    B_syn_nonhyper(j,k)=B_syn_nonhyper(j,k)+1;
                else
                    B_nonsyn_nonhyper(j,k)=B_nonsyn_nonhyper(j,k)+1;
                end            
            end
        end
    end     
end


save('Path_To\Data_Extracted_files\PCAWG_3mer_synNonsyn_Random5Percent_SkinColonReduced','B_syn_hyper','B_nonsyn_hyper','B_syn_nonhyper','B_nonsyn_nonhyper')

