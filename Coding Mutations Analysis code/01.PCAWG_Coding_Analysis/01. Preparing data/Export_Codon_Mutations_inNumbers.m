%Converting genome coding content data into numerical form to reduce data memory size
%and extract mutations only in cancer genes and those only in passenger genes

fileID = fopen("Path_to\coding_sites_annotations_effects_3mers.txt");
% overall # rows = 94,132,026
% Column1  chromosome
% Column2  position
% Column3  reference ntd
% Column4  reference 3-mer
% Column5  alternative ntd
% Column6  strand           ignore for now
% Column7  codon            ignore for now
% Column8  gene name        ignore for now
% Column9  Gene ID          ignore for now
% Column10  transcript ID   ignore for now
% Column11  cancer gene         -> column 10
% Column12  degeneracy          -> column 7
% Column13  mutation effect     -> column 6

C = textscan(fileID,'%f %f %s %s %s %*f %*s %*s %*s %*s %s %f %s','CommentStyle','chr'); 

fclose(fileID);

h=cellfun(@height,C(1));

%%

% delete the rows that have unknown 3mer context
for i=1:8
    C{i}(72709657,:)=[];
    C{i}(48832043,:)=[];
    C{i}(1775038,:)=[];
end

h=cellfun(@height,C(1));


%%
% Create a large matrix R and save the data numerically to save space
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
% Column10 0:noncancer gene | 1:cancer gene


R=zeros(h,10);

for j=1:h
    R(j,1)=C{1}(j);
    R(j,2)=C{2}(j);
    R(j,3)=let2num(char(C{3}(j)));
    R(j,4)=let2num(char(C{5}(j)));
    a=char(C{4}(j));
    R(j,5)=let2num(a(3))+10*let2num(a(2))+100*let2num(a(1));
    R(j,6)=effect2num(string(C{8}(j)));
    R(j,7)=C{7}(j);   
    if or(or(char(C{3}(j))=='A' & char(C{5}(j))=='G' ,char(C{3}(j))=='G' & char(C{5}(j))=='A'),or(char(C{3}(j))=='C' & char(C{5}(j))=='T' ,char(C{3}(j))=='T' & char(C{5}(j))=='C'))
        R(j,8)=1;
    end
    if or(string(C{8}(j))=="MISSENSE",string(C{8}(j))=="NONSENSE")
        R(j,9)=1;
    end   
    if string(C{6}(j))=="yes"
        R(j,10)=1;
    end
end

%%

R = sortrows(R,2);
R = sortrows(R,1);

%%

save("Path_To\Data_Extracted_files\Codon_Content_Mutations",'R', '-v7.3')


%% Extract only cancer genes and only noncancer genes

a=find(R(:,10)==0);
b=find(R(:,10)==1);

RC=R;
RC(a,:)=[];

RNC=R;
RNC(b,:)=[];


save("Path_To\Data_Extracted_files\Codon_Content_CancerGenes_Mutations",'RC','RNC', '-v7.3')


%%

Syn_Tv=0;
Syn_Ts=0;
NonSyn_Tv=0;
NonSyn_Ts=0;

for j=1:h
    if (R(j,8)==0 && R(j,9)==0)
        Syn_Tv=Syn_Tv+1;
    elseif (R(j,8)==1 && R(j,9)==0)
        Syn_Ts=Syn_Ts+1;    
    elseif (R(j,8)==0 && R(j,9)==1)
        NonSyn_Tv=NonSyn_Tv+1;
    elseif (R(j,8)==1 && R(j,9)==1)
        NonSyn_Ts=NonSyn_Ts+1;        
    end
end

all_mut=Syn_Tv+Syn_Ts+NonSyn_Tv+NonSyn_Ts;
Syn=Syn_Tv+Syn_Ts;
NonSyn=NonSyn_Tv+NonSyn_Ts;
Tv=Syn_Tv+NonSyn_Tv;
Ts=Syn_Ts+NonSyn_Ts;

frac_syn=Syn/h;
frac_nonsyn=NonSyn/h;
frac_tv=Tv/h;
frac_ts=Ts/h;

frac_tv_of_syn=Syn_Tv/Syn;
frac_ts_of_syn=Syn_Ts/Syn;
frac_tv_of_nonsyn=NonSyn_Tv/NonSyn;
frac_ts_of_nonsyn=NonSyn_Ts/NonSyn;

frac_syn_of_tv=Syn_Tv/Tv;
frac_syn_of_ts=Syn_Ts/Ts;
frac_nonsyn_of_tv=NonSyn_Tv/Tv;
frac_nonsyn_of_ts=NonSyn_Ts/Ts;

%%

figure
subplot(1,2,1)
X = categorical(["All Mutations","Synonymous","Nonsynonymous"]);
X = reordercats(X,["All Mutations","Synonymous","Nonsynonymous"]);
bar(X,[frac_tv,frac_ts;frac_tv_of_syn,frac_ts_of_syn;frac_tv_of_nonsyn,frac_ts_of_nonsyn],'stacked')
title("trv:trs")
legend("Transversions","Transitions")

subplot(1,2,2)
X = categorical(["All Mutations","Transversions","Transitions"]);
X = reordercats(X,["All Mutations","Transversions","Transitions"]);
bar(X,[frac_syn,frac_nonsyn;frac_syn_of_tv,frac_nonsyn_of_tv;frac_syn_of_ts,frac_nonsyn_of_ts],'stacked')
title("syn:nonsyn")
legend("Synonymous","Nonsynonymous")


%%

% Test code

syn_ACT_A=[];
syn_ACT_G=[];
syn_ATT_A=[];
syn_ATT_G=[];

syn_AGT_T=[];
syn_AGT_C=[];
syn_AAT_T=[];
syn_AAT_C=[];

for i=1:h
    if string(C{7}(i))=="SILENT"
        if string(C{4}(i))=="ACT" && string(C{5}(i))=="A"
            syn_ACT_A=[syn_ACT_A,i];
        elseif string(C{4}(i))=="ACT" && string(C{5}(i))=="G"
            syn_ACT_G=[syn_ACT_G,i];
        elseif string(C{4}(i))=="ATT" && string(C{5}(i))=="A"
            syn_ATT_A=[syn_ATT_A,i];
        elseif string(C{4}(i))=="ATT" && string(C{5}(i))=="G"
            syn_ATT_G=[syn_ATT_G,i];
        elseif string(C{4}(i))=="AGT" && string(C{5}(i))=="T"
            syn_AGT_T=[syn_AGT_T,i];
        elseif string(C{4}(i))=="AGT" && string(C{5}(i))=="C"
            syn_AGT_C=[syn_AGT_C,i];
        elseif string(C{4}(i))=="AAT" && string(C{5}(i))=="T"
            syn_AAT_T=[syn_AAT_T,i];
        elseif string(C{4}(i))=="AAT" && string(C{5}(i))=="C"
            syn_AAT_C=[syn_AAT_C,i];            
        end
    end
end

%%

syn_ACT_A=[];
syn_ACT_G=[];
syn_ATT_A=[];
syn_ATT_G=[];

r=length(R);

for i=1:r
    if R(i,9)==0
        if R(i,5)==124 && R(i,4)==1
            syn_ACT_A=[syn_ACT_A,i];
        elseif R(i,5)==124 && R(i,4)==3
            syn_ACT_G=[syn_ACT_G,i];
        elseif R(i,5)==144 && R(i,4)==1
            syn_ATT_A=[syn_ATT_A,i];
        elseif R(i,5)==144 && R(i,4)==3
            syn_ATT_G=[syn_ATT_G,i];        
        end
    end
end

%% fixing 3mer columns
% Unifying the complement mutations by using function comp_3mer

load("Path_To\Data_Extracted_files\Codon_Content_Mutations",'R')

r=length(R);

for i=1:r
    if or(R(i,3)==1,R(i,3)==3)
        [comp_ref,comp_alt,comp_mer]=comp_3mer(R(i,3),R(i,4),R(i,5));
        R(i,3)=comp_ref;
        R(i,4)=comp_alt;
        R(i,5)=comp_mer;
    end
end

save("Path_To\Data_Extracted_files\Codon_Content_Mutations_3merFixed",'R', '-v7.3')
