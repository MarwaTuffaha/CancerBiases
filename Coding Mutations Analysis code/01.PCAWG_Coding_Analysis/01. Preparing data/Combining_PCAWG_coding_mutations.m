%Combining all PCAWG coding mutations in one file and save numerically
%Fix the 3mer column to be similar to COSMIC (pyrimidine in the middle, otherwise find the complement)


files=load('Path_to\FileDetails');
files=files.details;
filenames=files(:,1);
tissue=files(:,3);

n=length(filenames);
tissue_types=unique(tissue);

myFolder = 'Path_to_AnnotatedPCAWG_DataFiles';

C_all = cell(1,21);

%Columns:
% 1-chr	    2-stop	3-ID	4-REF	5-ALT	6-QUAL	7-MUSE	8-AD_REF	
% 9-AD_ALT	10-DP	11-total_cn	    12-mappability_score	13-icgc_donor_id	
% 14-dcc_project_code	15-geneID	16-transcriptID	 17-codon	
% 18-degeneracy	 19-geneName	20-cancer_gene

for i=1:n
    filePattern = fullfile(myFolder, filenames(i));
    fileID = fopen(filePattern);
    C = textscan(fileID,'%f %f %s %s %s %s %s %f %f %f %s %s %s %s %s %s %s %s %s %s','CommentStyle','chr');
    fclose(fileID);

    h=cellfun(@height,C(1));
    C=horzcat(C,C(20));
    for j=1:h
        C{21}(j)=cellstr(tissue(i));
    end

    % Removing noncoding genes
    removeIndex = strcmp(C{20},"NA");
    for j=1:21
        C{j}(removeIndex) = [];
    end
    
    for j=1:21
        C_all{j}=[C_all{j};C{j}];
    end

end

%% count Donors for each tissue

Donor_count=zeros(20,1);
for i=1:20
    index=find(string(C_all{21})==string(tissue_types(i)));
    Donor_count(i)=length(unique(C_all{13}(index)));
end


%%
h=cellfun(@height,C_all(1));
P=zeros(h,8);

% Column1  chromosome
% Column2  position
% Column3  reference ntd
% Column4  alternative ntd
% Column5  degeneracy
% Column6  tissue type (using index of tissue type as a proxy)
% Column7  Cancer gene? (yes 1, no 0)
% Column8  Donor ID
for j=1:h
    P(j,1)=C_all{1}(j);
    P(j,2)=C_all{2}(j);
    P(j,3)=let2num(char(C_all{4}(j)));
    P(j,4)=let2num(char(C_all{5}(j)));
    P(j,5)=str2double(C_all{18}(j));
    Donor=char(C_all{13}(j)); Donor(1:2)=[];
    P(j,8)=str2double(Donor);
    for k=1:20
        if C_all{21}(j)==tissue_types(k)
            P(j,6)=k;
        end
    end
    if C_all{20}(j)=="yes"
        P(j,7)=1;
    end
end

%%
P = sortrows(P,2);
P = sortrows(P,1);

save("Path_To\Data_Extracted_files\PCAWG_Coding_Mutations.mat","P")

% count mutations in each tissue type
counts=zeros(20,1);
for i=1:20
    counts(i)=sum(P(:,6) == i);
end

save("Path_To\Data_Extracted_files\PCAWG_coding_counts","counts")

% Find proportions of each tissue type
prop=zeros(20,1);
for i=1:20
    prop(i)=counts(i)/h;
end

X = categorical(tissue_types);
X = reordercats(X,tissue_types);

figure
bar(X,counts)
title("Number of mutations - coding PCAWG")

figure
bar(X,prop)
title("Proportion of mutations - coding PCAWG")







