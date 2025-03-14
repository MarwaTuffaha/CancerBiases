
%clear


fileID = fopen("Path_To\Germline Dataset\germinal_ultimate_dataset.txt");
% # rows = 679,547
% Column1  chromosome
% Column2  start
% Column3  stop
% Column4  reference ntd
% Column5  alternative ntd
% Column6  study            ignore for now
% Column7  mutation type    ignore for now
% Column8  cases types      ignore for now


C = textscan(fileID,'%s %*f %f %s %s %*s %*s %*s'); 

fclose(fileID);

h=cellfun(@height,C(1));

%%
% Create a large matrix G and save the data numerically to save space
% A-1   C-2   G-3   T-4

% Column1  chromosome
% Column2  position
% Column3  reference ntd
% Column4  alternative ntd

G=zeros(h,4);

for j=1:h
    a=char(C{1}(j));
    G(j,1)=str2double(a(4:end));
    G(j,2)=C{2}(j);
    G(j,3)=let2num(char(C{3}(j)));
    G(j,4)=let2num(char(C{4}(j)));
end

%%

G = sortrows(G,2);
G = sortrows(G,1);

%%

save("Path_To\Coding_Analysis\Data_Extracted_files\Germline_Mutations_Numbers",'G', '-v7.3')



%%
