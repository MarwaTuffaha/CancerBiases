% analysis of the combined file on the 3mer level: Labeled coding PCAWG mutations

% Column1  chromosome
% Column2  position
% Column3  reference ntd
% Column4  alternative ntd
% Column5  reference 3-mer
% Column6  mutation effect
% Column7  0:trv | 1:trs
% Column8  0:syn | 1:nonsyn
% Column9 Cancer gene? (yes 1, no 0)

H=load('Path_To\Data_Extracted_files\Labeled_coding_Germline.mat');
H=H.H;
h=length(H);

tissue_types=["Bladder","Blood","Bone","Brain","Breast","Cervix","Colorectal","Esophagus","Gall Bladder","Head & neck","Kidney","Liver","Lung","Mesenchymal","Ovary","Pancreas","Prostate","Skin","Stomach","Uterus"];

mers=load('3mers.mat');
U=mers.U;
n=length(U);

% Counting ts/tv - syn/nonsyn mutations for each 3mer mutation
G_syn=zeros(n,1);
G_nonsyn=zeros(n,1);
G_tv=zeros(n,1);
G_ts=zeros(n,1);

G_syn_c=zeros(n,1); % syn muts in cancer genes
G_nonsyn_c=zeros(n,1); % nonsyn muts in cancer genes
G_syn_nc=zeros(n,1); % syn muts in noncancer genes
G_nonsyn_nc=zeros(n,1); % nonsyn muts in noncancer gene

for i=1:h
    for j=1:n
        if H(i,5)==U(j,1) && H(i,4)==U(j,2)
            if H(i,7)==0
                G_tv(j)=G_tv(j)+1;
            else
                G_ts(j)=G_ts(j)+1;
            end
            if H(i,8)==0
                G_syn(j)=G_syn(j)+1;
                if H(i,9)==1
                    G_syn_c(j)=G_syn_c(j)+1;
                else
                    G_syn_nc(j)=G_syn_nc(j)+1;
                end
            else
                G_nonsyn(j)=G_nonsyn(j)+1;
                if H(i,9)==1
                    G_nonsyn_c(j)=G_nonsyn_c(j)+1;
                else
                    G_nonsyn_nc(j)=G_nonsyn_nc(j)+1;
                end
            end            
        end
    end
end

save('Path_To\Data_Extracted_files\Germline_synNonsyn','G_syn','G_nonsyn')
save('Path_To\Data_Extracted_files\Germline_synNonsyn_CancerGenes','G_syn_c','G_nonsyn_c','G_syn_nc','G_nonsyn_nc')

%% Find 1mer spectra

ind1=1:16; ind2=17:32; ind3=33:48; ind4=49:64; ind5=65:80; ind6=81:96;

onemer_germline=zeros(1,6); onemer_c_germline=zeros(1,6); onemer_nc_germline=zeros(1,6);
onemer_germline(1)=(sum(G_syn(ind1))+sum(G_nonsyn(ind1)));
onemer_germline(2)=(sum(G_syn(ind2))+sum(G_nonsyn(ind2)));
onemer_germline(3)=(sum(G_syn(ind3))+sum(G_nonsyn(ind3)));
onemer_germline(4)=(sum(G_syn(ind4))+sum(G_nonsyn(ind4)));
onemer_germline(5)=(sum(G_syn(ind5))+sum(G_nonsyn(ind5)));
onemer_germline(6)=(sum(G_syn(ind6))+sum(G_nonsyn(ind6)));

onemer_germline_c(1)=(sum(G_syn_c(ind1))+sum(G_nonsyn_c(ind1)));
onemer_germline_c(2)=(sum(G_syn_c(ind2))+sum(G_nonsyn_c(ind2)));
onemer_germline_c(3)=(sum(G_syn_c(ind3))+sum(G_nonsyn_c(ind3)));
onemer_germline_c(4)=(sum(G_syn_c(ind4))+sum(G_nonsyn_c(ind4)));
onemer_germline_c(5)=(sum(G_syn_c(ind5))+sum(G_nonsyn_c(ind5)));
onemer_germline_c(6)=(sum(G_syn_c(ind6))+sum(G_nonsyn_c(ind6)));

onemer_germline_nc(1)=(sum(G_syn_nc(ind1))+sum(G_nonsyn_nc(ind1)));
onemer_germline_nc(2)=(sum(G_syn_nc(ind2))+sum(G_nonsyn_nc(ind2)));
onemer_germline_nc(3)=(sum(G_syn_nc(ind3))+sum(G_nonsyn_nc(ind3)));
onemer_germline_nc(4)=(sum(G_syn_nc(ind4))+sum(G_nonsyn_nc(ind4)));
onemer_germline_nc(5)=(sum(G_syn_nc(ind5))+sum(G_nonsyn_nc(ind5)));
onemer_germline_nc(6)=(sum(G_syn_nc(ind6))+sum(G_nonsyn_nc(ind6)));

Freq=zeros(1,6);
Freq(1)=onemer_germline(1)/h;
Freq(2)=onemer_germline(2)/h;
Freq(3)=onemer_germline(3)/h;
Freq(4)=onemer_germline(4)/h;
Freq(5)=onemer_germline(5)/h;
Freq(6)=onemer_germline(6)/h;

SE=zeros(1,6);
for i=1:6
    SE(i)=sqrt(Freq(i)*(1-Freq(i))/h);
end

save("Path_To\Data_Extracted_files\Germline_1mer_coding_freq","Freq")

% Find TsTv overall
Ts_count_germline=(sum(G_syn(ind3))+sum(G_nonsyn(ind3))+sum(G_syn(ind5))+sum(G_nonsyn(ind5)));
Tv_count_germline=h-Ts_count_germline;
Ts_freq_germline=Ts_count_germline/h;

Ts_count_germline_c=(sum(G_syn_c(ind3))+sum(G_nonsyn_c(ind3))+sum(G_syn_c(ind5))+sum(G_nonsyn_c(ind5)));
Tv_count_germline_c=sum(G_syn_c)+sum(G_nonsyn_c)-Ts_count_germline_c;

Ts_count_germline_nc=(sum(G_syn_nc(ind3))+sum(G_nonsyn_nc(ind3))+sum(G_syn_nc(ind5))+sum(G_nonsyn_nc(ind5)));
Tv_count_germline_nc=sum(G_syn_nc)+sum(G_nonsyn_nc)-Ts_count_germline_nc;

save("Path_To\Data_Extracted_files\Germline_counts_TsTv1mer_coding","Ts_count_germline","Tv_count_germline","onemer_germline")

save("Path_To\Data_Extracted_files\Germline_counts_TsTv1mer_CancerGenes","Ts_count_germline_c","Tv_count_germline_c","onemer_germline_c","Ts_count_germline_nc","Tv_count_germline_nc","onemer_germline_nc")


