% analysis of the combined file: Labeled coding PCAWG mutations
% Lowest level analysis: syn/nonsyn

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

H=load('Path_To\Data_Extracted_files\Labeled_Coding_PCAWG.mat');
H=H.H;

tissue_types=["Bladder","Blood","Bone","Brain","Breast","Cervix","Colorectal","Esophagus","Gall Bladder","Head & neck","Kidney","Liver","Lung","Mesenchymal","Ovary","Pancreas","Prostate","Skin","Stomach","Uterus"];

h=length(H);

% Counting all mutations equally first
Syn_Tv=0;
Syn_Ts=0;
NonSyn_Tv=0;
NonSyn_Ts=0;

for j=1:h
    if (H(j,8)==0 && H(j,9)==0)
        Syn_Tv=Syn_Tv+1;
    elseif (H(j,8)==1 && H(j,9)==0)
        Syn_Ts=Syn_Ts+1;    
    elseif (H(j,8)==0 && H(j,9)==1)
        NonSyn_Tv=NonSyn_Tv+1;
    elseif (H(j,8)==1 && H(j,9)==1)
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

sgtitle("All mutations - Dominated by skin and colon")


%%
% Now counting for each tissue type alone
n=length(tissue_types);

Syn_Tv_tissues=zeros(n,1);
Syn_Ts_tissues=zeros(n,1);
NonSyn_Tv_tissues=zeros(n,1);
NonSyn_Ts_tissues=zeros(n,1);
counts=zeros(n,1); % mutation number of each tissue type

for i=1:n
    counts(i)=sum(H(:,10) == i);
end

save("Path_To\Data_Extracted_files\tissue_mutations_count", "counts")

for j=1:h
    if (H(j,8)==0 && H(j,9)==0)
        Syn_Tv_tissues(H(j,10))=Syn_Tv_tissues(H(j,10))+1;
    elseif (H(j,8)==1 && H(j,9)==0)
        Syn_Ts_tissues(H(j,10))=Syn_Ts_tissues(H(j,10))+1;    
    elseif (H(j,8)==0 && H(j,9)==1)
        NonSyn_Tv_tissues(H(j,10))=NonSyn_Tv_tissues(H(j,10))+1;
    elseif (H(j,8)==1 && H(j,9)==1)
        NonSyn_Ts_tissues(H(j,10))=NonSyn_Ts_tissues(H(j,10))+1;        
    end
end

all_mut_tissues=Syn_Tv_tissues+Syn_Ts_tissues+NonSyn_Tv_tissues+NonSyn_Ts_tissues;
Syn_tissues=Syn_Tv_tissues+Syn_Ts_tissues;
NonSyn_tissues=NonSyn_Tv_tissues+NonSyn_Ts_tissues;
Tv_tissues=Syn_Tv_tissues+NonSyn_Tv_tissues;
Ts_tissues=Syn_Ts_tissues+NonSyn_Ts_tissues;

frac_syn_tissues=Syn_tissues./counts;
frac_nonsyn_tissues=NonSyn_tissues./counts;
frac_tv_tissues=Tv_tissues./counts;
frac_ts_tissues=Ts_tissues./counts;

frac_tv_of_syn_tissues=Syn_Tv_tissues./Syn_tissues;
frac_ts_of_syn_tissues=Syn_Ts_tissues./Syn_tissues;
frac_tv_of_nonsyn_tissues=NonSyn_Tv_tissues./NonSyn_tissues;
frac_ts_of_nonsyn_tissues=NonSyn_Ts_tissues./NonSyn_tissues;

frac_syn_of_tv_tissues=Syn_Tv_tissues./Tv_tissues;
frac_syn_of_ts_tissues=Syn_Ts_tissues./Ts_tissues;
frac_nonsyn_of_tv_tissues=NonSyn_Tv_tissues./Tv_tissues;
frac_nonsyn_of_ts_tissues=NonSyn_Ts_tissues./Ts_tissues;

mean_all_mut=mean(all_mut_tissues);
mean_Syn=mean(Syn_tissues);
mean_NonSyn=mean(NonSyn_tissues);
mean_Tv=mean(Tv_tissues);
mean_Ts=mean(Ts_tissues);

mean_frac_syn=mean(frac_syn_tissues);
mean_frac_nonsyn=mean(frac_nonsyn_tissues);
mean_frac_tv=mean(frac_tv_tissues);
mean_frac_ts=mean(frac_ts_tissues);

mean_frac_tv_of_syn=mean(frac_tv_of_syn_tissues);
mean_frac_ts_of_syn=mean(frac_ts_of_syn_tissues);
mean_frac_tv_of_nonsyn=mean(frac_tv_of_nonsyn_tissues);
mean_frac_ts_of_nonsyn=mean(frac_ts_of_nonsyn_tissues);

mean_frac_syn_of_tv=mean(frac_syn_of_tv_tissues);
mean_frac_syn_of_ts=mean(frac_syn_of_ts_tissues);
mean_frac_nonsyn_of_tv=mean(frac_nonsyn_of_tv_tissues);
mean_frac_nonsyn_of_ts=mean(frac_nonsyn_of_ts_tissues);


median_all_mut=median(all_mut_tissues);
median_Syn=median(Syn_tissues);
median_NonSyn=median(NonSyn_tissues);
median_Tv=median(Tv_tissues);
median_Ts=median(Ts_tissues);

median_frac_syn=median(frac_syn_tissues);
median_frac_nonsyn=median(frac_nonsyn_tissues);
median_frac_tv=median(frac_tv_tissues);
median_frac_ts=median(frac_ts_tissues);

median_frac_tv_of_syn=median(frac_tv_of_syn_tissues);
median_frac_ts_of_syn=median(frac_ts_of_syn_tissues);
median_frac_tv_of_nonsyn=median(frac_tv_of_nonsyn_tissues);
median_frac_ts_of_nonsyn=median(frac_ts_of_nonsyn_tissues);

median_frac_syn_of_tv=median(frac_syn_of_tv_tissues);
median_frac_syn_of_ts=median(frac_syn_of_ts_tissues);
median_frac_nonsyn_of_tv=median(frac_nonsyn_of_tv_tissues);
median_frac_nonsyn_of_ts=median(frac_nonsyn_of_ts_tissues);


figure
subplot(1,2,1)
X = categorical(["All Mutations","Synonymous","Nonsynonymous"]);
X = reordercats(X,["All Mutations","Synonymous","Nonsynonymous"]);
bar(X,[mean_frac_tv,mean_frac_ts;mean_frac_tv_of_syn,mean_frac_ts_of_syn;mean_frac_tv_of_nonsyn,mean_frac_ts_of_nonsyn],'stacked')
title("trv:trs")
legend("Transversions","Transitions")

subplot(1,2,2)
X = categorical(["All Mutations","Transversions","Transitions"]);
X = reordercats(X,["All Mutations","Transversions","Transitions"]);
bar(X,[mean_frac_syn,mean_frac_nonsyn;mean_frac_syn_of_tv,mean_frac_nonsyn_of_tv;mean_frac_syn_of_ts,mean_frac_nonsyn_of_ts],'stacked')
title("syn:nonsyn")
legend("Synonymous","Nonsynonymous")

sgtitle("Tissue types weighted equally - Means")


figure
subplot(1,2,1)
X = categorical(["All Mutations","Synonymous","Nonsynonymous"]);
X = reordercats(X,["All Mutations","Synonymous","Nonsynonymous"]);
bar(X,[median_frac_tv,median_frac_ts;median_frac_tv_of_syn,median_frac_ts_of_syn;median_frac_tv_of_nonsyn,median_frac_ts_of_nonsyn],'stacked')
title("trv:trs")
legend("Transversions","Transitions")

subplot(1,2,2)
X = categorical(["All Mutations","Transversions","Transitions"]);
X = reordercats(X,["All Mutations","Transversions","Transitions"]);
bar(X,[median_frac_syn,median_frac_nonsyn;median_frac_syn_of_tv,median_frac_nonsyn_of_tv;median_frac_syn_of_ts,median_frac_nonsyn_of_ts],'stacked')
title("syn:nonsyn")
legend("Synonymous","Nonsynonymous")

sgtitle("Tissue types weighted equally - Medians")

