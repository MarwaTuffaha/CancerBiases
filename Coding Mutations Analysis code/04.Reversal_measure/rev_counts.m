function d = rev_counts(spect,germ_spect,unbiased_spect) 
% takes a spectrum and returns its bias reversal measure from germline 
% towards the uniform spectrum

d=0;
for i=1:length(germ_spect)
    s=spect(i); g=germ_spect(i); u=unbiased_spect(i);
    if u<g
        if s<g
            add=abs(s-g);
        else % s>=g
            add=-abs(s-g);
        end
    else % u>=g
        if s>g
            add=abs(s-g);
        else % s<=g
            add=-abs(s-g);
        end
    end
    d=d+add;
end

end  


