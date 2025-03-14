function num = effect2num(effect)
    if effect=="SILENT" 
        num = 0; 
    elseif effect=="MISSENSE" 
        num = 1; 
    elseif effect=="NONSENSE" 
        num = 2; 
    end
end
