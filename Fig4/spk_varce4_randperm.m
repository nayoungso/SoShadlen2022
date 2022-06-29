function [res_cnt, mean_cnt, trialid, raw_cnt] = spk_varce4_randperm(spk, addresses, align_time, time_array, id_start, permno)

%% same as spk_varce4, but spk_varce4_randperm generates shuffled surrogates 

trialnocnt = 0;
startid = 0;
for k = 1:length(addresses) % for each condition (coherence)
    
    
    for i = 1:length(addresses{k})
        
        trialnocnt = trialnocnt+1;
        
        for j = 1:size(spk,2)
            
            if spk(addresses{k}(i),j) > 0
                
                spk_ts(i,j) = spk(addresses{k}(i),j)-align_time(addresses{k}(i));
                
            else
                spk_ts(i,j) = nan;
                
            end
            
        end
        
        spk_ts_cnt(i,1) = histcounts(spk_ts(i,:),time_array);
        
          
    end

    try
    mean_spk_cnt = nanmean(spk_ts_cnt,1);  % mean spk count for the given condition k
    catch
        mean_spk_cnt = nan;
    end
    
    for j = 1:permno
        
        clear tempid
        tempid = randperm(length(addresses{k}));
        
        for i = 1:length(addresses{k})
            raw_cnt(startid+i,j) = spk_ts_cnt(tempid(i),1);
            res_cnt(startid+i,j) = spk_ts_cnt(tempid(i),1)-mean_spk_cnt;
            mean_cnt(startid+i,j) = mean_spk_cnt;
            trialid(startid+i,j) = id_start+k;
        end
    end
    
    startid = startid+length(addresses{k});
    
    clear spk_ts_cnt_temp spk_ts_cnt spk_ts mean_spk_cnt   
    
end


return

end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


