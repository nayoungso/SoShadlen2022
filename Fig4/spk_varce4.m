function [res_cnt, mean_cnt, trialid, raw_cnt] = spk_varce4(spk, addresses, align_time, time_array, id_start)

%% for VarCE computation, returns the residual spk count, mean count, and the trialid
%% similar to spk_varce3, but spk_varce4 returns activity from one epoch (spk_varce3 returns an array of activity measured from different time windows) 


trialnocnt = 0;
startid = 0;
for k = 1:length(addresses)  % for each motion coherence
    
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
    
    for i = 1:length(addresses{k})
        raw_cnt(startid+i,1) = spk_ts_cnt(i,1);
        res_cnt(startid+i,1) = spk_ts_cnt(i,1)-mean_spk_cnt;    % residual spk counts 
        mean_cnt(startid+i,1) = mean_spk_cnt;                   % mean spk counts (trials of the same signed motion coherence)
        trialid(startid+i,1) = id_start+k;                      % to concatenate the data from different coherences 
    end
    
    startid = startid+length(addresses{k});
    
    clear spk_ts_cnt_temp spk_ts_cnt spk_ts mean_spk_cnt   
    
end

return

end
    
    
