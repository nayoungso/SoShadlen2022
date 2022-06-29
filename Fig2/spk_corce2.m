function [res_cnt, mean_cnt, trialid] = spk_corce2(spk, addresses, align_time, time_array, actual_time_array, binwidth, id_start)

%% for CorCE computation, returns the residual spk count, mean count, and trialid (for concatenation)


totalno = 0;
for i = 1:length(addresses)   % total number of trials across different conditions
    totalno = totalno+length(addresses{i});
end


% constructing cnt_vector (determines if a certain time bin needs to be excluded for the analysis, e.g., attrition rule for the decision-related activity )
trialnocnt = 0;
startid = 0;
cnt_array = time_array(1)-binwidth/2:binwidth:time_array(2)-binwidth/2;

for k = 1:length(addresses)
    
    cnt_vector = nan(length(addresses{k}),length(cnt_array));  % default for cnt_vector is nan (not including the time bin)

    for i = 1:length(addresses{k})
        
        trialnocnt = trialnocnt+1;
        
        for j = 1:size(spk,2)
            
            if spk(addresses{k}(i),j) > 0
                
                spk_ts(i,j) = spk(addresses{k}(i),j)-align_time(addresses{k}(i));
                
            else
                spk_ts(i,j) = nan;
                
            end
            
        end
        
        for j = 1:length(cnt_array)
            
            spk_ts_cnt(i,j) = histcounts(spk_ts(i,:),[cnt_array(j) cnt_array(j)+binwidth]);
            
            if (cnt_array(j)+binwidth/2 >= actual_time_array(addresses{k}(i),1))&(cnt_array(j)+binwidth/2 <=actual_time_array(addresses{k}(i),2))
                
                cnt_vector(i,j) = 1;    % cnt_vector is flagged such that the time bin is included in the analysis
            end
        end
    
        spk_ts_cnt(i,:) = spk_ts_cnt(i,:).*cnt_vector(i,:);
        
    end
    
    try
    mean_spk_cnt = nanmean(spk_ts_cnt,1);  % mean spk count for the given condition k
    catch
        mean_spk_cnt = nan(1,length(cnt_array));
    end
    
    for i = 1:length(addresses{k})
        raw_cnt(startid+i,:) = spk_ts_cnt(i,:);
        res_cnt(startid+i,:) = spk_ts_cnt(i,:)-mean_spk_cnt;
        mean_cnt(startid+i,:) = mean_spk_cnt;
        trialid(startid+i,1) = id_start+k;
    end
    
    startid = startid+length(addresses{k});
    
    clear spk_ts_cnt_temp spk_ts_cnt spk_ts mean_spk_cnt   
    
end



return

end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


