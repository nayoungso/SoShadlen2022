function PSTH = spk_psth_v2(spk, addresses, align_time, time_array, actual_time_array, binwidth, step)


cnt_array = (time_array(1)-binwidth/2):step:(time_array(2)-binwidth/2); % start time of the bins 
cnt_vector = nan(length(addresses),length(cnt_array));  % units of bin

for i = 1:length(addresses)
    
    
    %align the spike time stamps w.r.t the align_time
    for j = 1:size(spk,2)
        
        if spk(addresses(i),j) > 0
            
            spk_ts(i,j) = spk(addresses(i),j)-align_time(addresses(i));
        
        else
            spk_ts(i,j) = nan;
        
        end
        
    end
    
    
    % constructing cnt vector (to apply the attrition rule for variable duration dots), and counting spikes
    for j = 1:length(cnt_array)
        
        spk_ts_cnt(i,j) = histcounts(spk_ts(i,:),[cnt_array(j) cnt_array(j)+binwidth]);
        
        % actual_time array is to handle the trials of different lengths 
        % (e.g., attrition rule for the motion-viewing epochs in the variable duration task)
        if (cnt_array(j)+binwidth >= actual_time_array(addresses(i),1))&(cnt_array(j) <=actual_time_array(addresses(i),2))   
            cnt_vector(i,j) = 1;
        end
    end
    
    spk_ts_cnt(i,:) = spk_ts_cnt(i,:).*cnt_vector(i,:);
            
end


try
PSTH_temp = nanmean(spk_ts_cnt,1);
catch
    PSTH_temp = nan;
end

PSTH = PSTH_temp.*(1000/binwidth);


return

end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


