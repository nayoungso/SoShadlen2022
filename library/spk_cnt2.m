function [spk_ts_cnt, fr] = spk_cnt2(spk, addresses, align_time, time_array, actual_time_array)

% note that the firing rate (fr) is calculated in units of spikes/ms (not spikes/s)

for i = 1:length(addresses)
    
    for j = 1:size(spk,2)
        
        if spk(addresses(i),j) > 0
            
            spk_ts(i,j) = spk(addresses(i),j)-align_time(addresses(i));
        
        else
            spk_ts(i,j) = nan;
        
        end
        
    end
    
end

for i = 1:length(addresses)
    
    spk_ts_cnt(i) = sum(histc(spk_ts(i,:),actual_time_array(addresses(i),:)));
    duration = actual_time_array(addresses(i),2)-actual_time_array(addresses(i),1);
    fr(i) = spk_ts_cnt(i)./duration;
    
end


return

end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


