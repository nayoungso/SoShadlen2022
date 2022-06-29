function Pulse_1stVar_deltaFR_P1P2_hist_FigS4


path = '~/Dropbox/plxdata/Pulse/Vprobe/prelim_groups/coh_dep_dynamics_based_consec3_p_05_base_adjusted/main/classic_leader/';


% delta FR 
binwidth = 100; 
pre_dots = 100; 
post_dots = 400; 

outlier_remv = 1; 

% criteria for deltaFR_{other}
z_thr_array = [1 1.5 2];

filename_a = dir(strcat(path,'30/r/','*_DSP*.mat'));
filename_b = dir(strcat(path,'30/g/','*_DSP*.mat'));
filename_c = dir(strcat(path,'31/r/','*_DSP*.mat'));
filename_d = dir(strcat(path,'31/g/','*_DSP*.mat'));

trialno_total = 0; 
trial_spike_ts = zeros(1,1);

coh_set = [-640 -320 -160 -80 -40 -1 1 40 80 160 320 640];

for p = 1:4
    
    clear filename
    
    if p==1
        taskid = 30; filename = filename_a; filepath = [path,'/30/r/'];
    elseif p==2
        taskid = 30; filename = filename_b; filepath = [path,'/30/g/'];
    elseif p==3
        taskid = 31; filename = filename_c; filepath = [path,'/31/r/'];
    elseif p==4
        taskid = 31; filename = filename_d; filepath = [path,'/31/g/'];
    end
    
    for i = 1:length(filename)
        
        clear trials trial_spike_ts dots_time1 dots_time2 D1_spks_start D2_spks_start dir_sign 
        clear D1_spks_end D2_spks_end 
        clear D1_raw_spks D2_raw_spks 
        clear t_array_cnt G1address dir_sign choice dirsign
        filename(i).name
        load([filepath,filename(i).name]);
        
        trialno_cell = 0;
        
        for j = 1:length(trials)
            
            if (trials(j).taskid == taskid)&&(trials(j).response>=0)&&(trials(j).mean_maxFR >= 10)
                
                if mod(p,2)==0
                    trials(j).t1_dir = trials(j).t1_dir+180;  % coherence sign is based on the leader's preference
                end
                
                trialno_cell = trialno_cell+1;
                
                dspsize = find(DSP011(j,:) > trials(j).time_sacc(end)+300,1);                
                trial_spike_ts(trialno_cell,1:dspsize) = DSP011(j,1:dspsize);
                
                dots_time1(trialno_cell) = trials(j).time_dots_on(1);
                dots_time2(trialno_cell) = trials(j).time_dots_on(2);
                
                
                G1address{1}(trialno_cell) = trialno_cell;
                dir_sign(trialno_cell) = (trials(j).dot_dir ~= trials(j).t1_dir); % 0 for t1_dir / 1 for t2_dir
                
                t_array_cnt(trialno_cell,:) = [-binwidth/2 binwidth/2];
                
                if ((dir_sign(trialno_cell) == 0)&&(trials(j).response == 1))||((dir_sign(trialno_cell) == 1)&&(trials(j).response == 0))
                    choice(trialno_cell) = 1;   % t1 choice
                else
                    choice(trialno_cell) = 0;   % t2 choice
                end
                
                
            end
            
            
        end
        
        if trialno_cell > 3
        [D1_spks_start, ~] = spk_cnt2(trial_spike_ts,G1address{1},[dots_time1(:)+pre_dots],[-binwidth/2 binwidth/2],t_array_cnt);
        [D1_spks_end, ~] = spk_cnt2(trial_spike_ts,G1address{1},[dots_time1(:)+post_dots],[-binwidth/2 binwidth/2],t_array_cnt);
        D1_spks_start_det = nanmean(D1_spks_start);
        D1_spks_end_det = nanmean(D1_spks_end);
        
        [D2_spks_start, ~] = spk_cnt2(trial_spike_ts,G1address{1},[dots_time2(:)+pre_dots],[-binwidth/2 binwidth/2],t_array_cnt);
        [D2_spks_end, ~] = spk_cnt2(trial_spike_ts,G1address{1},[dots_time2(:)+post_dots],[-binwidth/2 binwidth/2],t_array_cnt);
        D2_spks_start_det = nanmean(D2_spks_start);
        D2_spks_end_det = nanmean(D2_spks_end);
        
        D1_raw_spks = (D1_spks_end -D1_spks_end_det) - (D1_spks_start - D1_spks_start_det);
        D2_raw_spks = (D2_spks_end -D2_spks_end_det) - (D2_spks_start - D2_spks_start_det);
        
        D1_start_det_temp = (D1_spks_start - D1_spks_start_det).*(1000/binwidth); 
        std_D1_start_det = nanstd(D1_start_det_temp);
        
        D2_start_det_temp = (D2_spks_start - D2_spks_start_det).*(1000/binwidth); 
        std_D2_start_det = nanstd(D2_start_det_temp);
        
        D1_delta_spks(trialno_total+1:trialno_total+trialno_cell) = (D1_raw_spks).*(1000/binwidth);
        z_D1_delta_spks(trialno_total+1:trialno_total+trialno_cell) = (D1_delta_spks(trialno_total+1:trialno_total+trialno_cell))/std_D1_start_det;
        
        D2_delta_spks(trialno_total+1:trialno_total+trialno_cell) = (D2_raw_spks).*(1000/binwidth);
        z_D2_delta_spks(trialno_total+1:trialno_total+trialno_cell) = (D2_delta_spks(trialno_total+1:trialno_total+trialno_cell))/std_D2_start_det;
        
        choice_grand(trialno_total+1:trialno_total+trialno_cell) = choice;
        trialno_total = trialno_total+trialno_cell;
        
        end
        
    end
    
    
end % end of p-loop


% handling outliers 
outind = [];
crit_z = 5;
if outlier_remv == 1    
    outind = [outind, find(z_D1_delta_spks > crit_z)];
    outind = [outind, find(z_D2_delta_spks > crit_z)];
    outind = [outind, find(z_D1_delta_spks < -crit_z)];
    outind = [outind, find(z_D2_delta_spks < -crit_z)];
end

% nan removal 
temp_nandelta = isnan(z_D2_delta_spks);
temp_infdelta = isinf(z_D2_delta_spks);
outind = [outind, find(temp_nandelta==1)];
outind = [outind, find(temp_infdelta==1)];

choice_grand(outind) = []; 
D1_delta_spks(outind) = []; D2_delta_spks(outind) = [];
z_D1_delta_spks(outind) = []; z_D2_delta_spks(outind) = []; 

for i = 1:length(z_D1_delta_spks)
    
    for_z_D1_delta_spks(i) = z_D1_delta_spks(i);
    for_z_D2_delta_spks(i) = z_D2_delta_spks(i);
    
    if choice_grand(i) == 0 % t2 choice   %dir_sign_grand(i) == 1 % t2_dir %choice_grand(i) == 0 % T2 choice
        
        for_z_D1_delta_spks(i) = -z_D1_delta_spks(i);
        for_z_D2_delta_spks(i) = -z_D2_delta_spks(i);
    end
end


figure(6); clf;
for i = 1:length(z_thr_array)
    
    clear idxall1 idxall2 deltaFR12_array
    idxall1 = find(for_z_D2_delta_spks > z_thr_array(i));
    idxall2 = find(for_z_D1_delta_spks > z_thr_array(i));
    
    deltaFR12_array = for_z_D1_delta_spks(idxall1);
    deltaFR12_array = [deltaFR12_array for_z_D2_delta_spks(idxall2)];
    
    subplot(length(z_thr_array),1,i);
    
    histogram(deltaFR12_array,'NumBins',60,'Normalization','probability');
    title(['mean deltaFR when (Z(deltaFR_{other}) > ',num2str(z_thr_array(i)),' )  = ',num2str(nanmean(deltaFR12_array))],'fontsize',14)
    line([0 0],[0 0.15],'color','k');
    xlim([-5 5])
    set(gca,'TickDir','out','FontSize',14); box off;
    
    [zh_delta12_array_2tail(i) zp_delta12_array_2tail(i)] = ttest(deltaFR12_array);
end


end
