%function Pulse_1stVar_leader_autocorr_dynamics_FigS4

path = '~/Dropbox/plxdata/Pulse/Vprobe/prelim_groups/coh_dep_dynamics_based_consec3_p_05_base_adjusted/main/classic_leader/';


filename_a = dir(strcat(path,'30/r/','*_DSP*.mat'));
filename_b = dir(strcat(path,'30/g/','*_DSP*.mat'));
filename_c = dir(strcat(path,'31/r/','*_DSP*.mat'));
filename_d = dir(strcat(path,'31/g/','*_DSP*.mat'));

% decision-related activity before IEM is measured +200 ~ 500 ms after dots
pre_dec1 = 200; post_dec1 = 500;
pre_targ =  0; post_targ = 300; % post-target activity to be controlled in partial corr

moving_binwidth = 300;
moving_step = 20;


sacc1_start = -250;  sacc1_end = 550;
slib_start = -350;  slib_end = 350;
fix2_start = -400;  fix2_end = 200;
dots2_start = -400;  dots2_end = 400;

sacc1_center = sacc1_start:moving_step:sacc1_end;
slib_center = slib_start:moving_step:slib_end;
fix2_center = fix2_start:moving_step:fix2_end;
dots2_center = dots2_start:moving_step:dots2_end;

cellno = 0;

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
        
        clear trials DSP011
        clear targ_times_ti dots_times_ti dots2_times_ti sacc_times_ti slib_times_ti fix2_times_ti choice_sacc_times_ti
        clear t_array_dec_compare_ti t_array_t_ti t_array_s1_ti t_array_sl_ti t_array_f2_ti t_array_t_ti t_array_d2_ti t_array_s_ti
        clear trial_spike_ts Gaddress
        
        filename(i).name
        load([filepath,filename(i).name]);
        
        cellno =  cellno+1;
        trial_spike_ts = zeros(1,1);
        
        trialno = 0;
        
        for j = 1:length(trials)
            
            if ((trials(j).taskid == taskid))&&(trials(j).response >= 0)&&(trials(j).mean_maxFR >= 10) 
                
                if mod(p,2)==0
                    trials(j).t1_dir = trials(j).t1_dir+180;  % coherence sign is based on the leader's preference
                end
                
                trialno = trialno+1;
                dspsize = find(DSP011(j,:) > trials(j).time_sacc(end)+300,1);
                trial_spike_ts(trialno,1:dspsize) = DSP011(j,1:dspsize);
                
                dir_sign = (trials(j).dot_dir ~= trials(j).t1_dir); 
                cho{cellno}(trialno,1) = abs(dir_sign-trials(j).response); % 1 for T+, 0 for T- choice
                coh{cellno}(trialno,1)= trials(j).dot_coh(1)*(1-2*dir_sign);
                
                Gaddress{1}(trialno) = trialno; % 
                
                targ_times_ti(trialno) = trials(j).time_targ_on(1);
                dots_times_ti(trialno) = trials(j).time_dots_on(1);
                sacc_times_ti(trialno) = trials(j).time_sacc(1);
                slib_times_ti(trialno) = trials(j).time_slib_start(1);
                fix2_times_ti(trialno) = trials(j).time_fp_acq(2);
                dots2_times_ti(trialno) = trials(j).time_dots_on(2);
                
                t_array_dec_compare_ti(trialno,:) = [pre_dec1 post_dec1];   % decision-related activity during P1
                t_array_t_ti(trialno,:) = [pre_targ post_targ];             % visual response to the targets (for control)
                
                for k = 1:length(sacc1_center)
                    t_array_s1_ti{k}(trialno,:) = [sacc1_center(k)-moving_binwidth/2 sacc1_center(k)+moving_binwidth/2];
                end
                
                for k = 1:length(slib_center)
                    t_array_sl_ti{k}(trialno,:) = [slib_center(k)-moving_binwidth/2 slib_center(k)+moving_binwidth/2];
                end
                
                for k = 1:length(fix2_center)
                    t_array_f2_ti{k}(trialno,:) = [fix2_center(k)-moving_binwidth/2 fix2_center(k)+moving_binwidth/2];
                end
                
                for k = 1:length(dots2_center)
                    t_array_d2_ti{k}(trialno,:) = [dots2_center(k)-moving_binwidth/2 dots2_center(k)+moving_binwidth/2];
                end
                
            end
        end
        
        clear res_Dec1 std_Dec1 z_Dec1 res_FIX2  std_FIX2 z_FIX2
        clear res_SACC1y std_SACC1y z_SACC1y res_SLIB std_SLIB z_SLIB
        clear res_SACC std_SACC z_SACC res_DOTS2 std_DOTS2 z_DOTS2
        clear res_Targ std_Targ z_Targ
        
        % activity standardized for each cell (1xnTrials column vector)
        [res_Targ, std_Targ] = spk_z2(trial_spike_ts,Gaddress,[targ_times_ti(:)],[pre_targ post_targ],t_array_t_ti);
        z_Targ = res_Targ./std_Targ;
        
        [res_Dec1, std_Dec1] = spk_z2(trial_spike_ts,Gaddress,[dots_times_ti(:)],[pre_dec1 post_dec1],t_array_dec_compare_ti);
        z_Dec1 = res_Dec1./std_Dec1;
        
        % column vectors from individual cells are concatenated 
        if cellno == 1
            z_Dec1_all = z_Dec1;
            z_Targ_all = z_Targ;
            coh_all = coh{cellno};
            cho_all = cho{cellno};
        else
            z_Dec1_all = [z_Dec1_all; z_Dec1];
            z_Targ_all = [z_Targ_all; z_Targ];
            coh_all = [coh_all; coh{cellno}];
            cho_all = [cho_all; cho{cellno}];
        end
        
        % activity aligned onto the first saccade
        for k = 1:length(sacc1_center)
            
            [res_SACC1y{k}, std_SACC1y{k}] = spk_z2(trial_spike_ts,Gaddress,[sacc_times_ti(:)],t_array_s1_ti{k}(1,:),t_array_s1_ti{k});
            z_SACC1y{k} = res_SACC1y{k}./std_SACC1y{k};
            
            if cellno == 1
                z_SACC1y_all{k} = z_SACC1y{k};
            else
                z_SACC1y_all{k} = [z_SACC1y_all{k}; z_SACC1y{k}];
            end
            
        end
        
        % activity aligned onto the pursuit onset
        for k = 1:length(slib_center)
            
            
            [res_SLIB{k}, std_SLIB{k}] = spk_z2(trial_spike_ts,Gaddress,[slib_times_ti(:)],t_array_sl_ti{k}(1,:),t_array_sl_ti{k});
            z_SLIB{k} = res_SLIB{k}./std_SLIB{k};
            
            if cellno == 1
                z_SLIB_all{k} = z_SLIB{k};
            else
                z_SLIB_all{k} = [z_SLIB_all{k}; z_SLIB{k}];
            end
            
        end
        
        % activity aligned onto the resumed fixation
        for k = 1:length(fix2_center)
            
            [res_FIX2{k}, std_FIX2{k}] = spk_z2(trial_spike_ts,Gaddress,[fix2_times_ti(:)],t_array_f2_ti{k}(1,:),t_array_f2_ti{k});
            z_FIX2{k} = res_FIX2{k}./std_FIX2{k};
            
            if cellno == 1
                z_FIX2_all{k} = z_FIX2{k};
            else
                z_FIX2_all{k} = [z_FIX2_all{k}; z_FIX2{k}];
            end
            
        end
        
        % activity aligned onto P2
        for k = 1:length(dots2_center)
            
            [res_DOTS2{k}, std_DOTS2{k}] = spk_z2(trial_spike_ts,Gaddress,[dots2_times_ti(:)],t_array_d2_ti{k}(1,:),t_array_d2_ti{k});
            z_DOTS2{k} = res_DOTS2{k}./std_DOTS2{k};
            
            if cellno == 1
                z_DOTS2_all{k} = z_DOTS2{k};
            else
                z_DOTS2_all{k} = [z_DOTS2_all{k}; z_DOTS2{k}];
            end
            
        end
        
    end % end of cell-specific part
end


% get rid of nan
nanind = [];
for k = 1:length(dots2_center) % removing nan entries
    nanind = [nanind; find(isnan(z_DOTS2_all{k}))];     % for some reason, nan trials found only when aligned on P2
end
nanind = unique(nanind);
length(nanind)

coh_all(nanind) = []; cho_all(nanind) = []; z_Dec1_all(nanind) = []; z_Targ_all(nanind) = [];
for k = 1:length(sacc1_center)
    z_SACC1y_all{k}(nanind) = [];
end
for k = 1:length(slib_center)
    z_SLIB_all{k}(nanind) = [];
end
for k = 1:length(fix2_center)
    z_FIX2_all{k}(nanind) = [];
end
for k = 1:length(dots2_center)
    z_DOTS2_all{k}(nanind) = [];
end

for k = 1:length(sacc1_center)
    [simpleR_all_temp, ~] = partialcorr(z_Dec1_all,z_SACC1y_all{k},ones(size(coh_all))); % simple correlation
    [r_all_temp, p_all_temp] = partialcorr(z_Dec1_all,z_SACC1y_all{k},[coh_all cho_all z_Targ_all]); % partial correlation controlling motion strength, choice, and visual response on each trial 
    r_sacc1_all(k) = r_all_temp;
    p_sacc1_all(k) = p_all_temp;
    simpleR_sacc1(k) = simpleR_all_temp;
end

for k = 1:length(slib_center)
    [simpleR_all_temp, ~] = partialcorr(z_Dec1_all,z_SLIB_all{k},ones(size(coh_all)));
    [r_all_temp, p_all_temp] = partialcorr(z_Dec1_all,z_SLIB_all{k},[coh_all cho_all z_Targ_all]);
    r_slib_all(k) = r_all_temp;
    p_slib_all(k) = p_all_temp;
    simpleR_slib(k) = simpleR_all_temp;
    
end

for k = 1:length(fix2_center)
    [simpleR_all_temp, ~] = partialcorr(z_Dec1_all,z_FIX2_all{k},ones(size(coh_all)));
    [r_all_temp, p_all_temp] = partialcorr(z_Dec1_all,z_FIX2_all{k},[coh_all cho_all z_Targ_all]);
    r_fix2_all(k) = r_all_temp;
    p_fix2_all(k) = p_all_temp;
    simpleR_fix2(k) = simpleR_all_temp;
end

for k = 1:length(dots2_center)
    [simpleR_all_temp, ~] = partialcorr(z_Dec1_all,z_DOTS2_all{k},ones(size(coh_all)));
    [r_all_temp, p_all_temp] = partialcorr(z_Dec1_all,z_DOTS2_all{k},[coh_all cho_all z_Targ_all]);
    r_dots2_all(k) = r_all_temp;
    p_dots2_all(k) = p_all_temp;
    simpleR_dots2(k) = simpleR_all_temp;
    
end


figure(84)

subplot(141)
plot(sacc1_center,simpleR_sacc1,sacc1_center,r_sacc1_all,'LineWidth',1.5)
line([0 0],[0 1],'Color','k');
xlabel('Sacc to T^0')
xlim([sacc1_center(1) sacc1_center(end)]); ylim([0.12 0.29]);
p_max_sacc1 = max(p_sacc1_all)
set(gca, 'TickDir','out');
box off

subplot(142)
plot(slib_center,simpleR_slib,slib_center,r_slib_all,'LineWidth',1.5)
line([0 0],[0 0.4],'Color','k');
xlabel('Pursuit on')
xlim([slib_center(1) slib_center(end)]); ylim([0.13 0.35]);
p_max_slib = max(p_slib_all)
set(gca, 'TickDir','out');
box off

subplot(143)
plot(fix2_center,simpleR_fix2,fix2_center,r_fix2_all,'LineWidth',1.5)
line([0 0],[0 0.4],'Color','k');
xlabel('Resumed fixation')
xlim([fix2_center(1) fix2_center(end)]); ylim([0.13 0.35]);
p_max_fix2 = max(p_fix2_all)
set(gca, 'TickDir','out');
box off

subplot(144)
plot(dots2_center,simpleR_dots2,dots2_center,r_dots2_all,'LineWidth',1.5)
line([0 0],[0 0.4],'Color','k');
xlabel('P2 onset')
xlim([dots2_center(1) dots2_center(end)]); ylim([0.13 0.35]);
p_max_d2 = max(p_dots2_all)
set(gca, 'TickDir','out');
box off
legend('Simple r','partial corr r')

suptitle('Autocorrelation of leader neuron activity, controlling P1 coh, choice, & visual response to the targets');


