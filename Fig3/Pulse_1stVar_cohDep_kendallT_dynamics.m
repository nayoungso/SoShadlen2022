%function Pulse_1stVar_cohDep_kendallT_dynamics


path = ['~/Dropbox/plxdata/Pulse/Vprobe/AllCells/'];
filename = dir(strcat(path,'*DSP*.mat'))

binwidth = 300;
stepsize = 50;

CohDep_30{1,1} = 'Filename'; CohDep_30{1,2} = 'P1'; CohDep_30{1,3} = 'P2-coh2'; CohDep_30{1,4} = 'P2-sumcoh';
CohDep_30{1,5} = 'sacc1'; CohDep_30{1,6} = 'pursuit'; CohDep_30{1,7} = 'FP2';  CohDep_30{1,8} = 'choice'; 

CohDep_31{1,1} = 'Filename'; CohDep_31{1,2} = 'P1'; CohDep_31{1,3} = 'P2-coh2'; CohDep_31{1,4} = 'P2-sumcoh';
CohDep_31{1,5} = 'sacc1'; CohDep_31{1,6} = 'pursuit'; CohDep_31{1,7} = 'FP2'; CohDep_31{1,8} = 'choice'; %CohDep_31{1,9} = 'P2-coh1';

p_30{1,1} = 'Filename'; p_30{1,2} = 'P1'; p_30{1,3} = 'P2-coh2'; p_30{1,4} = 'P2-sumcoh';
p_30{1,5} = 'sacc1'; p_30{1,6} = 'pursuit'; p_30{1,7} = 'FP2';  p_30{1,8} = 'choice';  %p_30{1,9} = 'P2-coh1';

p_31{1,1} = 'Filename'; p_31{1,2} = 'P1'; p_31{1,3} = 'P2-coh2'; p_31{1,4} = 'P2-sumcoh';
p_31{1,5} = 'sacc1'; p_31{1,6} = 'pursuit'; p_31{1,7} = 'FP2'; p_31{1,8} = 'choice'; %p_31{1,9} = 'P2-coh1';

base_time = [-300 0]; % activity before dots onset will be subtracted (a device for correcting drift in the multi-channel recording)

% time windows of the survey -- note this 'survey' windows are wider than the 'criteria' windows for cell categorization 
P1_time = [-50 500]; sacc1_time = [-400 400]; pursuit_time = [-400 400];
FP2_time = [-400 400]; P2_time = [-400 400]; choice_time = [-400 200];

p1_t_start = P1_time(1)-binwidth/2:stepsize:P1_time(end)-binwidth/2;  
p2_t_start = P2_time(1)-binwidth/2:stepsize:P2_time(end)-binwidth/2;
sacc1_t_start = sacc1_time(1)-binwidth/2:stepsize:sacc1_time(end)-binwidth/2;
pursuit_t_start = pursuit_time(1)-binwidth/2:stepsize:pursuit_time(end)-binwidth/2;
FP2_t_start = FP2_time(1)-binwidth/2:stepsize:FP2_time(end)-binwidth/2;
choice_t_start = choice_time(1)-binwidth/2:stepsize:choice_time(end)-binwidth/2;

TimeArray{1,1} = 'P1 time'; TimeArray{1,2} = p1_t_start; 
TimeArray{2,1} = 'Sacc1 time'; TimeArray{2,2} = sacc1_t_start;
TimeArray{3,1} = 'Pursuit time'; TimeArray{3,2} = pursuit_t_start;
TimeArray{4,1} = 'FP2 time'; TimeArray{4,2} = FP2_t_start;
TimeArray{5,1} = 'P2 time'; TimeArray{5,2} = p2_t_start;
TimeArray{6,1} = 'Choice time'; TimeArray{6,2} = choice_t_start;
TimeArray{7,1} = 'binwidth'; TimeArray{7,2} = binwidth;
TimeArray{8,1} = 'stepsize'; TimeArray{8,2} = stepsize;

for ii = 1:length(filename)
    clear trials
    
    clear coh_30 sum_coh_30 P1_30 P2_30 sacc1_30 pursuit_30 FP2_30 choice_30
    clear coh_31 sum_coh_31 P1_31 P2_31 sacc1_31 pursuit_31 FP2_31 choice_31
    
    
    fullfilename = filename(ii).name;
    load([path,filename(ii).name]);
    
    CohDep_30{ii+1,1} = fullfilename;   
    CohDep_31{ii+1,1} = fullfilename;   
    
    p_30{ii+1,1} = fullfilename;
    p_31{ii+1,1} = fullfilename;
    
    dspsize = size(DSP011,2);
    
    trialno_30 = 0;
    trialno_31 = 0;
    
    for j = 1:length(trials)
        
        clear spktimes spktimes_base spktimes_P1 spktimes_P2 spktimes_sacc1 spktimes_pursuit spktimes_FP2 spktimes_choice
        
        % Target configuration A
        if (trials(j).taskid==30)&&(trials(j).response>=0)&&(trials(j).mean_maxFR >= 10)    
            
            trialno_30 = trialno_30+1;
            spktimes = DSP011(j,1:dspsize);
            
            % spike time stamps re-aligned onto the event of interest
            spktimes_base = spktimes - trials(j).time_dots_on(1);            
            spktimes_P1 = spktimes - trials(j).time_dots_on(1);
            spktimes_P2 = spktimes - trials(j).time_dots_on(2);
            spktimes_sacc1 = spktimes - trials(j).time_sacc(1);
            spktimes_pursuit = spktimes - trials(j).time_slib_start;
            spktimes_FP2 = spktimes - trials(j).time_fp_acq(2);
            spktimes_choice = spktimes - trials(j).time_sacc(end);
            
            pre_dots_spk = histcounts(spktimes_base,base_time);
            
            dir_sign = (trials(j).dot_dir ~= trials(j).t1_dir); 
            coh_30(trialno_30,:) = trials(j).dot_coh(:)*(1-2*dir_sign); % signed coh ([C1 C2])
            sum_coh_30(trialno_30,1) = sum(coh_30(trialno_30,:));       % summed coh (regressor for the activity during P2)
            
            for i = 1:length(p1_t_start)
                P1_30(trialno_30,i) = histcounts(spktimes_P1,[p1_t_start(i) p1_t_start(i)+binwidth]) - pre_dots_spk;
            end
            
            for i = 1:length(p2_t_start)
                P2_30(trialno_30,i) = histcounts(spktimes_P2,[p2_t_start(i) p2_t_start(i)+binwidth]) - pre_dots_spk;
            end
            
            for i = 1:length(sacc1_t_start)
                sacc1_30(trialno_30,i) = histcounts(spktimes_sacc1,[sacc1_t_start(i) sacc1_t_start(i)+binwidth]) - pre_dots_spk;
            end
            for i = 1:length(pursuit_t_start)
                pursuit_30(trialno_30,i) = histcounts(spktimes_pursuit,[pursuit_t_start(i) pursuit_t_start(i)+binwidth]) - pre_dots_spk;
            end
            for i = 1:length(FP2_t_start)
                FP2_30(trialno_30,i) = histcounts(spktimes_FP2,[FP2_t_start(i) FP2_t_start(i)+binwidth]) - pre_dots_spk;
            end
            
            for i = 1:length(choice_t_start)
                choice_30(trialno_30,i) = histcounts(spktimes_choice,[choice_t_start(i) choice_t_start(i)+binwidth]) - pre_dots_spk;
            end
            
        %Target configuration B
        elseif (trials(j).taskid==31)&&(trials(j).response>=0)&&(trials(j).mean_maxFR >= 10)    
            
            trialno_31 = trialno_31+1;
            spktimes = DSP011(j,1:dspsize);
            
            % spike time stamps re-aligned onto the event of interest
            spktimes_base = spktimes - trials(j).time_dots_on(1);
            spktimes_P1 = spktimes - trials(j).time_dots_on(1);
            spktimes_P2 = spktimes - trials(j).time_dots_on(2);
            spktimes_sacc1 = spktimes - trials(j).time_sacc(1);
            spktimes_pursuit = spktimes - trials(j).time_slib_start;
            spktimes_FP2 = spktimes - trials(j).time_fp_acq(2);
            spktimes_choice = spktimes - trials(j).time_sacc(end);
            
            pre_dots_spk = histcounts(spktimes_base,base_time);
            
            dir_sign = (trials(j).dot_dir ~= trials(j).t1_dir); 
            coh_31(trialno_31,:) = trials(j).dot_coh(:)*(1-2*dir_sign);
            sum_coh_31(trialno_31,1) = sum(coh_31(trialno_31,:));
            
            for i = 1:length(p1_t_start)
                P1_31(trialno_31,i) = histcounts(spktimes_P1,[p1_t_start(i) p1_t_start(i)+binwidth]) - pre_dots_spk;
            end
            
            for i = 1:length(p2_t_start)
                P2_31(trialno_31,i) = histcounts(spktimes_P2,[p2_t_start(i) p2_t_start(i)+binwidth]) - pre_dots_spk;
            end
            
            for i = 1:length(sacc1_t_start)
                sacc1_31(trialno_31,i) = histcounts(spktimes_sacc1,[sacc1_t_start(i) sacc1_t_start(i)+binwidth]) - pre_dots_spk;
            end
            for i = 1:length(pursuit_t_start)
                pursuit_31(trialno_31,i) = histcounts(spktimes_pursuit,[pursuit_t_start(i) pursuit_t_start(i)+binwidth]) - pre_dots_spk;
            end
            for i = 1:length(FP2_t_start)
                FP2_31(trialno_31,i) = histcounts(spktimes_FP2,[FP2_t_start(i) FP2_t_start(i)+binwidth]) - pre_dots_spk;
            end
            
            for i = 1:length(choice_t_start)
                choice_31(trialno_31,i) = histcounts(spktimes_choice,[choice_t_start(i) choice_t_start(i)+binwidth]) - pre_dots_spk;
            end
            
            
        end
    end
    
    % calculating kendall-tau and its p-value (config A)
    if trialno_30 > 50
        
        [r_30_P1, p_30_P1] = corr([coh_30(:,1) P1_30],'type','kendall'); 
        [r_30_P2_coh2, p_30_P2_coh2] = corr([coh_30(:,2) P2_30],'type','kendall'); [r_30_P2_cohsum, p_30_P2_cohsum] = corr([sum_coh_30 P2_30],'type','kendall');
        [r_30_sacc1, p_30_sacc1] = corr([coh_30(:,1) sacc1_30],'type','kendall'); [r_30_pursuit, p_30_pursuit] = corr([coh_30(:,1) pursuit_30],'type','kendall');
        [r_30_FP2, p_30_FP2] = corr([coh_30(:,1) FP2_30],'type','kendall'); [r_30_choice, p_30_choice] = corr([sum_coh_30 choice_30],'type','kendall');
        
        CohDep_30{ii+1,2} = r_30_P1(1,2:end);   %P1
        p_30{ii+1,2} = p_30_P1(1,2:end);
        
        CohDep_30{ii+1,3} = r_30_P2_coh2(1,2:end);   % P2-coh2 
        p_30{ii+1,3} = p_30_P2_coh2(1,2:end);
        
        CohDep_30{ii+1,4} = r_30_P2_cohsum(1,2:end);   % P2-sum_coh
        p_30{ii+1,4} = p_30_P2_cohsum(1,2:end);
        
        CohDep_30{ii+1,5} = r_30_sacc1(1,2:end);   % first saccade 
        p_30{ii+1,5} = p_30_sacc1(1,2:end);
        
        CohDep_30{ii+1,6} = r_30_pursuit(1,2:end);  % pursuit
        p_30{ii+1,6} = p_30_pursuit(1,2:end);
        
        CohDep_30{ii+1,7} = r_30_FP2(1,2:end);      % resumed fixation
        p_30{ii+1,7} = p_30_FP2(1,2:end);
        
        CohDep_30{ii+1,8} = r_30_choice(1,2:end);   % choice saccade
        p_30{ii+1,8} = p_30_choice(1,2:end);
        
        
    end
    
    % calculating kendall-tau and its p-value (config B)
    if trialno_31 > 50
        
        [r_31_P1, p_31_P1] = corr([coh_31(:,1) P1_31],'type','kendall'); 
        [r_31_P2_coh2, p_31_P2_coh2] = corr([coh_31(:,2) P2_31],'type','kendall'); [r_31_P2_cohsum, p_31_P2_cohsum] = corr([sum_coh_31 P2_31],'type','kendall');
        [r_31_sacc1, p_31_sacc1] = corr([coh_31(:,1) sacc1_31],'type','kendall'); [r_31_pursuit, p_31_pursuit] = corr([coh_31(:,1) pursuit_31],'type','kendall');
        [r_31_FP2, p_31_FP2] = corr([coh_31(:,1) FP2_31],'type','kendall'); [r_31_choice, p_31_choice] = corr([sum_coh_31 choice_31],'type','kendall');
        
        CohDep_31{ii+1,2} = r_31_P1(1,2:end);   %P1
        p_31{ii+1,2} = p_31_P1(1,2:end);
        
        CohDep_31{ii+1,3} = r_31_P2_coh2(1,2:end);   % P2-coh2 
        p_31{ii+1,3} = p_31_P2_coh2(1,2:end);
        
        CohDep_31{ii+1,4} = r_31_P2_cohsum(1,2:end);   % P2-sum_coh
        p_31{ii+1,4} = p_31_P2_cohsum(1,2:end);
        
        CohDep_31{ii+1,5} = r_31_sacc1(1,2:end);    % first saccade 
        p_31{ii+1,5} = p_31_sacc1(1,2:end);
        
        CohDep_31{ii+1,6} = r_31_pursuit(1,2:end);  % pursuit
        p_31{ii+1,6} = p_31_pursuit(1,2:end);
        
        CohDep_31{ii+1,7} = r_31_FP2(1,2:end);      % resumed fixation
        p_31{ii+1,7} = p_31_FP2(1,2:end);
        
        CohDep_31{ii+1,8} = r_31_choice(1,2:end);   % choice saccade
        p_31{ii+1,8} = p_31_choice(1,2:end);
        
        
    end
    
end

%end

    
    
    
    
    
    
