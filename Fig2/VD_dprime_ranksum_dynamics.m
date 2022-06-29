%function VD_dprime_ranksum_dynamics

path = ['~/Dropbox/plxdata/RDM_IEM/AllCells/']
filename = dir(strcat(path,'*DSP*.mat'))

binwidth = 300;
stepsize = 50;

Dprime_20{1,1} = 'File Name'; Dprime_20{1,2} = 'dots'; Dprime_20{1,3} = 'sacc1'; 
Dprime_20{1,4} = 'pursuit'; Dprime_20{1,5} = 'FP2';  Dprime_20{1,6} = 'choice'; 

Dprime_21{1,1} = 'File Name'; Dprime_21{1,2} = 'dots'; Dprime_21{1,3} = 'sacc1'; 
Dprime_21{1,4} = 'pursuit'; Dprime_21{1,5} = 'FP2'; Dprime_21{1,6} = 'choice'; 

RanksumP_20{1,1} = 'File Name'; RanksumP_20{1,2} = 'dots'; RanksumP_20{1,3} = 'sacc1'; 
RanksumP_20{1,4} = 'pursuit'; RanksumP_20{1,5} = 'FP2';  RanksumP_20{1,6} = 'choice'; 

RanksumP_21{1,1} = 'File Name'; RanksumP_21{1,2} = 'dots'; RanksumP_21{1,3} = 'sacc1'; 
RanksumP_21{1,4} = 'pursuit'; RanksumP_21{1,5} = 'FP2'; RanksumP_21{1,6} = 'choice'; 

base_time = [-300 0]; % activity before dots onset will be subtracted (a device for correcting drift in the multi-channel recording)

% time windows of the survey -- note this 'survey' windows are wider than the 'criteria' windows for cell categorization 
dots_time = [-50 600]; sacc1_time = [-400 500]; pursuit_time = [-400 600];
FP2_time = [-400 300]; choice_time = [-500 200];

dots_t_start = dots_time(1)-binwidth/2:stepsize:dots_time(end)-binwidth/2;  
sacc1_t_start = sacc1_time(1)-binwidth/2:stepsize:sacc1_time(end)-binwidth/2;
pursuit_t_start = pursuit_time(1)-binwidth/2:stepsize:pursuit_time(end)-binwidth/2;
FP2_t_start = FP2_time(1)-binwidth/2:stepsize:FP2_time(end)-binwidth/2;
choice_t_start = choice_time(1)-binwidth/2:stepsize:choice_time(end)-binwidth/2;

% storing the time windows of the survey 
TimeArray{1,1} = 'Dots time'; TimeArray{1,2} = dots_t_start; 
TimeArray{2,1} = 'Sacc1 time'; TimeArray{2,2} = sacc1_t_start;
TimeArray{3,1} = 'Pursuit time'; TimeArray{3,2} = pursuit_t_start;
TimeArray{4,1} = 'FP2 time'; TimeArray{4,2} = FP2_t_start;
TimeArray{5,1} = 'Choice time'; TimeArray{5,2} = choice_t_start;
TimeArray{6,1} = 'binwidth'; TimeArray{6,2} = binwidth;
TimeArray{7,1} = 'stepsize'; TimeArray{7,2} = stepsize;

for ii = 1:length(filename)
    clear trials
    filename(ii).name
    
    clear coh_20 cho_20 dots_20 sacc1_20 pursuit_20 FP2_20 choice_20
    clear coh_21 cho_21 dots_21 sacc1_21 pursuit_21 FP2_21 choice_21
    clear dots_20_det sacc1_20_det pursuit_20_det FP2_20_det choice_20_det
    clear dots_21_det sacc1_21_det pursuit_21_det FP2_21_det choice_21_det
    
    
    fullfilename = filename(ii).name;
    load([path,filename(ii).name]);
    
    Dprime_20{ii+1,1} = fullfilename;   % giving up storing the simpler cell name here
    Dprime_21{ii+1,1} = fullfilename;   
    
    RanksumP_20{ii+1,1} = fullfilename;
    RanksumP_21{ii+1,1} = fullfilename;
    
    dspsize = size(DSP011,2);
    
    trialno_20 = 0;
    trialno_21 = 0;
    
    for j = 1:length(trials)
        
        clear spktimes spktimes_base spktimes_dots spktimes_sacc1 spktimes_pursuit spktimes_FP2 spktimes_choice
        try 
        % configuration A    
        if (trials(j).taskid==20)&&(trials(j).response>=0)&&(trials(j).mean_maxFR >= 10)    
            
            trialno_20 = trialno_20+1;
            spktimes = DSP011(j,1:dspsize);
            
            spktimes_base = spktimes - trials(j).time_dots_on(1);
            
            spktimes_dots = spktimes - trials(j).time_dots_on(1);
            spktimes_sacc1 = spktimes - trials(j).time_sacc(1);
            spktimes_pursuit = spktimes - trials(j).time_slib_start;
            spktimes_FP2 = spktimes - trials(j).time_fp_acq(2);
            spktimes_choice = spktimes - trials(j).time_sacc(end);
            
            pre_dots_spk = histcounts(spktimes_base,base_time); % for baseline adjustment
            
            % assign a sign for the motion coherence 
            dir_sign = (trials(j).dot_dir ~= trials(j).t1_dir); 
            coh_20(trialno_20,:) = trials(j).dot_coh*(1-2*dir_sign);
            
            % choice information (1 for right/up motion choice, 0 for left/down motion choice)
            if ((dir_sign == 0)&&(trials(j).response == 1))||((dir_sign == 1)&&(trials(j).response == 0))
                cho_20(trialno_20,:) = 1;
            else
                cho_20(trialno_20,:) = 0;
            end
            
            % activity around motion onset
            for i = 1:length(dots_t_start)
                dots_20(trialno_20,i) = histcounts(spktimes_dots,[dots_t_start(i) dots_t_start(i)+binwidth]) - pre_dots_spk;
            end
            
            % activity around first saccade
            for i = 1:length(sacc1_t_start)
                sacc1_20(trialno_20,i) = histcounts(spktimes_sacc1,[sacc1_t_start(i) sacc1_t_start(i)+binwidth]) - pre_dots_spk;
            end
            
            % activity around pursuit onset
            for i = 1:length(pursuit_t_start)
                pursuit_20(trialno_20,i) = histcounts(spktimes_pursuit,[pursuit_t_start(i) pursuit_t_start(i)+binwidth]) - pre_dots_spk;
            end
            
            % activity around resumed fixation
            for i = 1:length(FP2_t_start)
                FP2_20(trialno_20,i) = histcounts(spktimes_FP2,[FP2_t_start(i) FP2_t_start(i)+binwidth]) - pre_dots_spk;
            end
            
            % activity around choice
            for i = 1:length(choice_t_start)
                choice_20(trialno_20,i) = histcounts(spktimes_choice,[choice_t_start(i) choice_t_start(i)+binwidth]) - pre_dots_spk;
            end
            
            
            % repeat for config B
        elseif (trials(j).taskid==21)&&(trials(j).response>=0)&&(trials(j).mean_maxFR >= 10)   
            
            trialno_21 = trialno_21+1;
            spktimes = DSP011(j,1:dspsize);
            
            spktimes_base = spktimes - trials(j).time_dots_on(1);
            
            spktimes_dots = spktimes - trials(j).time_dots_on(1);
            spktimes_sacc1 = spktimes - trials(j).time_sacc(1);
            spktimes_pursuit = spktimes - trials(j).time_slib_start;
            spktimes_FP2 = spktimes - trials(j).time_fp_acq(2);
            spktimes_choice = spktimes - trials(j).time_sacc(end);
            
            pre_dots_spk = histcounts(spktimes_base,base_time);
            
            dir_sign = (trials(j).dot_dir ~= trials(j).t1_dir); 
            coh_21(trialno_21,:) = trials(j).dot_coh*(1-2*dir_sign);
            
            if ((dir_sign == 0)&&(trials(j).response == 1))||((dir_sign == 1)&&(trials(j).response == 0))
                cho_21(trialno_21,:) = 1;
            else
                cho_21(trialno_21,:) = 0;
            end
            
            for i = 1:length(dots_t_start)
                dots_21(trialno_21,i) = histcounts(spktimes_dots,[dots_t_start(i) dots_t_start(i)+binwidth]) - pre_dots_spk;
            end
            
            for i = 1:length(sacc1_t_start)
                sacc1_21(trialno_21,i) = histcounts(spktimes_sacc1,[sacc1_t_start(i) sacc1_t_start(i)+binwidth]) - pre_dots_spk;
            end
            for i = 1:length(pursuit_t_start)
                pursuit_21(trialno_21,i) = histcounts(spktimes_pursuit,[pursuit_t_start(i) pursuit_t_start(i)+binwidth]) - pre_dots_spk;
            end
            for i = 1:length(FP2_t_start)
                FP2_21(trialno_21,i) = histcounts(spktimes_FP2,[FP2_t_start(i) FP2_t_start(i)+binwidth]) - pre_dots_spk;
            end
            
            for i = 1:length(choice_t_start)
                choice_21(trialno_21,i) = histcounts(spktimes_choice,[choice_t_start(i) choice_t_start(i)+binwidth]) - pre_dots_spk;
            end
            
            
        end
        catch
        end
    end
    
    
    % calculating the ranksum p values & d' for each time point 
    if trialno_20 > 50
        
        clear cho1ind cho2ind
        cho1ind = cho_20==1; cho2ind = cho_20==0;
        
        for i = 1:length(dots_t_start)
            
            RanksumP_20{ii+1,2}(i) = ranksum(dots_20(cho1ind,i),dots_20(cho2ind,i)); % dots
            
            diff_mu = nanmean(dots_20(cho1ind,i))-nanmean(dots_20(cho2ind,i));
            sum_var = nanvar(dots_20(cho1ind,i))+nanvar(dots_20(cho2ind,i));
            Dprime_20{ii+1,2}(i) = diff_mu/sqrt(sum_var/2);
            
        end
        
        for i = 1:length(sacc1_t_start)
            
            RanksumP_20{ii+1,3}(i) = ranksum(sacc1_20(cho1ind,i),sacc1_20(cho2ind,i)); % dots
            
            diff_mu = nanmean(sacc1_20(cho1ind,i))-nanmean(sacc1_20(cho2ind,i));
            sum_var = nanvar(sacc1_20(cho1ind,i))+nanvar(sacc1_20(cho2ind,i));
            Dprime_20{ii+1,3}(i) = diff_mu/sqrt(sum_var/2);
            
        end
            
        for i = 1:length(pursuit_t_start)
            
            RanksumP_20{ii+1,4}(i) = ranksum(pursuit_20(cho1ind,i),pursuit_20(cho2ind,i)); % dots
            
            diff_mu = nanmean(pursuit_20(cho1ind,i))-nanmean(pursuit_20(cho2ind,i));
            sum_var = nanvar(pursuit_20(cho1ind,i))+nanvar(pursuit_20(cho2ind,i));
            Dprime_20{ii+1,4}(i) = diff_mu/sqrt(sum_var/2);
            
            
        end
        
        for i = 1:length(FP2_t_start)
            
            RanksumP_20{ii+1,5}(i) = ranksum(FP2_20(cho1ind,i),FP2_20(cho2ind,i)); % dots
            
            diff_mu = nanmean(FP2_20(cho1ind,i))-nanmean(FP2_20(cho2ind,i));
            sum_var = nanvar(FP2_20(cho1ind,i))+nanvar(FP2_20(cho2ind,i));
            Dprime_20{ii+1,5}(i) = diff_mu/sqrt(sum_var/2);
            
            
        end
        
        for i = 1:length(choice_t_start)
            
            RanksumP_20{ii+1,6}(i) = ranksum(choice_20(cho1ind,i),choice_20(cho2ind,i)); % dots
            
            diff_mu = nanmean(choice_20(cho1ind,i))-nanmean(choice_20(cho2ind,i));
            sum_var = nanvar(choice_20(cho1ind,i))+nanvar(choice_20(cho2ind,i));
            Dprime_20{ii+1,6}(i) = diff_mu/sqrt(sum_var/2);
            
            
        end
        
    end
    
    % config B
    if trialno_21 > 50
        
        clear cho1ind cho2ind
        cho1ind = cho_21==1; cho2ind = cho_21==0;
        
        for i = 1:length(dots_t_start)
            
            RanksumP_21{ii+1,2}(i) = ranksum(dots_21(cho1ind,i),dots_21(cho2ind,i)); % dots
            
            diff_mu = nanmean(dots_21(cho1ind,i))-nanmean(dots_21(cho2ind,i));
            sum_var = nanvar(dots_21(cho1ind,i))+nanvar(dots_21(cho2ind,i));
            Dprime_21{ii+1,2}(i) = diff_mu/sqrt(sum_var/2);
            
        end
        
        for i = 1:length(sacc1_t_start)
            
            RanksumP_21{ii+1,3}(i) = ranksum(sacc1_21(cho1ind,i),sacc1_21(cho2ind,i)); % dots
            
            diff_mu = nanmean(sacc1_21(cho1ind,i))-nanmean(sacc1_21(cho2ind,i));
            sum_var = nanvar(sacc1_21(cho1ind,i))+nanvar(sacc1_21(cho2ind,i));
            Dprime_21{ii+1,3}(i) = diff_mu/sqrt(sum_var/2);
            
            
        end
            
        for i = 1:length(pursuit_t_start)
            
            RanksumP_21{ii+1,4}(i) = ranksum(pursuit_21(cho1ind,i),pursuit_21(cho2ind,i)); % dots
            
            diff_mu = nanmean(pursuit_21(cho1ind,i))-nanmean(pursuit_21(cho2ind,i));
            sum_var = nanvar(pursuit_21(cho1ind,i))+nanvar(pursuit_21(cho2ind,i));
            Dprime_21{ii+1,4}(i) = diff_mu/sqrt(sum_var/2);
            
        end
        
        for i = 1:length(FP2_t_start)
            
            RanksumP_21{ii+1,5}(i) = ranksum(FP2_21(cho1ind,i),FP2_21(cho2ind,i)); % dots
            
            diff_mu = nanmean(FP2_21(cho1ind,i))-nanmean(FP2_21(cho2ind,i));
            sum_var = nanvar(FP2_21(cho1ind,i))+nanvar(FP2_21(cho2ind,i));
            Dprime_21{ii+1,5}(i) = diff_mu/sqrt(sum_var/2);
            
        end
        
        for i = 1:length(choice_t_start)
            
            RanksumP_21{ii+1,6}(i) = ranksum(choice_21(cho1ind,i),choice_21(cho2ind,i)); % dots
            
            diff_mu = nanmean(choice_21(cho1ind,i))-nanmean(choice_21(cho2ind,i));
            sum_var = nanvar(choice_21(cho1ind,i))+nanvar(choice_21(cho2ind,i));
            Dprime_21{ii+1,6}(i) = diff_mu/sqrt(sum_var/2);
            
        end
        
    end
    
end


    
    
    
    
    
    
