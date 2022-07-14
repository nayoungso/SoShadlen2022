%function Pulse_2ndVar_P2_regression_eq_8_9


% correct trials only and remove choice term
% R1 = activity at P2 onset
% R2 = activity after P2 (200,250,300,350,...,550)
% C1 = coherence of P1, C2 = coherence of P2

% two separate regressions: 
% R1 = a0+a1C1          (eq.8)
% R2 = (a0+a1C1)+b1C2   (eq.9)


path = ['~/Two-pulse task, 2nd variant/Supporter']

filename_a = dir(strcat(path,'/A/r/','*_DSP*.mat'));    % configuration A, up/right-preferring neurons
filename_b = dir(strcat(path,'/A/g/','*_DSP*.mat'));    % configuration A, down/left-preferring neurons
filename_c = dir(strcat(path,'/B/r/','*_DSP*.mat'));    % configuration B, up/right-preferring neurons
filename_d = dir(strcat(path,'/B/g/','*_DSP*.mat'));    % configuration B, down/left-preferring neurons


binwidth = 100;
step = 10; 
crit_fr = 6;
trialno = 0;
cellno = 0;

trial_spike_ts = zeros(1,1);


glm_binwidth = 100;
glm_init_t = 0;
glm_time_array = 200:50:550; % different values of the end time for P2-viewing epoch (when R2 is measured). Results hold for all values tested.


for p=1:4
    
    clear filename
    
    if p==1
        taskid = 34; filename = filename_a; filepath = [path,'/A/r/'];
    elseif p==2
        taskid = 34; filename = filename_b; filepath = [path,'/A/g/'];
    elseif p==3
        taskid = 35; filename = filename_c; filepath = [path,'/B/r/'];
    elseif p==4
        taskid = 35; filename = filename_d; filepath = [path,'/B/g/'];
    end
    
    for i = 1:length(filename)
        
        clear trials
        filename(i).name
        load([filepath,filename(i).name]);
        
        cellno = cellno+1;
        cell_trialno = 0;
        
        for j = 1:length(trials)
            
            if (trials(j).taskid == taskid)&&(trials(j).response==1)&&(trials(j).mean_maxFR >= crit_fr) % correct trials only
                
                if mod(p,2)==0
                    trials(j).t1_dir = trials(j).t1_dir+180;  
                end
                
                trialno = trialno+1;
                dspsize = find(DSP011(j,:) > trials(j).time_sacc(end)+300,1);
                trial_spike_ts(trialno,1:dspsize) = DSP011(j,1:dspsize);
                
                dots_time2(trialno) = trials(j).time_dots_on(2);
                
                dir_sign = (trials(j).dot_dir ~= trials(j).t1_dir); 
                coh = trials(j).dot_coh(:)*(1-2*dir_sign);  
                
                
                if ((dir_sign == 0)&&(trials(j).response == 1))||((dir_sign == 1)&&(trials(j).response == 0))
                    choice(trialno,1) = 1;
                else
                    choice(trialno,1) = 0;
                end
                
                
                coh1(trialno,1) = coh(1); coh2(trialno,1) = coh(2);
                
                % calculating deltaFR
                cell_trialno = cell_trialno+1;
               
                deltaFR_temp = histcounts(trial_spike_ts(trialno,:)-dots_time2(trialno),[glm_init_t-glm_binwidth/2 glm_init_t+glm_binwidth/2]);
                FR0{cellno}(cell_trialno,1) = deltaFR_temp(1);   % FR before P2
                    

                for k = 1:length(glm_time_array)
                    
                    deltaFR_temp = histcounts(trial_spike_ts(trialno,:)-dots_time2(trialno),[glm_time_array(k)-glm_binwidth/2 glm_time_array(k)+glm_binwidth/2]);
                    FRpre{cellno}(cell_trialno,k) = deltaFR_temp(1);   % FR before P2
                    
                end
                
            end
            
            
        end % end of the trials loop 
        
        % z-scoring FR for each cell
        if cellno == 1
            FR0_z = FR0{cellno}-nanmean(FR0{cellno},1);   % residual FR
            FR0_z = FR0_z./nanstd(FR0{cellno},1);
            
            FRpre_z = FRpre{cellno}-nanmean(FRpre{cellno},1);   % residual FR
            FRpre_z = FRpre_z./nanstd(FRpre{cellno},1);
            
            
        else
            
            clear FRpre_z_temp FR0_z_temp 
            
            FR0_z_temp = FR0{cellno}-nanmean(FR0{cellno},1);   % residual FR
            FR0_z_temp = FR0_z_temp./nanstd(FR0{cellno},1);
            FR0_z = [FR0_z; FR0_z_temp];
            
            FRpre_z_temp = FRpre{cellno}-nanmean(FRpre{cellno},1);   % residual FR
            FRpre_z_temp = FRpre_z_temp./nanstd(FRpre{cellno},1);
            FRpre_z = [FRpre_z; FRpre_z_temp];
            
        end
        
        
    end
end

clear trials

coh1 = coh1/1000; coh2 = coh2/1000;

coh1_p = coh1(coh1+coh2>=0);
FR0_p = FR0_z(coh1+coh2>=0); % correct T+ choice only
coh2_p = coh2(coh1+coh2>=0);

coh1_n = coh1(coh1+coh2<=0);
FR0_n = FR0_z(coh1+coh2<=0); % correct T- choice only
coh2_n = coh2(coh1+coh2<=0);

for k = 1:length(glm_time_array)
    FRpre_p(:,k) = FRpre_z(coh1+coh2>=0,k); 
    FRpre_n(:,k) = FRpre_z(coh1+coh2<=0,k); 
end

% eq. 3
[a,~,stat_a] = glmfit(coh1,FR0_z,'normal','link','identity');
[a_p,~,stat_a_p] = glmfit(coh1_p,FR0_p,'normal','link','identity');
[a_n,~,stat_a_n] = glmfit(coh1_n,FR0_n,'normal','link','identity');

Fitted_base = a(1)+a(2)*coh1;   
Fitted_base_p = a_p(1)+a_p(2)*coh1_p;
Fitted_base_n = a_n(1)+a_n(2)*coh1_n;

% eq. 4
for k = 1:length(glm_time_array)
    clear toBeFit toBeFit_p toBeFit_n
    toBeFit = FRpre_z(:,k)-Fitted_base;
    [b, ~, stat] = glmfit(coh2,toBeFit,'normal','link','identity','Constant','off');
    
    toBeFit_p = FRpre_p(:,k)-Fitted_base_p;
    [b_p, ~, stat_p] = glmfit(coh2_p,toBeFit_p,'normal','link','identity','Constant','off');
    
    toBeFit_n = FRpre_n(:,k)-Fitted_base_n;
    [b_n, ~, stat_n] = glmfit(coh2_n,toBeFit_n,'normal','link','identity','Constant','off');
    
    bb(:,k) = b; bb_p(:,k) = b_p; bb_n(:,k) = b_n;
    pp(:,k) = stat.p; pp_p(:,k) = stat_p.p; pp_n(:,k) = stat_n.p;

end

stat_a.p
stat_a_p.p
stat_a_n.p

