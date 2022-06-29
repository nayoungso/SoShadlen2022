function Pulse_2ndVar_PSTH_P2coh

% Fig 3G

path = ['~/Dropbox/plxdata/Pulse/Vprobe/prelim_groups/coh_dep_dynamics_based_consec3_p_05_base_adjusted/allo/follower/']

filename_a = dir(strcat(path,'/34/r/','*_DSP*.mat'));
filename_b = dir(strcat(path,'/34/g/','*_DSP*.mat'));
filename_c = dir(strcat(path,'/35/r/','*_DSP*.mat'));
filename_d = dir(strcat(path,'/35/g/','*_DSP*.mat'));


binwidth = 100;
step = 10; 
trialno = 0;

trial_spike_ts = zeros(1,1);

pre_dots2 = -150; post_dots2 = 400;  
base_start = -150; base_end = 150; % for adjusting the baseline FR 


coh_set = [-640 -320 -160 -80 -40 -1 1 40 80 160 320 640];

G1no_1st = zeros(1,12); % coh based on the first pulse 
G1no_2nd = zeros(1,12); % coh based on the 2nd pulse 

for p=1:4
    
    clear filename
    
    if p==1
        taskid = 34; filename = filename_a; filepath = [path,'/34/r/'];
    elseif p==2
        taskid = 34; filename = filename_b; filepath = [path,'/34/g/'];
    elseif p==3
        taskid = 35; filename = filename_c; filepath = [path,'/35/r/'];
    elseif p==4
        taskid = 35; filename = filename_d; filepath = [path,'/35/g/'];
    end
    
    for i = 1:length(filename)
        
        clear trials
        filename(i).name
        load([filepath,filename(i).name]);
        
        for j = 1:length(trials)
            
            if (trials(j).taskid == taskid)&&(trials(j).response>=0)&&(trials(j).mean_maxFR >= 10)
                
                if mod(p,2)==0
                    trials(j).t1_dir = trials(j).t1_dir+180;  
                end
                
                trialno = trialno+1;
                dspsize = find(DSP011(j,:) > trials(j).time_sacc(end)+300,1);
                
                trial_spike_ts(trialno,1:dspsize) = DSP011(j,1:dspsize);
                
                dots_time2(trialno) = trials(j).time_dots_on(2);
                
                t_array_d2(trialno,:) = [pre_dots2 post_dots2];
                
                
                dir_sign = (trials(j).dot_dir ~= trials(j).t1_dir); % 0 for t1_dir / 1 for t2_dir
                coh = trials(j).dot_coh(:)*(1-2*dir_sign);  % coh depending on the 1st pulse
                
                if ((dir_sign == 0)&&(trials(j).response == 1))||((dir_sign == 1)&&(trials(j).response == 0))
                    choice(trialno,1) = 1;
                else
                    choice(trialno,1) = 0;
                end
                
                for k = 1:2
                    if coh(k) == 0
                        if (dir_sign == 0)
                            coh(k) = 1;
                        else
                            coh(k) = -1;
                        end
                    end
                    
                    indext(k) = find(coh_set == coh(k));
                    
                end
                
                G1no_1st(indext(1)) = G1no_1st(indext(1))+1;
                G1address_1st{indext(1)}(G1no_1st(indext(1))) = trialno;
                
                G1no_2nd(indext(2)) = G1no_2nd(indext(2))+1;
                G1address_2nd{indext(2)}(G1no_2nd(indext(2))) = trialno;
                
            end
            
            
        end % end of the trials loop - time to normalize the activity for each cell
        
        
    end
end

clear trials


% plotting

figure(43)
clf;

timex_d2 = pre_dots2:step:post_dots2;

timex_dstart = find(timex_d2==base_start);
timex_dend = find(timex_d2==base_end);

for k = 1:12  
    
    if sum(G1no_1st) > 10
        DotsOn2_2nd(k,:) = spk_psth_v2(trial_spike_ts,G1address_2nd{k},[dots_time2(:)],[pre_dots2 post_dots2],t_array_d2,binwidth,step);
        DotsOn2_1st(k,:) = spk_psth_v2(trial_spike_ts,G1address_1st{k},[dots_time2(:)],[pre_dots2 post_dots2],t_array_d2,binwidth,step);
        
        DotsOn2_2nd_base(k) = nanmean(DotsOn2_2nd(k,timex_dstart:timex_dend));
        
    end
    
end


%for detrending: calculating the average activity across all coherences for subtractaction
DotsOn2(1,:) = spk_psth_v2(trial_spike_ts,[1:trialno],[dots_time2(:)],[pre_dots2 post_dots2],t_array_d2,binwidth,step);
DotsOn2_base = nanmean(DotsOn2(1,timex_dstart:timex_dend));    

maxy = max([max(max(DotsOn2_1st)) max(max(DotsOn2_2nd))])+3;


ColorMapCoh = [0 0.5 0; 0 0.7 0.3; 0.4 0.8 0.6; 0.4 0.8 0.6; 0.6 0.9 0.8; 0.8 1 0.8; ...
                1 1 0.8 ;1 1 0.6; 1 0.9 0.5; 1 0.8 0.4; 1 0.6 0; 1 0 0];

ColorMapCoh3 = [0 0 0.4; 0 0 0.8; 0.2 0.4 0.8; 0.2 0.6 1; 0.2 1 1; 0.8 1 1; ...
                1 0.8 1; 1 0.5 1 ; 1 0.4 0.6; 1 0 0.8; 0.8 0 0.6; 0.6 0 0.4];


% plotting 
ax1 = subplot(131);   
hold on
for i = 1:12
    plot(timex_d2, DotsOn2_1st(i,:), 'Color', ColorMapCoh(i,:), 'LineWidth', 2);
end
line([0 0],[-1 maxy],'Color','k');
xlim([-100 post_dots2]);
set(gca, 'TickDir','out');
ylabel('Firing rates, sorted by P1 coh')
title('aligned on P2')

ax2 = subplot(132);
hold on
for i = 1:12
    plot(timex_d2, DotsOn2_2nd(i,:), 'Color', ColorMapCoh3(i,:), 'LineWidth', 2);
end
line([0 0],[-1 maxy],'Color','k');
xlim([-100 post_dots2]);
set(gca, 'TickDir','out');
ylabel('Firing rates, sorted by P2 coh')
title('aligned on P2')

linkaxes([ax1,ax2],'y');


ax3 = subplot(133);
hold on
for i = 1:12
    plot(timex_d2, DotsOn2_2nd(i,:)-DotsOn2(1,:)-(DotsOn2_2nd_base(i)-DotsOn2_base), 'Color', ColorMapCoh3(i,:), 'LineWidth', 2);
 end
%hold off
%line([0 0],[-1 maxy],'Color','k');
grid on
xlim([-100 post_dots2]);
ylim([-3 3]);
set(gca, 'TickDir','out');
ylabel('Detrended firing rates, sorted by P2 coh')
title('P2 onset')


end