function Pulse_1stVar_PSTH(celltype)

if celltype==1
    groupname = 'Leader';
elseif celltype==2
    groupname = 'Supporter';
end
    
path = ['~/Two-pulse task, 1st variant/']

filename_a = dir(strcat(path,'/A/r/','*_DSP*.mat'));    % configuration A, up/right-preferring neurons
filename_b = dir(strcat(path,'/A/g/','*_DSP*.mat'));    % configuration A, down/left-preferring neurons
filename_c = dir(strcat(path,'/B/r/','*_DSP*.mat'));    % configuration B, up/right-preferring neurons
filename_d = dir(strcat(path,'/B/g/','*_DSP*.mat'));    % configuration B, down/left-preferring neurons


binwidth = 100;
step = 20;

pre_dots_det = -300; post_dots_det = 0;

pre_dots = -50;  post_dots = 400;
pre_sacc1 = -300;  post_sacc1 = 350;
pre_slib = -350; post_slib = 350;
pre_fp2 = -400; post_fp2 = 400;
pre_dots2 = -200;  post_dots2 = 400;
pre_sacc = -200; post_sacc = 200;

trialno = 0; trialno_choice = 0;
trial_spike_ts = zeros(1,1);
trial_spike_choice_ts = zeros(1,1);

coh_set = [-640 -320 -160 -80 -40 -1 1 40 80 160 320 640];  
G1no = zeros(1,12); G1no_choice = zeros(1,12);

for p = 1:4
    
    clear filename
    
    if p==1
        taskid = 30; filename = filename_a; filepath = [path,'/A/r/']; % neurons preferring right/up target in config A
    elseif p==2
        taskid = 30; filename = filename_b; filepath = [path,'/A/g/']; % neurons preferring left/down target in config A
    elseif p==3
        taskid = 31; filename = filename_c; filepath = [path,'/B/r/']; % neurons preferring right/up target in config B
    elseif p==4
        taskid = 31; filename = filename_d; filepath = [path,'/B/g/']; % neurons preferring left/down target in config B
    end
    
    for i = 1:length(filename)
        
        clear trials
        filename(i).name
        load([filepath,filename(i).name]);
        
        for j = 1:length(trials)
            
            if (trials(j).taskid == taskid)&&(trials(j).response >=0)&&(trials(j).mean_maxFR >= 10)
                
                if mod(p,2)==0
                    trials(j).t1_dir = trials(j).t1_dir+180;  % flip the sign for left/down preferring neurons (p=2 or 4) to be combined with right/up preferring neurons
                end
                
                
                trialno = trialno+1;
                dspsize = find(DSP011(j,:) > trials(j).time_sacc(end)+300,1);
                trial_spike_ts(trialno,1:dspsize) = DSP011(j,1:dspsize);
                dir_sign = (trials(j).dot_dir ~= trials(j).t1_dir); % 0 for t1_dir / 1 for t2_dir
                coh = trials(j).dot_coh(:)*(1-2*dir_sign);
                
                for k = 1:length(coh)
                    if coh(k) == 0
                        if (dir_sign == 0)
                            coh(k) = 1;
                        else
                            coh(k) = -1;
                        end
                    end
                    
                    indext(k) = find(coh_set == coh(k));
                    
                end
                
                dots_time1(trialno) = trials(j).time_dots_on(1);
                sacc1_times(trialno) = trials(j).time_sacc(1);
                slib_times(trialno) = trials(j).time_slib_start(1);
                fix2_times(trialno) = trials(j).time_fp_acq(2);
                dots_time2(trialno) = trials(j).time_dots_on(2);
                
                t_array_d(trialno,:) = [pre_dots post_dots];
                t_array_s1(trialno,:) = [pre_sacc1, post_sacc1];
                t_array_sl(trialno,:) = [pre_slib post_slib];
                t_array_fp2(trialno,:) = [pre_fp2 post_fp2];
                t_array_d2(trialno,:) = [pre_dots2 post_dots2];
                
                %for baseline adjustment to mitigate the drift effect (subtracting the pre-motion activity)
                t_array_det(trialno,:) = [pre_dots_det post_dots_det];
                
                G1no(indext(1)) = G1no(indext(1))+1;
                G1address{indext(1)}(G1no(indext(1))) = trialno;
                
                
                if trials(j).response ==1
                    trialno_choice=trialno_choice+1;
                    trial_spike_choice_ts(trialno_choice,1:dspsize) = DSP011(j,1:dspsize);
                
                    choice_sacc_times(trialno_choice) = trials(j).time_sacc(end);
                    t_array_s(trialno_choice,:) = [pre_sacc post_sacc];
                    t_array_choice_det(trialno_choice,:) = [pre_dots_det post_dots_det];
                
                    G1no_choice(indext(1)) = G1no_choice(indext(1))+1;
                    G1address_choice{indext(1)}(G1no_choice(indext(1))) = trialno_choice;
                end
                
                
                
                
            end
            
        end
        
    end
end


for i = 1:length(coh_set) % to calculate residual FR to subtract
    
    [spk_cnt_temp, fr_temp] = spk_cnt2(trial_spike_ts,G1address{i},[dots_time1(:)],[pre_dots_det post_dots_det],t_array_det);
    FR_grand{i} = fr_temp*1000; % currently, fr_temp = spk cnt/duration (duration in ms unit)
    
    [spk_cnt_temp, fr_temp] = spk_cnt2(trial_spike_choice_ts,G1address_choice{i},[dots_time1(:)],[pre_dots_det post_dots_det],t_array_choice_det);
    FR_grand_choice{i} = fr_temp*1000; % currently, fr_temp = spk cnt/duration (duration in ms unit)
end


for i = 1:length(coh_set)
    
    % baseline adjustment (gathering the residual FR for each condition)
    Det(i) = nanmean(FR_grand{i})-nanmean(cell2mat(FR_grand));
    Det_choice(i) = nanmean(FR_grand_choice{i})-nanmean(cell2mat(FR_grand_choice));
    
    DotsOn{i} = spk_psth_v2(trial_spike_ts,G1address{i},[dots_time1(:)],[pre_dots post_dots],t_array_d,binwidth,step) - Det(i);
    
    Sacc1{i} = spk_psth_v2(trial_spike_ts,G1address{i},[sacc1_times(:)],[pre_sacc1 post_sacc1],t_array_s1,binwidth,step) - Det(i);
    Slib{i} = spk_psth_v2(trial_spike_ts,G1address{i},[slib_times(:)],[pre_slib post_slib],t_array_sl,binwidth,step) - Det(i);
    LastFix{i} = spk_psth_v2(trial_spike_ts,G1address{i},[fix2_times(:)],[pre_fp2 post_fp2],t_array_fp2,binwidth,step) - Det(i);
    DotsOn2{i} = spk_psth_v2(trial_spike_ts,G1address{i},[dots_time2(:)],[pre_dots2 post_dots2],t_array_d2,binwidth,step) - Det(i);
    ChoiceSacc{i} = spk_psth_v2(trial_spike_choice_ts,G1address_choice{i},[choice_sacc_times(:)],[pre_sacc post_sacc],t_array_s,binwidth,step) - Det_choice(i);
    
end


timex_s = pre_sacc:step:post_sacc;
timex_s1 = pre_sacc1:step:post_sacc1;
timex_sl = pre_slib:step:post_slib;
timex_d = pre_dots:step:post_dots;
timex_d2 = pre_dots2:step:post_dots2;
timex_fp2 = pre_fp2:step:post_fp2;


maxy1 = max([Sacc1{12} DotsOn{12} ChoiceSacc{12} Slib{12}]);
maxy1 = maxy1+5;

% separating 0% coh trials

ColorMapCoh = [0 0.5 0; 0 0.7 0.3; 0.4 0.8 0.6; 0.4 0.8 0.6; 0.6 0.9 0.8; 0.8 1 0.8; ...
    1 1 0.8 ;1 1 0.6; 1 0.9 0.5; 1 0.8 0.4; 1 0.6 0; 1 0 0];

figure(11)
clf;

%% dots
ax1 = subplot(161);
title('P1 on')
hold on
for i = 1:length(coh_set)
    plot(timex_d, DotsOn{i}, 'Color', ColorMapCoh(i,:),'LineWidth',2);
end
line([0 0],[-1 maxy1],'Color','k');
hold off
xlim([pre_dots post_dots]);
ylim([-1 maxy1]);
set(gca, 'TickDir','out');


%% first sacc
ax2=subplot(162);
title('sacc to T^0')
hold on
for i = 1:length(coh_set)
    plot(timex_s1, Sacc1{i}, 'Color', ColorMapCoh(i,:),'LineWidth',2);
end
line([0 0],[-1 maxy1],'Color','k');
hold off
xlim([pre_sacc1 post_sacc1]);
ylim([-1 maxy1]);
set(gca, 'TickDir','out');

%% pursuit
ax3 = subplot(163);
title('pursuit starts')
hold on
for i = 1:length(coh_set)
    plot(timex_sl, Slib{i}, 'Color', ColorMapCoh(i,:),'LineWidth',2);
end
line([0 0],[-1 maxy1],'Color','k');
hold off
xlim([pre_slib post_slib]);
ylim([-1 maxy1]);
set(gca, 'TickDir','out');


%% regained fix
ax4 = subplot(164);
title('regained fix');
hold on
for i = 1:length(coh_set)
    plot(timex_fp2, LastFix{i}, 'Color', ColorMapCoh(i,:),'LineWidth',2);
end
line([0 0],[-1 maxy1],'Color','k');
hold off
xlim([pre_fp2 post_fp2]);
ylim([-1 maxy1]);
set(gca, 'TickDir','out');

%% dots #2
ax5 = subplot(165);
title('P2 on')
hold on
for i = 1:length(coh_set)
    plot(timex_d2, DotsOn2{i}, 'Color', ColorMapCoh(i,:),'LineWidth',2);
end
line([0 0],[-1 maxy1],'Color','k');
hold off
xlim([pre_dots2 post_dots2]);
ylim([-1 maxy1]);
set(gca, 'TickDir','out');

%% choice sacc
ax6 = subplot(166);
title('choice sacc');
hold on
for i = 1:length(coh_set)
    plot(timex_s, ChoiceSacc{i}, 'Color',ColorMapCoh(i,:),'LineWidth',2);
end
line([0 0],[-1 maxy1],'Color','k');
hold off
xlim([pre_sacc post_sacc]);
ylim([-1 maxy1]);
set(gca, 'TickDir','out');

totalcells_used = length(filename_a)+length(filename_b)+length(filename_c)+length(filename_d)

toc
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'y');


suptitle([groupname,'  n = ',num2str(totalcells_used)])

set(figure(11),'Position',[100 300 1200 300])

clear all

end

