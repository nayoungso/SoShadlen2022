function VD_PSTH(celltype)

if celltype==1
    groupname = 'classic_leader';
elseif celltype==2
    groupname = 'classic_supporter';
end
    
path = ['~/Dropbox/plxdata/RDM_IEM/groups/ranksum_dynamics_based_consec3_p_05_base_adjusted/' groupname '/']

filename_a = dir(strcat(path,'20/r/','*DSP*.mat'))
filename_b = dir(strcat(path,'20/g/','*DSP*.mat'))
filename_c = dir(strcat(path,'21/r/','*DSP*.mat'))
filename_d = dir(strcat(path,'21/g/','*DSP*.mat'))


cutoff = 250;    % 250ms after dots off (for attrition rule application)

binwidth = 100;
step = 20;

pre_dots_base = -300; post_dots_base = 0;  % baseline adjustment to mitigate the effect of drift (in a multi-channel recording, the electrode is not moved)

pre_dots = -50;  post_dots = 450;
pre_sacc1 = -350;  post_sacc1 = 250;
pre_slib = -250; post_slib = 350;
pre_fp2 = -400; post_fp2 = 200;
pre_sacc = -400; post_sacc = 200;

trialno_dots = 0; trialno = 0;
trial_spike_ts = zeros(1,1);
coh_set = [-640 -320 -160 -80 -40 -1 1 40 80 160 320 640];

G1no_dots = zeros(1,12); G1no = zeros(1,12);

for p = 1:4
    
    clear filename
    
    if p==1
        taskid = 20; filename = filename_a; filepath = [path,'/20/r/']; % neurons preferring right/up target in config A
    elseif p==2
        taskid = 20; filename = filename_b; filepath = [path,'/20/g/']; % neurons preferring left/down target in config A
    elseif p==3
        taskid = 21; filename = filename_c; filepath = [path,'/21/r/']; % neurons preferring right/up target in config B
    elseif p==4
        taskid = 21; filename = filename_d; filepath = [path,'/21/g/']; % neurons preferring left/down target in config B
    end
    
    
    for i = 1:length(filename)
        
        clear trials
        load([filepath,filename(i).name]);
        
        for j = 1:length(trials)
            
            if (trials(j).taskid == taskid)&&(trials(j).response >=0)&&(trials(j).mean_maxFR >= 10)
                
                if mod(p,2)==0
                    trials(j).t1_dir = trials(j).t1_dir+180;  % flip the sign for left/down preferring neurons (p=2 or 4) to be combined with right/up preferring neurons
                end
                
                trialno_dots = trialno_dots+1;
                dspsize = find(DSP011(j,:) > trials(j).time_sacc(end)+300,1);
                dir_sign = (trials(j).dot_dir ~= trials(j).t1_dir);
                coh = trials(j).dot_coh*(1-2*dir_sign); % signed coherence
                
                if coh == 0
                    if ((dir_sign == 0)&&(trials(j).response == 1))||((dir_sign==1)&&(trials(j).response == 0))
                        coh = 1;
                    else
                        coh = -1;
                    end
                end
                
                indext = find(coh_set == coh);
                
                
                trial_spike_dots_ts(trialno_dots,1:dspsize) = DSP011(j,1:dspsize);
                dots_time1(trialno_dots) = trials(j).time_dots_on(1);
                tempt = trials(j).time_dots_off(1)-trials(j).time_dots_on(1);
                t_array_d(trialno_dots,:) = [pre_dots min(post_dots,tempt+cutoff)];
                %to adjust the baseline
                t_array_base_d(trialno_dots,:) = [pre_dots_base post_dots_base];
                G1no_dots(indext) = G1no_dots(indext)+1;
                G1address_dots{indext}(G1no_dots(indext)) = trialno_dots;
                
                if trials(j).response == 1  % correct trials only
                    trialno = trialno+1;
                    trial_spike_ts(trialno,1:dspsize) = DSP011(j,1:dspsize);
                    
                    sacc1_times(trialno) = trials(j).time_sacc(1);
                    slib_times(trialno) = trials(j).time_slib_start(1);
                    fix2_times(trialno) = trials(j).time_fp_acq(2);
                    choice_sacc_times(trialno) = trials(j).time_sacc(end);
                    
                    t_array_s1(trialno,:) = [pre_sacc1, post_sacc1];
                    t_array_sl(trialno,:) = [pre_slib post_slib];
                    t_array_fp2(trialno,:) = [pre_fp2 post_fp2];
                    t_array_s(trialno,:) = [pre_sacc post_sacc];
                    
                    %to adjust the baseline
                    t_array_base(trialno,:) = [pre_dots_base post_dots_base];
                    
                    G1no(indext) = G1no(indext)+1;
                    G1address{indext}(G1no(indext)) = trialno;
                end
                
            end
            
        end
        
    end
    
end


for i = 1:length(coh_set) % to calculate residual FR to adjust the baseline
    
    % for motion-viewing epoch, where all trials are shown
    [spk_cnt_temp, fr_temp] = spk_cnt2(trial_spike_dots_ts,G1address_dots{i},[dots_time1(:)],[pre_dots_base post_dots_base],t_array_base_d);
    FR_grand_dots{i} = fr_temp*1000; % currently, fr_temp = spk cnt/duration (duration in ms unit)
    
    % for the rest of the epochs, where only the correct trials are shown
    [spk_cnt_temp, fr_temp] = spk_cnt2(trial_spike_ts,G1address{i},[dots_time1(:)],[pre_dots_base post_dots_base],t_array_base);
    FR_grand{i} = fr_temp*1000; % currently, fr_temp = spk cnt/duration (duration in ms unit)
    
end


for i = 1:length(coh_set)
    
    % for each condition, the residual activity before motion onset is calculated for baseline adjustment
    BaseAdj_dots(i) = nanmean(FR_grand_dots{i})-nanmean(cell2mat(FR_grand_dots));
    BaseAdj(i) = nanmean(FR_grand{i})-nanmean(cell2mat(FR_grand));
    
    DotsOn{i} = spk_psth_v2(trial_spike_dots_ts,G1address_dots{i},[dots_time1(:)],[pre_dots post_dots],t_array_d,binwidth,step) - BaseAdj_dots(i);
    LibSacc{i} = spk_psth_v2(trial_spike_ts,G1address{i},[sacc1_times(:)],[pre_sacc1 post_sacc1],t_array_s1,binwidth,step) - BaseAdj(i);
    Slib{i} = spk_psth_v2(trial_spike_ts,G1address{i},[slib_times(:)],[pre_slib post_slib],t_array_sl,binwidth,step) - BaseAdj(i);
    LastFix{i} = spk_psth_v2(trial_spike_ts,G1address{i},[fix2_times(:)],[pre_fp2 post_fp2],t_array_fp2,binwidth,step) - BaseAdj(i);
    ChoiceSacc{i} = spk_psth_v2(trial_spike_ts,G1address{i},[choice_sacc_times(:)],[pre_sacc post_sacc],t_array_s,binwidth,step) - BaseAdj(i);
    
end


timex_s = pre_sacc:step:post_sacc;
timex_s1 = pre_sacc1:step:post_sacc1;
timex_sl = pre_slib:step:post_slib;
timex_d = pre_dots:step:post_dots;
timex_fp2 = pre_fp2:step:post_fp2;


maxy1 = max([LibSacc{12} DotsOn{12} ChoiceSacc{12} Slib{12}]);
maxy1 = maxy1+10;

ColorMapCoh = [0 0.5 0; 0 0.7 0.3; 0.4 0.8 0.6; 0.4 0.8 0.6; 0.6 0.9 0.8; 0.8 1 0.8; ...
    1 1 0.8 ;1 1 0.6; 1 0.9 0.5; 1 0.8 0.4; 1 0.6 0; 1 0 0];

figure(10)
clf;


%% dots
ax1 = subplot(151);
title('motion on')
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
ax2=subplot(152);
title('sacc to T^0')
hold on
for i = 1:length(coh_set)
    plot(timex_s1, LibSacc{i}, 'Color', ColorMapCoh(i,:),'LineWidth',2);
end
line([0 0],[-1 maxy1],'Color','k');
hold off
xlim([pre_sacc1 post_sacc1]);
ylim([-1 maxy1]);
set(gca, 'TickDir','out');


%% pursuit
ax3 = subplot(153);
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
ax4 = subplot(154);
title('resumed fix');
hold on
for i = 1:length(coh_set)
    plot(timex_fp2, LastFix{i}, 'Color', ColorMapCoh(i,:),'LineWidth',2);
end
line([0 0],[-1 maxy1],'Color','k');
hold off
xlim([pre_fp2 post_fp2]);
ylim([-1 maxy1]);
set(gca, 'TickDir','out');

%% choice sacc
ax5 = subplot(155);
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

totalcells_used = length(filename_a)+length(filename_b)+length(filename_c)+length(filename_d);

toc
linkaxes([ax1,ax2,ax3,ax4,ax5],'y');


suptitle([groupname,'  n = ',num2str(totalcells_used)])
set(figure(10),'Position',[100 300 1200 300])




end

