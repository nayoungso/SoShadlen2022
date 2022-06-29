function VD_behav_ME0_FigS1

figure(2); clf;
for monkey_id = 1:2
    
    clear filename
    if monkey_id==1
        path = ['~/Dropbox/plxdata/RDM_IEM/BrBehav/ME/'];
        titlestring = 'Monkey B'
    else
        path = ['~/Dropbox/plxdata/RDM_IEM/HBehav/ME/'];
        titlestring = 'Monkey H'
    end
    
filename = dir(strcat(path,'*behav_me.mat'))

clear trials_temp
trialno = 0;


for i = 1:length(filename)
    
    clear trials
    %filename(i).name
    load([path,filename(i).name]);
    
    for j = 1:length(trials)
        
        if (trials(j).taskid <= 21)&(trials(j).response>=0)&(trials(j).dot_coh == 0)&~isempty(trials(j).cME)
            trialno = trialno+1;
            trials_temp(trialno) = trials(j);            
        end
        
    end
    
end

clear trials 
trials = trials_temp;
clear trials_temp 
length(trials)

for i = 1:length(trials)
    
    
    dir_sign = 1-floor(mod(trials(i).dot_dir,360)/180);  % 1 for right/up motion (0 or 90 deg); 0 for left/down motion (180 or 270 deg) 
    if (trials(i).t1_dir == 270)||(trials(i).t1_dir == -90)     % for some 1ch recording session where I set the t1_dir to downward
        trials(i).mME{1} = -1*trials(i).mME{1};
    end
    
    if ((trials(i).response == 1)&&(dir_sign == 1))||((trials(i).response == 0)&&(dir_sign == 0))       % T+ choice
        
        ME_pref(i,1:length(trials(i).mME{1})) = trials(i).mME{1};
        
    else  % ME sign flipped for T- choice trials s.t. positive ME means the direction that's consistent with the monkey's choice
        
        ME_pref(i,1:length(trials(i).mME{1})) = -1*trials(i).mME{1};
        
    end
        
    
end

ME_pref_avg = nanmean(ME_pref,1);
ME_pref_se = nanstd(ME_pref,1)./sqrt(sum(~isnan(ME_pref),1));

timex = 0:1000/75:(length(ME_pref_avg)-1)*1000/75;



figure(2)

subplot(2,1,monkey_id)
e10 = shadedErrorBar(timex,ME_pref_avg,ME_pref_se,{'-','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5]});
e10.Color = [0.5 0.5 0.5];
hold on
line([0 700],[0 0],'Color','k','LineWidth',1);
hold off
set(gca,'TickDir', 'out')

xlim([0 700]);
ylim([-3 6]);
xlabel('Time from motion on (ms)');
ylabel('Motion energy (a.u.)')
box off
title(titlestring)

end

end






    
    
    
    
    
    