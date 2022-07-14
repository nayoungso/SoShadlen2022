function VD_behav_pmf

% Figure 1B (Variable duration PMF)

figure(10); clf;
for monkey_id = 1:2
    
    clear filename
    path = '~/Variable duration task/BehaviorOnly/';
    
    if monkey_id==1
        filename = dir(strcat(path,'Br*behav_me.mat'));
        titlestring = 'Monkey B'
    else
        filename = dir(strcat(path,'H*behav_me.mat'));
        titlestring = 'Monkey H'
    end
    


trialno = 0;

for i = 1:length(filename)    
    clear trials trials_vd
    load([path,filename(i).name]);
    
    for j = 1:length(trials)   
        
        dir_sign = (trials(j).dot_dir ~= trials(j).t1_dir); % 0 for t1_dir (T+) / 1 for t2_dir (T-)
        
        % Motion coherence is convereted to percent(%) from per mille, which is the unit used in Rex 
        coh = (trials(j).dot_coh/10)*(1-2*dir_sign);  % positive coh for T+, negative for T-   
        
        if (trials(j).taskid <= 21)&(trials(j).response>=0)&(trials(j).dot_coh>=0)
            
            trialno = trialno+1;
            iem_coh(trialno) = coh;
            iem_resp(trialno) = abs(dir_sign-trials(j).response);    % 1 for t1 choice, 0 for t2 choice
            
                        
        end
        
    end
    
end


coh_set = unique(iem_coh);

plotx = linspace(coh_set(1)-5,coh_set(end)+5,30);

%logistic regression
b_lr_iem = glmfit(iem_coh',[iem_resp' ones(length(iem_resp),1)],'binomial','link','logit')
ploty_lr_iem = glmval(b_lr_iem,plotx,'logit');

for i = 1:length(coh_set)
    
    idx = (iem_coh==coh_set(i));
    avg_iem_resp(i) = nanmean(iem_resp(idx));
    se_iem_resp(i) = nanstd(iem_resp(idx))/sqrt(sum(double(idx)));
    
end


figure(10)
subplot(1,2,monkey_id)
hold on
plot(plotx,ploty_lr_iem,'-','Color',[0 0.5 0],'LineWidth',2);
e1 = errorbar(coh_set,avg_iem_resp,se_iem_resp,'o','MarkerEdgeColor',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'MarkerSize',6,'capsize',0);
e1.Color = 'r';


ylabel('P_{T+}');
xlabel('Motion strength (%coh)')
set(gca,'TickDir', 'out')
title([titlestring,',  b_{1} = ',num2str(b_lr_iem(2))]);
axis square
xlim([-70 70])


end
