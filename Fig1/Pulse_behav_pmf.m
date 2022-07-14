function Pulse_behav_pmf(task_var)

% task_var = 1 (1st variant), 2 (2nd variant)
% Figure 1C, 1D 

figure(32); clf;

for monkey_id = 1:2
    
    clear filename
    
    if monkey_id==1
        
        if task_var==1
            path = '~/Two-pulse task, 1st variant/BehaviorOnly/';
            filename = dir(strcat(path,'Br*behav.mat'))
            allo_ind = 0;
            titlestring = 'Monkey B, 1^{st} variant';
        else
            path = '~/Two-pulse task, 2nd variant/BehaviorOnly/';
            filename = dir(strcat(path,'Br*allo*behav.mat'))
            allo_ind =1;
            titlestring = 'Monkey B, 2^{nd} variant';
        end
    else
        
        if task_var==1
            path = '~/Two-pulse task, 1st variant/BehaviorOnly/';
            filename = dir(strcat(path,'H*behav.mat'))
            allo_ind = 0;
            titlestring = 'Monkey H, 1^{st} variant';
        else
            path = '~/Two-pulse task, 2nd variant/BehaviorOnly/';
            filename = dir(strcat(path,'H*allo*behav.mat'))
            allo_ind =1;
            titlestring = 'Monkey H, 2^{nd} variant';
        end
    end
    

dp_trialno = 0;
sp_trialno = 0;

for i = 1:length(filename)
    
    includ_idx = 0;
    
    if (isempty(strfind(filename(i).name,'allo')))&&(allo_ind==0)
        includ_idx = 1;
    elseif (~isempty(strfind(filename(i).name,'allo')))&&(allo_ind==1)
        includ_idx = 1;
    end
    
    if includ_idx
        
        clear trials
        load([path,filename(i).name]);
        
        for j = 1:length(trials)
            
            if ((allo_ind==0)&&(trials(j).taskid <= 31)&&(trials(j).response>=0))||((allo_ind==1)&&((trials(j).taskid == 34)||(trials(j).taskid == 35))&&(trials(j).response>=0))
                
                dp_trialno = dp_trialno+1;
                
                dir_id = 1-floor(mod(trials(j).dot_dir,360)/180)*2; % 1 for T+ / -1 for T-
                coh_dp(dp_trialno,:) = trials(j).dot_coh*dir_id/10;
                
                Rarray_dp(dp_trialno) = ((trials(j).response==1)&&(dir_id==1))||((trials(j).response==0)&&(dir_id==-1));
                
                
                
            elseif ((allo_ind==0)&&(trials(j).taskid >=36)&&(trials(j).response>=0))||((allo_ind==1)&&(trials(j).taskid >=38)&&(trials(j).response>=0))            
                
                sp_trialno = sp_trialno+1;
                
                dir_id = 1-floor(mod(trials(j).dot_dir,360)/180)*2;
                coh_sp(sp_trialno,:) = trials(j).dot_coh*dir_id/10;
                
                Rarray(sp_trialno) = ((trials(j).response==1)&&(dir_id==1))||((trials(j).response==0)&&(dir_id==-1));
                
            end
        end
        
    end
    
end

dp_trialno
sp_trialno

clear trials
pos_coh_set = unique([coh_sp]);
neg_coh_set = sort(pos_coh_set,'descend').*-1;
coh_set = unique([neg_coh_set pos_coh_set]);

coh_first = coh_dp(:,1);
coh_second = coh_dp(:,2);
coh_comb = sum(coh_dp,2)/2;
coh_max = max(coh_first,coh_second);
coh_min = min(coh_first,coh_second);


pos_comb_coh_set = unique(coh_comb);
neg_comb_coh_set = sort(pos_comb_coh_set,'descend').*-1;
comb_coh_set = unique([neg_comb_coh_set pos_comb_coh_set]);


trialno_comb = zeros(1,length(comb_coh_set));
trialno_sp = zeros(1,length(coh_set));

for i = 1:dp_trialno
    
    coh_id_comb = find(comb_coh_set == sum(coh_dp(i,:))/2);
    trialno_comb(coh_id_comb) = trialno_comb(coh_id_comb)+1;
    Rarray_comb{coh_id_comb,trialno_comb(coh_id_comb)} = Rarray_dp(i);
    
end

for i = 1:sp_trialno
    
    coh_id_sp = find(coh_set == coh_sp(i));
    trialno_sp(coh_id_sp) = trialno_sp(coh_id_sp)+1;
    Rarray_sp{coh_id_sp, trialno_sp(coh_id_sp)} = Rarray(i);
    
end

for i = 1:length(coh_set)

    Pr_mean_sp(i) = nanmean([Rarray_sp{i,:}]);
    Pr_se_sp(i) = nanstd([Rarray_sp{i,:}])./sqrt(trialno_sp(i));
    
end

for i = 1:length(comb_coh_set)
    
    if trialno_comb(i) > 0
        Pr_mean_comb(i) = nanmean([Rarray_comb{i,:}]);
        Pr_se_comb(i) = nanstd([Rarray_comb{i,:}])./sqrt(trialno_comb(i));
    else
        Pr_mean_comb(i) = nan; Pr_se_comb(i) = nan;
    end
    
    
end


coh_plot = -65:65;
comb_coh_plot = coh_plot;


Rarray_dp_fit = [Rarray_dp' ones(length(Rarray_dp),1)];
size(Rarray_dp_fit)

Rarray_sp_fit = [Rarray' ones(length(Rarray),1)];
size(Rarray_sp_fit)

% fitting
[b_comb,~,stat_comb] = glmfit(coh_comb,Rarray_dp_fit,'binomial','link','logit');        % eq. 2
Pr_plot_comb = glmval(b_comb,comb_coh_plot,'logit');

[b_sp,~,stat_sp] = glmfit(coh_sp,Rarray_sp_fit,'binomial','link','logit');          % eq. 1 (single pulse control trials)
Pr_plot_sp = glmval(b_sp,coh_plot,'logit');


[b_comb2,~,stat2] = glmfit([coh_first, coh_second],Rarray_dp_fit,'binomial','link','logit');    % eq. 15
[b_max_min,~,stat_max_min] = glmfit([coh_max coh_min],Rarray_dp_fit,'binomial','link','logit'); % eq. 16
text1 = ([num2str(b_comb2(2)),' P1 + ',num2str(b_comb2(3)),' P2']);
text2 = ([num2str(b_max_min(2)),' P_{max} + ',num2str(b_max_min(3)),' P_{min}']);


title_addstr = [titlestring ' ;'];




subplot(2,3,3*(monkey_id-1)+1);
hold on
plot(coh_plot,Pr_plot_sp,'Color',[0.5 0.5 0.5],'LineWidth',1.2);
plot(comb_coh_plot,Pr_plot_comb,'k-','LineWidth',1.2);
errorbar(coh_set,Pr_mean_sp(:),Pr_se_sp(:),'o','LineWidth',1,'Color',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',6,'capsize',0);
errorbar(comb_coh_set,Pr_mean_comb(:),Pr_se_comb(:),'ko','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6,'capsize',0);
xlim([coh_set(1)-1 coh_set(end)+1]);
ylim([0 1]);
grid on
ylabel('P_{T+}'); xlabel('Average motion strength (%coh)');
title(strcat(title_addstr, ' SPslope = ', num2str(b_sp(2)), '  DPslope = ', num2str(b_comb(2))));
hold off
box off
set(gca, 'TickDir','out');


subplot(2,3,[3*(monkey_id-1)+2 3*(monkey_id-1)+3]);
xlim([0 10]); ylim([0 5]);

text(1,3.3,text1,'FontSize',16);
text(1,2.3,text2,'FontSize',16);


set(gca,'visible','off')
set(gca,'xtick',[]); set(gca,'ytick',[]);


end
end
