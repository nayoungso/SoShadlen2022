%function Pulse_behav_modelSimul_pmf_FigS1

clear all;

monkey_id = 'H';   % {'B','H'} for monkey B, H
task_var_id = 2;    % [1, 2] for 1st / 2nd variant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coh_set = [0 4 8 16 32 64]';
dir_set = [-1 1]';
coh_plot = -65:65;

simulno = 10;


%% actual data
path = ['~/Dropbox/plxdata/Pulse/Vprobe/Behav/'];

if task_var_id==1
    %filename = dir(strcat(path,'*',monkey_id,'*Pulse_behav.mat'));
    filename = dir(strcat(path,monkey_id,'*Pulse_behav.mat'));
elseif task_var_id==2
    %filename = dir(strcat(path,'*',monkey_id,'*allo*Pulse_behav.mat'));
    filename = dir(strcat(path,monkey_id,'*allo*Pulse_behav.mat'));
end

dp_trialno = 0;
for i = 1:length(filename)
    
    clear trials
    filename(i).name
    load([path,filename(i).name],'trials');
     
    for j = 1:length(trials)        

        if ((task_var_id==1)&(trials(j).taskid <= 31))|((task_var_id==2)&(trials(j).taskid <= 35)&(trials(j).taskid >= 34))  
        %if (trials(j).taskid <= 31)&&(trials(j).response>=0)
        
        if trials(j).response>=0
            dp_trialno = dp_trialno+1;
            
            dir_id = 1-floor(mod(trials(j).dot_dir,360)/180)*2;
            
            dir_id_dp(dp_trialno) = dir_id;
            coh_dp(dp_trialno,:) = trials(j).dot_coh*dir_id/10; % signed coh in units of % 
            coh_avg(dp_trialno) = sum(coh_dp(dp_trialno,:))/2;
            
            Rarray_dp(dp_trialno) = ((trials(j).response==1)&&(dir_id==1))||((trials(j).response==0)&&(dir_id==-1));
       
        end
        end
    end
    
end

clear trials

coh_avg_set = unique(coh_avg);

trialno_comb = zeros(1,length(coh_avg_set));

for i = 1:dp_trialno
    
    coh_id = find(coh_avg_set == sum(coh_dp(i,:))/2);
    trialno_comb(coh_id) = trialno_comb(coh_id)+1;
    
    Rarray_data{coh_id,trialno_comb(coh_id)} = Rarray_dp(i);
    
end

for i = 1:length(coh_avg_set)
    
    if trialno_comb(i) > 0
        Pr_mean_data(i) = nanmean([Rarray_data{i,:}]);
        Pr_se_data(i) = nanstd([Rarray_data{i,:}])./sqrt(trialno_comb(i));
    else
        Pr_mean_data(i) = nan; Pr_se_data(i) = nan;
    end
    
end

[b_data,~,~] = glmfit(coh_avg',[Rarray_dp' ones(length(Rarray_dp),1)],'binomial','link','logit');
Pr_plot_data = glmval(b_data,coh_plot,'logit');

b_data

trialno = dp_trialno;   % matching the simulated trial number to the one from the actual data



%% simulated data

b0 = b_data(1);
b1_array = b_data(2);


for k = 1:length(b1_array)
    
    b1 = b1_array(k);
    b = [b0 b1];
    
    for j = 1:simulno
        
        %% generating trials
        dir_array = randsample(dir_set,trialno,true);
        coh_array(:,1) = randsample(coh_set,trialno,true);  % P1 coh
        coh_array(:,2) = randsample(coh_set,trialno,true);  % P2 coh
        maxcoh_array = max(coh_array,[],2);                 % C_stronger (stronger motrion coh between P1 and P2)
        mincoh_array = min(coh_array,[],2);                 % C_weaker
        random_idx = randsample(2,trialno,true);            % one randomly picked pulse for each trial
        
        signed_coh_array = dir_array.*coh_array;
        maxcoh_array = dir_array.*maxcoh_array;
        mincoh_array = dir_array.*mincoh_array;

        p1p2_comb_array = (signed_coh_array(:,1) + signed_coh_array(:,2))/2;    % average motion strength of the two pulses--using avg coh to compare with the data (Figure 1C,D) 
        
        avg_coh_set = unique(p1p2_comb_array);
        trialno_comb = zeros(1,length(avg_coh_set));
        
        AllPulses_cp = calculate_cp(p1p2_comb_array,b0,b1);   % probability of choosing T+, if using both pulses (Model-2 in Supplementary Information)
         
        for i = 1:trialno
            
            randomOne_cp(i) = calculate_cp(2*signed_coh_array(i,random_idx(i)),b0,b1); % probability of choosing T+ if using one randomly selected pulse (Model-1 in Supplementary Information)   
            
            AllPulses_choice(i) = random('binomial',1,AllPulses_cp(i));     % generating choices based on Model-2
            randomOne_choice(i) = random('binomial',1,randomOne_cp(i));     % generating choices based on Model-1
            
            coh_id = find(avg_coh_set == p1p2_comb_array(i));
            trialno_comb(coh_id) = trialno_comb(coh_id)+1;
           
            Rarray_comb_rand{coh_id,trialno_comb(coh_id)} = randomOne_choice(i);
            Rarray_comb_all{coh_id,trialno_comb(coh_id)} = AllPulses_choice(i);
            
        end
        AllPulses_choice = AllPulses_choice(:);
        randomOne_choice = randomOne_choice(:);
        
        %% fitting each simulation to logistic functions (eq. 15, 16)
        
        % Model-1 (randomly selected one pulse)
        [b_p1p2_rand,~,~] = glmfit([signed_coh_array(:,1) signed_coh_array(:,2)],[randomOne_choice ones(trialno,1)],'binomial','link','logit');
        [b_pmaxpmin_rand,~,~] = glmfit([maxcoh_array mincoh_array],[randomOne_choice ones(trialno,1)],'binomial','link','logit');
        
        % storing coefficients from each iteration of simulation
        bp1_rand(k,j) = b_p1p2_rand(2); bp2_rand(k,j) = b_p1p2_rand(3);
        bmax_rand(k,j) = b_pmaxpmin_rand(2); bmin_rand(k,j) = b_pmaxpmin_rand(3);   
        
        % fitting to eq.2, to compare with the actual data fit in Figure 1C,D
        [b_pavg_rand,~,~] = glmfit([p1p2_comb_array],[randomOne_choice ones(trialno,1)],'binomial','link','logit');
        model1_b1(j) = b_pavg_rand(2);
        
        % Model-2 (both pulses)
        [b_p1p2_sum,~,~] = glmfit([signed_coh_array(:,1) signed_coh_array(:,2)],[AllPulses_choice ones(trialno,1)],'binomial','link','logit');
        [b_pmaxpmin_sum,~,~] = glmfit([maxcoh_array mincoh_array],[AllPulses_choice ones(trialno,1)],'binomial','link','logit');
        
        % storing coefficients from each iteration of simulation
        bp1_sum(k,j) = b_p1p2_sum(2); bp2_sum(k,j) = b_p1p2_sum(3);
        bmax_sum(k,j) = b_pmaxpmin_sum(2); bmin_sum(k,j) = b_pmaxpmin_sum(3);
        
        [b_pavg_sum,~,~] = glmfit([p1p2_comb_array],[AllPulses_choice ones(trialno,1)],'binomial','link','logit');
        model2_b1(j) = b_pavg_sum(2);
        
        % to generate the points for the pmf figure (Fig S1C)
        for i = 1:length(avg_coh_set)
            
            Pr_mean_rand(i) = nanmean([Rarray_comb_rand{i,:}]);
            Pr_se_rand(i) = nanstd([Rarray_comb_rand{i,:}])./sqrt(trialno_comb(i));
            
            Pr_mean_all(i) = nanmean([Rarray_comb_all{i,:}]);
            Pr_se_all(i) = nanstd([Rarray_comb_all{i,:}])./sqrt(trialno_comb(i));
            
        end
        
        
        
        
        
    end
end


[b_comb_rand,~,~] = glmfit(p1p2_comb_array,[randomOne_choice ones(length(randomOne_choice),1)],'binomial','link','logit');
Pr_plot_comb_rand = glmval(b_comb_rand,coh_plot,'logit');

[b_comb_all,~,~] = glmfit(p1p2_comb_array,[AllPulses_choice ones(length(AllPulses_choice),1)],'binomial','link','logit');
Pr_plot_comb_all = glmval(b_comb_all,coh_plot,'logit');

toc

figure(13)
clf;

plot(coh_plot,Pr_plot_comb_rand,'-','Color',[0.3010 0.7450 0.9330],'LineWidth',0.8);
hold on;
errorbar(avg_coh_set,Pr_mean_rand(:),Pr_se_rand(:),'o','Color',[0.3010 0.7450 0.9330],'LineWidth',1.5,'MarkerFaceColor',[0.3010 0.7450 0.9330],'MarkerSize',10,'capsize',0);

plot(coh_plot,Pr_plot_comb_all,'k-','LineWidth',0.8);
errorbar(avg_coh_set,Pr_mean_all(:),Pr_se_all(:),'ko','LineWidth',1.5,'MarkerFaceColor','k','MarkerSize',10,'capsize',0);

plot(coh_plot,Pr_plot_data,'r-','LineWidth',0.8);
errorbar(avg_coh_set,Pr_mean_data(:),Pr_se_data(:),'ro','LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',10,'capsize',0);

xlim([coh_plot(1)-1 coh_plot(end)+1]);
ylim([0 1]);
hold off; box off; axis square;
ylabel('P_R'); xlabel('Average motion strength (%coh)');
set(gca, 'TickDir','out');


    function cp = calculate_cp(motioncoh,b0,b1)
        cp = (1+exp(-(b0+b1*motioncoh))).^-1;
    end
 