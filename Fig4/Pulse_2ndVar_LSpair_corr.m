function Pulse_2ndVar_LSpair_corr(phi)

% phi=0.6 to create Figure 4B 

simpair_path = ['~/Two-pulse task, 2nd variant/SimPair/']
simpair_ref_file = [simpair_path,'Kendall_based_simpair_allo_p_05_consec3.mat'];

load(simpair_ref_file) % information for each pair (e.g., data path, the time when the supporter starts to show decision-related activity

binwidth = 300; 
coh_set = [-640 -320 -160 -80 -40 0 40 80 160 320 640];
pairno = 0;
idstart = 0;

permno = 1000;   % for the null distribution, random permutation number for IEM supporter surrogate


for p = 1:4
    
    clear Leader Supporter
    if p==1
        taskid = 34; Leader = RRpairLeaderA; Supporter = RRpairFollowerA;
    elseif p==2
        taskid = 34; Leader = GGpairLeaderA; Supporter = GGpairFollowerA;
    elseif p==3
        taskid = 35; Leader = RRpairLeaderB; Supporter = RRpairFollowerB;
    elseif p==4
        taskid = 35; Leader = GGpairLeaderB; Supporter = GGpairFollowerB;
    end
    
    for i = 2:size(Leader,1)
        
        % resetting the arrays for the next pair
        clear trials leader_times supporter_times
        clear t_array_leader t_array_supporter
        clear spk_N1 spk_N2 leader_trials supporter_trials Gaddress
        clear res_P1_leader_temp res_P1_supporter_temp res_IEM_leader_temp res_IEM_supporter_temp res_FP2_leader_temp res_FP2_supporter_temp
        clear mean_P1_leader_temp mean_P1_supporter_temp mean_IEM_leader_temp mean_IEM_supporter_temp mean_FP2_leader_temp mean_FP2_supporter_temp
        clear raw_P1_leader_temp raw_P1_supporter_temp raw_IEM_leader_temp raw_IEM_supporter_temp raw_FP2_leader_temp raw_FP2_supporter_temp
        clear trialid_P1_leader_temp trialid_P1_supporter_temp trialid_IEM_leader_temp trialid_IEM_supporter_temp trialid_FP2_leader_temp trialid_FP2_supporter_temp
        clear res_IEM_supporter_shuf_temp mean_IEM_supporter_shuf_temp raw_IEM_supporter_shuf_temp trialid_IEM_supporter_shuf_temp
        clear res_IEM_leader_shuf_temp mean_IEM_leader_shuf_temp raw_IEM_leader_shuf_temp trialid_IEM_leader_shuf_temp
        
        load(Leader{i,6})
        
        pairno = pairno+1;
        
        spk_leader = DSP011;
        leader_peaktime = 200;
        
        
        leader_trials = trials;
        clear DSP011 trials
        
        load(Supporter{i,6})
        spk_supporter = DSP011;
        supporter_peaktime = Supporter{i,4};  % collecting the time the supporter starts to show decision-related activity
        supporter_alignment = Supporter{i,5};
        supporter_trials = trials;
        clear DPS011
        trialno = 0;
        Gno = zeros(1,length(coh_set));
        
        for j = 1:length(trials)
            
            if mod(p,2)==0
                trials(j).t1_dir = trials(j).t1_dir+180;
            end
            
            if (trials(j).taskid == taskid)&&(trials(j).response>=0)
                
                dir_sign = (trials(j).dot_dir ~= trials(j).t1_dir); 
                coh = trials(j).dot_coh(1)*(1-2*dir_sign);  % signed coh
                
                if ((dir_sign == 0)&&(trials(j).response == 1))||((dir_sign == 1)&&(trials(j).response == 0))
                    choice = 2;     % T+ choice
                else
                    choice = 1;     % T- choice
                end
                
                if (leader_trials(j).mean_maxFR >= 10)&&(supporter_trials(j).mean_maxFR >= 10)
                    
                    trialno = trialno+1;
                    dspsize_n1 = find(spk_leader(j,:) > trials(j).time_sacc(end)+300,1);
                    dspsize_n2 = find(spk_supporter(j,:) > trials(j).time_sacc(end)+300,1);
                    spk_N1(trialno,1:dspsize_n1) = spk_leader(j,1:dspsize_n1);
                    spk_N2(trialno,1:dspsize_n2) = spk_supporter(j,1:dspsize_n2);
                    
                    leader_times(trialno) = trials(j).time_dots_on(1)+leader_peaktime;
                    t_array_leader(trialno,:) = [0 binwidth]; % l_{pre} is measured 200?500 from P1 onset   
                    
                    if strcmp(supporter_alignment,'Pursuit')
                        supporter_times(trialno) = trials(j).time_slib_start+supporter_peaktime;
                    elseif strcmp(supporter_alignment, 'FP2')
                        supporter_times(trialno) = trials(j).time_fp_acq(2)+supporter_peaktime;
                    elseif strcmp(supporter_alignment, 'Pulse2')
                        supporter_times(trialno) = trials(j).time_dots_on(2)+supporter_peaktime;                    
                    end
                    t_array_supporter(trialno,:) = [-binwidth/2 binwidth/2];  % s_{IEM} is measured from 300ms around its decision-related activity
                    
                    coh_cell{pairno}(trialno,1) = coh;
                    choice_cell{pairno}(trialno,1) = choice;
                    
                    indext = find(coh_set==coh_cell{pairno}(trialno));
                    Gno(indext) = Gno(indext)+1;
                    Gaddress{indext}(Gno(indext)) = trialno;     
                    
                end
            end
        end
        
    
        if trialno < 50   % excluding the pairs with 50 trials or less 
            
            coh_cell{pairno} = nan; choice_cell{pairno} = nan;
            res_P1_leader_temp = nan; mean_P1_leader_temp = nan; trialid_P1_leader_temp = nan; raw_P1_leader_temp = nan;
            res_P1_supporter_temp = nan; mean_P1_supporter_temp = nan; trialid_P1_supporter_temp = nan; raw_P1_supporter_temp = nan;
            res_IEM_leader_temp = nan; mean_IEM_leader_temp = nan; trialid_IEM_leader_temp = nan; raw_IEM_leader_temp = nan;
            res_IEM_supporter_temp = nan; mean_IEM_supporter_temp = nan; trialid_IEM_supporter_temp = nan; raw_IEM_supporter_temp = nan;
            
            res_IEM_supporter_shuf_temp = nan(1,permno); mean_IEM_supporter_shuf_temp = nan(1,permno); trialid_IEM_supporter_shuf_temp = nan(1,permno); raw_IEM_supporter_shuf_temp = nan(1,permno);
            res_IEM_leader_shuf_temp = nan(1,permno); mean_IEM_leader_shuf_temp = nan(1,permno); trialid_IEM_leader_shuf_temp = nan(1,permno); raw_IEM_leader_shuf_temp = nan(1,permno);
            
        else
            
            % collecting the activity of leader and supporter during the pre-IEM epoch
            [res_P1_leader_temp, mean_P1_leader_temp, trialid_P1_leader_temp, raw_P1_leader_temp] = spk_varce4(spk_N1,Gaddress,[leader_times(:)],t_array_leader(1,:),idstart);
            [res_P1_supporter_temp, mean_P1_supporter_temp, trialid_P1_supporter_temp, raw_P1_supporter_temp] = spk_varce4(spk_N2,Gaddress,[leader_times(:)],t_array_leader(1,:),idstart);
    
            % IEM epoch
            [res_IEM_leader_temp, mean_IEM_leader_temp, trialid_IEM_leader_temp, raw_IEM_leader_temp] = spk_varce4(spk_N1,Gaddress,[supporter_times(:)],t_array_supporter(1,:),idstart);
            [res_IEM_supporter_temp, mean_IEM_supporter_temp, trialid_IEM_supporter_temp, raw_IEM_supporter_temp] = spk_varce4(spk_N2,Gaddress,[supporter_times(:)],t_array_supporter(1,:),idstart);
    
            % (control) shuffling the trial orders (while preserving the coherence) for the activity during the IEM epoch
            [res_IEM_supporter_shuf_temp, mean_IEM_supporter_shuf_temp, trialid_IEM_supporter_shuf_temp, raw_IEM_supporter_shuf_temp] = spk_varce4_randperm(spk_N2,Gaddress,[supporter_times(:)],t_array_supporter(1,:),idstart,permno);
            [res_IEM_leader_shuf_temp, mean_IEM_leader_shuf_temp, trialid_IEM_leader_shuf_temp, raw_IEM_leader_shuf_temp] = spk_varce4_randperm(spk_N1,Gaddress,[supporter_times(:)],t_array_supporter(1,:),idstart,permno);
    
            idstart = trialid_P1_leader_temp(end);  % saving the collected trial numbers to concatenate data from all pairs 
            
        end
        
        % concatenating
        if pairno==1 % first pair
            
            coh_all = coh_cell{pairno}; choice_all = choice_cell{pairno};
            
            res_P1_leader_all = res_P1_leader_temp; res_P1_supporter_all = res_P1_supporter_temp;
            res_IEM_leader_all = res_IEM_leader_temp; res_IEM_supporter_all = res_IEM_supporter_temp;
            
            mean_P1_leader_all = mean_P1_leader_temp; mean_P1_supporter_all = mean_P1_supporter_temp;
            mean_IEM_leader_all = mean_IEM_leader_temp; mean_IEM_supporter_all = mean_IEM_supporter_temp;
            
            raw_P1_leader_all = raw_P1_leader_temp; raw_P1_supporter_all = raw_P1_supporter_temp;
            raw_IEM_leader_all = raw_IEM_leader_temp; raw_IEM_supporter_all = raw_IEM_supporter_temp;
            
            res_IEM_supporter_shuf = res_IEM_supporter_shuf_temp; raw_IEM_supporter_shuf = raw_IEM_supporter_shuf_temp;
            mean_IEM_supporter_shuf = mean_IEM_supporter_shuf_temp; 
            
            res_IEM_leader_shuf = res_IEM_leader_shuf_temp; raw_IEM_leader_shuf = raw_IEM_leader_shuf_temp;
            mean_IEM_leader_shuf = mean_IEM_leader_shuf_temp; 
            
        else
            coh_all = [coh_all; coh_cell{pairno}]; choice_all = [choice_all; choice_cell{pairno}];
            res_P1_leader_all = [res_P1_leader_all; res_P1_leader_temp]; res_P1_supporter_all = [res_P1_supporter_all; res_P1_supporter_temp];
            res_IEM_leader_all = [res_IEM_leader_all; res_IEM_leader_temp]; res_IEM_supporter_all = [res_IEM_supporter_all; res_IEM_supporter_temp];
            
            mean_P1_leader_all = [mean_P1_leader_all; mean_P1_leader_temp]; mean_P1_supporter_all = [mean_P1_supporter_all; mean_P1_supporter_temp];
            mean_IEM_leader_all = [mean_IEM_leader_all; mean_IEM_leader_temp]; mean_IEM_supporter_all = [mean_IEM_supporter_all; mean_IEM_supporter_temp];
            
            raw_P1_leader_all = [raw_P1_leader_all; raw_P1_leader_temp]; raw_P1_supporter_all = [raw_P1_supporter_all; raw_P1_supporter_temp];
            raw_IEM_leader_all = [raw_IEM_leader_all; raw_IEM_leader_temp]; raw_IEM_supporter_all = [raw_IEM_supporter_all; raw_IEM_supporter_temp];
            
            res_IEM_supporter_shuf = [res_IEM_supporter_shuf; res_IEM_supporter_shuf_temp]; 
            raw_IEM_supporter_shuf = [raw_IEM_supporter_shuf; raw_IEM_supporter_shuf_temp];
            mean_IEM_supporter_shuf = [mean_IEM_supporter_shuf; mean_IEM_supporter_shuf_temp]; 
            
            res_IEM_leader_shuf = [res_IEM_leader_shuf; res_IEM_leader_shuf_temp]; 
            raw_IEM_leader_shuf = [raw_IEM_leader_shuf; raw_IEM_leader_shuf_temp];
            mean_IEM_leader_shuf = [mean_IEM_leader_shuf; mean_IEM_leader_shuf_temp]; 
            
        end
    
    end % end of cell-specific part
    
    
end % end of p-loop (different types of pair, based on RF location and target configurations in which they are identified)


% get rid of nan
nanind = [];
nanind = [nanind; find(isnan(res_P1_leader_all))]; nanind = [nanind; find(isnan(res_P1_supporter_all))];
nanind = [nanind; find(isnan(res_IEM_leader_all))]; nanind = [nanind; find(isnan(res_IEM_supporter_all))];
nanind = unique(nanind);

coh_all(nanind) = []; choice_all(nanind) = [];
res_P1_leader_all(nanind) = []; res_P1_supporter_all(nanind) = [];
res_IEM_leader_all(nanind) = []; res_IEM_supporter_all(nanind) = [];
raw_P1_leader_all(nanind) = []; raw_P1_supporter_all(nanind) = [];
raw_IEM_leader_all(nanind) = []; raw_IEM_supporter_all(nanind) = [];
mean_P1_leader_all(nanind) = []; mean_P1_supporter_all(nanind) = [];
mean_IEM_leader_all(nanind) = []; mean_IEM_supporter_all(nanind) = [];
raw_IEM_supporter_shuf(nanind,:) = []; res_IEM_supporter_shuf(nanind,:) = [];
mean_IEM_supporter_shuf(nanind,:) = []; 
raw_IEM_leader_shuf(nanind,:) = []; res_IEM_leader_shuf(nanind,:) = [];
mean_IEM_leader_shuf(nanind,:) = []; 

% VarCE calculation (eq. 7 & 8)
resid_P1_leader = raw_P1_leader_all - mean_P1_leader_all;
VarCE_P1_leader = nanmean(resid_P1_leader.^2 - phi*mean_P1_leader_all);

resid_IEM_supporter = raw_IEM_supporter_all - mean_IEM_supporter_all;
VarCE_IEM_supporter = nanmean(resid_IEM_supporter.^2 - phi*mean_IEM_supporter_all);


% uninformative epochs
resid_P1_supporter = raw_P1_supporter_all - mean_P1_supporter_all;
VarCE_P1_supporter = nanmean(resid_P1_supporter.^2 - phi*mean_P1_supporter_all);

resid_IEM_leader = raw_IEM_leader_all - mean_IEM_leader_all;
VarCE_IEM_leader = nanmean(resid_IEM_leader.^2 - phi*mean_IEM_leader_all);


% CorCE calculation (eq. 9)
CorCE_P1_IEM_temp = nancov([resid_P1_leader resid_IEM_supporter],'pairwise');
CorCE_P1_IEM = CorCE_P1_IEM_temp(1,2)/sqrt(VarCE_P1_leader*VarCE_IEM_supporter);

CorCE_P1_IEM_cont_temp = nancov([resid_P1_supporter resid_IEM_leader],'pairwise');
CorCE_P1_IEM_cont = CorCE_P1_IEM_cont_temp(1,2)/sqrt(VarCE_P1_supporter*VarCE_IEM_leader);


% calculating correlation with the shuffled control 
for i = 1:permno
    resid_IEM_supporter_shuf(:,i) = raw_IEM_supporter_shuf(:,i) - mean_IEM_supporter_shuf(:,i);
    VarCE_IEM_supporter_shuf(i) = nanmean(resid_IEM_supporter_shuf(:,i).^2 - phi*mean_IEM_supporter_shuf(:,i));

    resid_IEM_leader_shuf(:,i) = raw_IEM_leader_shuf(:,i) - mean_IEM_leader_shuf(:,i);
    VarCE_IEM_leader_shuf(i) = nanmean(resid_IEM_leader_shuf(:,i).^2 - phi*mean_IEM_leader_shuf(:,i));
    
    CorCE_P1_IEM_null_temp = nancov([resid_P1_leader resid_IEM_supporter_shuf(:,i)],'pairwise');
    CorCE_P1_IEM_null(i) = CorCE_P1_IEM_null_temp(1,2)/sqrt(VarCE_P1_leader*VarCE_IEM_supporter_shuf(i));
    
    CorCE_P1_IEM_cont_null_temp = nancov([resid_P1_supporter resid_IEM_leader_shuf(:,i)],'pairwise');
    CorCE_P1_IEM_cont_null(i) = CorCE_P1_IEM_cont_null_temp(1,2)/sqrt(VarCE_P1_supporter*VarCE_IEM_leader_shuf(i));
    
end

CorCE_P1_IEM_null_mean = mean(CorCE_P1_IEM_null)
CorCE_P1_IEM_null_up95 = mean(CorCE_P1_IEM_null) + 2*std(CorCE_P1_IEM_null)

CorCE_P1_IEM_cont_null_mean = mean(CorCE_P1_IEM_cont_null)
CorCE_P1_IEM_cont_null_up95 = mean(CorCE_P1_IEM_cont_null) + 2*std(CorCE_P1_IEM_cont_null)

figure(7)
clf;

labelx{1} = 'l^{pre} & s^{IEM}'; labelx{2} = 's^{pre} & l^{IEM}'; 
line([0.5 2.5],[0 0],'Color','k');

for i = 1:2
    x = [i-0.2 i+0.2];
    if i==1
        y = [CorCE_P1_IEM CorCE_P1_IEM];
        ynull = [CorCE_P1_IEM_null_mean CorCE_P1_IEM_null_mean];
        err = [CorCE_P1_IEM_null_up95 CorCE_P1_IEM_null_up95];
    else
        y = [CorCE_P1_IEM_cont CorCE_P1_IEM_cont];
        ynull = [CorCE_P1_IEM_cont_null_mean CorCE_P1_IEM_cont_null_mean];
        err = [CorCE_P1_IEM_cont_null_up95 CorCE_P1_IEM_cont_null_up95];
    end
    
    hold on;
    shadedErrorBar(x,ynull,err,'y-');
    plot(x,y,'LineWidth',2);
end

xticks(1:2); xticklabels(labelx); 
ylabel('Correlation')
set(gca,'TickDir','out','FontSize',14); box off;
hold off;
ylim([-0.02 0.04])

% regression (eq.10)
[b2, ~, stat2] = glmfit([coh_all choice_all res_P1_leader_all], res_IEM_supporter_all,'normal','link','identity');
b2'
stat2.p'

% regression between the activities during the uninformative epochs
[bnull, ~, statnull] = glmfit([coh_all choice_all res_P1_supporter_all], res_IEM_leader_all,'normal','link','identity');
bnull'
statnull.p'
        
        
        
        
        
        
        
        
        
        
        
        
        
        
