function VD_dec_BU(celltype)

% celltype = 1 (leader) / 2 (supporter)

if celltype==1
    groupname = 'classic_leader';
elseif celltype==2
    groupname = 'classic_supporter';
end


path = ['~/Dropbox/plxdata/RDM_IEM/groups/ranksum_dynamics_based_consec3_p_05_base_adjusted/' groupname '/'];

filename_a = dir(strcat(path,'20/r/','*DSP*.mat'))
filename_b = dir(strcat(path,'20/g/','*DSP*.mat'))
filename_c = dir(strcat(path,'21/r/','*DSP*.mat'))
filename_d = dir(strcat(path,'21/g/','*DSP*.mat'))

% for buildup rate calculation as well as mean FR
initial_t = 160;    % based on the ranksum test (IEM_new_putativeDec_ranksum.m)
end_t = initial_t+200;        %360

spk_window = 20; % binwidth for counting spikes (nonoverlapping bins)
t_array = initial_t:spk_window:end_t;

trialno = 0;
coh_set = [-64 -32 -16 -8 -4 0 4 8 16 32 64];
G1no = zeros(1,length(coh_set));

for p = 1:4
    clear filename
    if p==1
        taskid = 20; filename = filename_a; filepath = [path,'/20/r/'];
    elseif p==2
        taskid = 20; filename = filename_b; filepath = [path,'/20/g/'];
    elseif p==3
        taskid = 21; filename = filename_c; filepath = [path,'/21/r/'];
    elseif p==4
        taskid = 21; filename = filename_d; filepath = [path,'/21/g/'];
    end
    
    
    for i = 1:length(filename)
        
        clear trials
        filename(i).name
        load([filepath,filename(i).name]);
        
        for j = 1:length(trials)
            
            if (trials(j).taskid == taskid)&&(trials(j).response >=0)&&(trials(j).mean_maxFR >= 10)
                
                if mod(p,2)==0
                    trials(j).t1_dir = trials(j).t1_dir+180;   % flip the sign for left/down preferring neurons (p=2 or 4) to be combined with right/up preferring neurons
                end
                
                trialno = trialno+1;
                dspsize = find(DSP011(j,:) > trials(j).time_sacc(1)+300,1);
                trial_spike_ts(trialno,1:dspsize) = DSP011(j,1:dspsize);
                dir_sign = (trials(j).dot_dir ~= trials(j).t1_dir); 
                coh = trials(j).dot_coh*(1-2*dir_sign); % signed coherence
                    
                dots_time1(trialno) = trials(j).time_dots_on(1);
                
                t_array_cnt(trialno,:) = [-spk_window/2 spk_window/2];
                
                indext = find(coh_set == coh/10);
                
                G1no(indext) = G1no(indext)+1;
                G1address{indext}(G1no(indext)) = trialno;
                
                
            end
            
            
        end
        
    end
    
end% end of p-loop

% get the average activity across all coherence (for detrending)
for k = 1:length(t_array)
    clear spk_cnt_array
    [spk_cnt_array, fr_temp] = spk_cnt2(trial_spike_ts,[1:trialno],[dots_time1(:)+t_array(k)],[-spk_window/2 spk_window/2],t_array_cnt);
    spk_cnt_det(k) = nanmean(spk_cnt_array);
end


% for each coherence, get the detrended spike counts at each time point (20ms wide bins)
for i = 1:length(coh_set)
    
    for k = 1:length(t_array)
        [spk_cnt{i,k}, fr_temp] = spk_cnt2(trial_spike_ts,G1address{i},[dots_time1(:)+t_array(k)],[-spk_window/2 spk_window/2],t_array_cnt);
        spk_cnt{i,k} = (spk_cnt{i,k}-spk_cnt_det(k))/(spk_window/1000);
        
    end
    
end


for i = 1:length(coh_set)
    
    clear glm_spk_cnt_array glm_t_array
    
    glm_spk_cnt_array = (spk_cnt{i,1})';
    glm_t_array = zeros(size(spk_cnt{i,1}))';
    
    for k = 2:length(t_array)
        
        glm_spk_cnt_array = vertcat(glm_spk_cnt_array,(spk_cnt{i,k})');
        glm_t_array = vertcat(glm_t_array,((t_array(k)-t_array(1))/1000)*ones(size((spk_cnt{i,k})')));
        
    end
    
    [b{i},dev,stat{i}] = glmfit(glm_t_array,glm_spk_cnt_array,'normal','link','identity');
    
    buildup_b(i,1) = b{i}(2);
    buildup_b_se(i,1) = stat{i}.se(2);
    
end


[bb,~,sstat] = glmfit(coh_set',buildup_b(:,1),'normal','link','identity');

coh_x = -64:1:64;
y_bdp = glmval(bb,coh_x,'identity');


figure()
clf;
errorbar(coh_set,buildup_b,buildup_b_se,'ko','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'capsize',0,'LineWidth',2);
hold on
plot(coh_x,y_bdp,'k-','LineWidth',2);
set(gca, 'TickDir','out');
axis square
box off
xlim([-65 65])
ylim([-20 20])
ylabel('Build up rate (sp/s^2)'); xlabel('Motion strength (% coh)');
axis square

bb
sstat.p
