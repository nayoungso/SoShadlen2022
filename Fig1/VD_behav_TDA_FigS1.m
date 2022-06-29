function VD_behav_TDA_FigS1


monkey_id = 1;   %1 for monkey B, 2 for monkey H
    
clear filename
if monkey_id==1
    path = ['~/Dropbox/plxdata/RDM_IEM/BrBehav/ME/'];
    titlestring = 'Monkey B'
else
    path = ['~/Dropbox/plxdata/RDM_IEM/HBehav/ME/'];
    titlestring = 'Monkey H'
end
    
filename = dir(strcat(path,'*behav_me.mat'))

trialno = 0;

for i = 1:length(filename)
    
    clear trials
    load([path,filename(i).name]);
    
    for j = 1:length(trials)
        
        
        if (trials(j).taskid <= 21)&(trials(j).response>=0)
            
            trialno = trialno+1;
            trials_temp(trialno) = trials(j);
            dotdur(trialno) = trials_temp(trialno).time_dots_off - trials_temp(trialno).time_dots_on;
            corr(trialno) = trials(j).response;
            
            dir_id = 1-floor(mod(trials(j).dot_dir,360)/180)*2;
            Rresponse(trialno) = ((trials(j).response==1)&&(dir_id==1))||((trials(j).response==0)&&(dir_id==-1));
            coh_info(trialno) = trials(j).dot_coh*dir_id;
            unsign_coh(trialno) = trials(j).dot_coh;
                        
        end
        
    end
    
end

clear trials 
trials = trials_temp;
clear trials_temp 

pos_coh_set = unique([trials(:).dot_coh]);
neg_coh_set = sort(pos_coh_set,'descend').*-1;
coh_set = unique([neg_coh_set pos_coh_set]);
coh_set = coh_set';

pos_coh_set = nonzeros(pos_coh_set);    

%% for running mean of accuracy (nonzero coherence trials only)
dotdur_steps = 15; % boxcar width
dotdur_matrix = 100:dotdur_steps:550;

for i = 1:length(pos_coh_set)
    for j = 1:length(dotdur_matrix)-1
        idx = (dotdur >= dotdur_matrix(j))&(dotdur < dotdur_matrix(j+1))&(unsign_coh==pos_coh_set(i));
        running_mean(i,j) = nanmean(corr(idx));
    end
end
        

binno = 10; % number of quantiles to divide the data (based on the dot duration)
t_q = quantile(dotdur,binno-1);

%% creating N and K matrix (for each signed coherence, number of total trials (N) and T+ choice trials (K) 
N = zeros(length(coh_set),binno);
K = zeros(length(coh_set),binno);

tax = zeros(binno,1);

for i = 1:length(trials)
    
    coh_id = find(coh_set == coh_info(i));
    
    %% for quantile-based scheme
    dd_index = find(t_q > dotdur(i),1);
    if isempty(dd_index)
        dd_index = binno;
    end
    
    tax(dd_index) = tax(dd_index)+dotdur(i);
    
    N(coh_id,dd_index) = N(coh_id,dd_index)+1;
    
    if Rresponse(i) % T+ choice trials
        K(coh_id,dd_index) = K(coh_id,dd_index)+1;
    end

end

for i = 1:binno
    
    tax(i) = tax(i)./sum(N(:,i));       % mean dot duration for each bin
    
end

%% fit data to diffusion
coh = coh_set./1000;
coh = coh';

bigCoh = repmat(coh,1,binno)
bigT = repmat(tax',length(coh),1)

data = [bigCoh(:) bigT(:) K(:) N(:)];
f = @(x) choiceVdurFromDiffusion_NchooseK_withBias(x,data,ceil(max(dotdur)),0);
thetaGuess = [.2 15 0.01];
errGuess = f(thetaGuess)
thetaFit = fminsearch(f,thetaGuess)
errFit = f(thetaFit)

%% make graph 
figure(10),clf, hold on,

% time-dependent accuracy plot (Figure S1A)
subplot(1,2,1)
tinyOffset = 3;
colors = {[0.7 0.7 0.7],[1 0.9 0.5],[1 0.8 0.4],[1 0.6 0],[1 0 0]};

plot_coh = nonzeros(coh);  % to calculate correct trials (nonzero coherence only)

dt = .05;
tvect = (0:dt:ceil(max(dotdur)))';
B = thetaFit(2);
k = thetaFit(1);
coh_bias = thetaFit(3);
Bup = B*ones(size(tvect));
Blo = -Bup;

cvectAll = plot_coh;
uvect = k * (cvectAll+coh_bias);
y0 = 0;
[pUpAbs, rtUp, rtLo, upDist, loDist,pLoAbs,dvNotAbs,pGT0,xtDist] = ...
    runFPonSeriesPar(uvect,tvect,Bup,Blo,y0); 
        
p = cell(length(plot_coh),1);

pos_coh_set = pos_coh_set/1000;

figure(10)
hold on
for i = 1:length(plot_coh)
    I = find(cvectAll==plot_coh(i));
    
    if coh(i)>=0
    p{i} = pGT0{I}(1:end-1)+cumsum(upDist{I}(1:end));   % at each time point, probability of up-bound crossing (T+ choice)
    else
        p{i} = 1-(pGT0{I}(1:end-1)+cumsum(upDist{I}(1:end)));  % at each time point, probability of up-bound NOT crossing (T- choice)
    end
    
end

set(gca,'xlim',[0 ceil(max(dotdur))+10]);
set(gca,'ylim',[.44 1]);

plot_corr_p = zeros(length(p{1}),length(pos_coh_set));

% Because of the potential bias of the subject, the predicted accuracy for T+ and T- choice might be different even when the motion is of the same strength.
% The following calculates the predicted accuracy for each "unsigned" coherence, 
%   by averaging the (diffusion-fit) accuracy pf the positive and negative motion of the same strength
for i = 1:length(plot_coh)
    I = find(pos_coh_set==abs(plot_coh(i)));
   
    plot_corr_p(:,I) = plot_corr_p(:,I)+p{i};       
end
plot_corr_p = plot_corr_p/2;    

% plotting the predicted accuracy along with the running mean from the data
for i = 1:length(pos_coh_set)
    if rem(i,2)==0
        t = tax + tinyOffset;
    else
        t = tax - tinyOffset;
    end
    
    hFit(i) = plot(tvect,plot_corr_p(:,i),'linewidth',2,'color',colors{i});
    hp(i) = plot(dotdur_matrix(1:end-1),running_mean(i,:),'color',colors{i},'LineWidth',1);
end

set(hFit,'linestyle','-','linewidth',3)
xlabel('Stimulus duration (ms)')
ylabel('Proportion correct')

% add text
for i=1:length(pos_coh_set)
    txtX = ceil(max(dotdur))+10;
    txtY = plot_corr_p(end,i);

    htxt(i) = text(txtX,txtY,sprintf('%.1f%% coh',100*pos_coh_set(i)),'color',colors{i});
    
end

try
set(htxt(5),'verticalAlignment','baseline')
set(htxt(4),'verticalAlignment','top')

set(htxt,'fontsize',14)
catch
end

set(gca, 'TickDir','out');


%% median decision time
cvectAll_length = length(cvectAll)

Qmed = 0;
pGuess = 0; % assuming no lapse

for ic = 1:length(cvectAll)
    
    [~,I]= min(abs(cumsum(upDist{ic}+loDist{ic})-(.5-pGuess))); % finding the time by which either the likelihood of crossing the bound reaches 50% (median decision time)
    
    tmedCoh(ic) = tvect(I); % median decision time for each coherence
    Qmed = Qmed+tvect(I);   
end
tmedCoh
tmed = Qmed/length(cvectAll)    % median decision time across all motion strengths 

    
for i = 1:round(length(cvectAll)/2)
    tmedCoh_plot(i) = tmedCoh(i);
    Coh_plot(i) = abs(cvectAll(i).*100); % coh in unit of %
end

figure(10)
subplot(1,2,2)

plot(Coh_plot, tmedCoh_plot, 'o-')
xlabel('motion strength (%coh)')
ylabel('median decision time (ms)')
box off
grid on

suptitle([titlestring,';  k = ',num2str(k),'; B = ',num2str(B),'; Bias = ',num2str(coh_bias),';  tmed = ',num2str(tmed), ' ms'])

end

