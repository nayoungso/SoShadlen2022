function Pulse_behav_modelSimul_FigS1



% simulating random one pulse only, with a set sensitivity (that also
% explains each monkey's pmf. using the coh_avg)
% then fit the fake data using the full logit function (P1+P2 / Pmax+Pmin)
% and compare the coefficients with the real one

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coh_set = [0 4 8 16 32 64]';  
dir_set = [-1 1]';


b0 = 0;
b1_array = 0.065:0.001:0.071;    % sensitivity ranges (pooled from two monkeys and two task variants). Note that this sensitivity is calculated using the avg motion strength (eq. 2 in Methods) 
simulno = 160;  % 160 simulations for each of 7 different sensitivities (b1_array) yields 1120 simulations 
trialno = 10000; 

    function cp = calculate_cp(motioncoh)
        cp = (1+exp(-(b0+b1*motioncoh))).^-1;   % to calculate the probability of choosing T+
    end


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
        AllPulses_cp = calculate_cp(p1p2_comb_array);   % probability of choosing T+, if using both pulses (Model-2 in Supplementary Information)
         
        for i = 1:trialno
            
            randomOne_cp(i) = calculate_cp(2*signed_coh_array(i,random_idx(i))); % probability of choosing T+ if using one randomly selected pulse (Model-1 in Supplementary Information)   
            
            AllPulses_choice(i) = random('binomial',1,AllPulses_cp(i));     % generating choices based on Model-2
            randomOne_choice(i) = random('binomial',1,randomOne_cp(i));     % generating choices based on Model-1
            
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
        
        % Model-2 (both pulses)
        [b_p1p2_sum,~,~] = glmfit([signed_coh_array(:,1) signed_coh_array(:,2)],[AllPulses_choice ones(trialno,1)],'binomial','link','logit');
        [b_pmaxpmin_sum,~,~] = glmfit([maxcoh_array mincoh_array],[AllPulses_choice ones(trialno,1)],'binomial','link','logit');
        
        % storing coefficients from each iteration of simulation
        bp1_sum(k,j) = b_p1p2_sum(2); bp2_sum(k,j) = b_p1p2_sum(3);
        bmax_sum(k,j) = b_pmaxpmin_sum(2); bmin_sum(k,j) = b_pmaxpmin_sum(3);
        
    end
end

mean_bp1_rand = nanmean(bp1_rand(:)); std_bp1_rand = nanstd(bp1_rand(:));
mean_bp2_rand = nanmean(bp2_rand(:)); std_bp2_rand = nanstd(bp2_rand(:));

mean_bp1_sum = nanmean(bp1_sum(:)); std_bp1_sum = nanstd(bp1_sum(:));
mean_bp2_sum = nanmean(bp2_sum(:)); std_bp2_sum = nanstd(bp2_sum(:));


mean_bmax_rand = nanmean(bmax_rand(:)); std_bmax_rand = nanstd(bmax_rand(:));
mean_bmin_rand = nanmean(bmin_rand(:)); std_bmin_rand = nanstd(bmin_rand(:));

mean_bmax_sum = nanmean(bmax_sum(:)); std_bmax_sum = nanstd(bmax_sum(:));
mean_bmin_sum = nanmean(bmin_sum(:)); std_bmin_sum = nanstd(bmin_sum(:));


figure(11)
clf;

hold on;


errorbar(mean_bmax_rand,mean_bmin_rand-0.055,2*std_bmin_rand,2*std_bmin_rand,2*std_bmax_rand,2*std_bmax_rand,'o','Color','b',...
    'MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor','b','Capsize',0);   % note that bmin (b1 of Model-1) is adjusted to be shown in the figure with axis bracket

errorbar(mean_bmax_sum,mean_bmin_sum,2*std_bmin_sum,2*std_bmin_sum,2*std_bmax_sum,2*std_bmax_sum,'o','Color','k',...
    'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k','Capsize',0);

% data
plot(0.0318,0.0389,'ro','MarkerFaceColor','r','MarkerSize',10);
plot(0.0327,0.0363,'go','MarkerFaceColor','g','MarkerSize',10);
plot(0.0301,0.0359,'r^','MarkerFaceColor','r','MarkerSize',10);
plot(0.0301,0.0370,'g^','MarkerFaceColor','g','MarkerSize',10);



xlabel('b_{max}','FontSize',16); ylabel('b_{min}','FontSize',16);
set(gca,'TickDir','out')
axis square

xticks([0.02 0.03 0.04 0.055 0.065 0.075 0.085])
xticklabels({'0.02','0.03','0.04','0.11','0.12','0.13','0.14'});    % to accomodate the axis break
yticks([0.02 0.03 0.04 0.055 0.065 0.075 0.085])
yticklabels({'0.02','0.03','0.04','0.11','0.12','0.13','0.14'});

xlim([0.01 0.09]);
ylim([0.01 0.09]);



end
