function VD_CorCE_VarCE_FigS2


for celltype = 1:2
    
    if celltype==1  
        groupname = 'Leader';
    elseif celltype==2
        groupname = 'Supporter';
    end
    
    path = ['~/Variable duration task/' groupname '/']
    
    
    clear trials_temp DSP011_temp
    clear res_DOTS_temp res_SACC_temp mean_DOTS_temp mean_SACC_temp
    clear res_DOTS res_SACC mean_DOTS mean_SACC
    
    
    filename_a = dir(strcat(path,'A/r/','*_DSP*.mat'));
    filename_b = dir(strcat(path,'A/g/','*_DSP*.mat'));
    filename_c = dir(strcat(path,'B/r/','*_DSP*.mat'));
    filename_d = dir(strcat(path,'B/g/','*_DSP*.mat'));
    
    cutoff = 250;    % activity after 250ms from the motion offset are excluded
    
    binwidth = 50;  
    binno = 5;
    phi_init = 0.5;
    crit_fr = 30;  
    
    pre_dots = 160; %based on IEM_new_putativeDec_ranksum.m, p<0.01
    post_dots = pre_dots + binwidth*(binno-1);
    
    pre_dots_varce = pre_dots; post_dots_varce = post_dots;
    
    idstart = 0;
    
    coh_set = [-80 -40 0 40 80];    % low coherence only to make sure that the diffusion process hasn't met the bound
    
    res_DOTS = nan;
    
    
    for p = 1:4
        
        clear filename
        
        if p==1
            taskid = 20; filename = filename_a; filepath = [path,'/A/r/']; % neurons preferring right/up target in config A
        elseif p==2
            taskid = 20; filename = filename_b; filepath = [path,'/A/g/']; % neurons preferring left/down target in config A
        elseif p==3
            taskid = 21; filename = filename_c; filepath = [path,'/B/r/']; % neurons preferring right/up target in config B
        elseif p==4
            taskid = 21; filename = filename_d; filepath = [path,'/B/g/']; % neurons preferring left/down target in config B
        end
        
        
        for i = 1:length(filename)
            
            clear trials trials_phi phi_times dot_times1 dot_off_times slib_times sacc1_times
            clear t_array_phi t_array_d t_array_de t_array_sl t_array_s
            clear trial_spike_ts DSP011 ds_Gaddress_coh
            clear phi_spk
            load([filepath,filename(i).name]);
            
            dspsize = size(DSP011,2);
            
            trialno = 0;
            Gno_coh = zeros(1,length(coh_set));
            
            for j = 1:length(trials)
                
                if (trials(j).taskid == taskid)&(trials(j).response>=0)&(trials(j).dot_coh <=max(coh_set))
                    if mod(p,2)==0
                        trials(j).t1_dir = trials(j).t1_dir+180;  % flip the sign for left/down preferring neurons (p=2 or 4) to be combined with right/up preferring neurons
                    end
                    
                    tempt = trials(j).time_dots_off(1)-trials(j).time_dots_on(1);
                    if (tempt+cutoff >= pre_dots+binwidth*binno)&&(trials(j).mean_maxFR >= crit_fr)
                        
                        trialno = trialno+1;
                        trial_spike_ts(trialno,1:dspsize) = DSP011(j,1:dspsize);
                        
                        dir_sign = (trials(j).dot_dir ~= trials(j).t1_dir); % 0 for tin, 1 for tout dir
                        
                        if trials(j).dot_coh >= 0
                            coh(trialno)= trials(j).dot_coh*(1-2*dir_sign);
                        else
                            coh_temp = 1;
                            coh(trialno) = coh_temp*(1-2*dir_sign);
                        end
                        
                        indext = find(coh_set==coh(trialno));
                        
                        Gno_coh(indext) = Gno_coh(indext)+1;
                        ds_Gaddress_coh{indext}(Gno_coh(indext)) = trialno;
                        
                        cho(trialno) = abs(dir_sign-trials(j).response); % 1 for tin, 0 for tout choice
                        
                        
                        dots_time1(trialno) = trials(j).time_dots_on(1);
                        tempt = trials(j).time_dots_off(1)-trials(j).time_dots_on(1);
                        t_array_d(trialno,:) = [pre_dots min(post_dots,tempt+cutoff)];
                        
                    end
                end
                
            end
            
            
            if trialno > 0 %~isempty(trial_spike_ts)
                [res_DOTS_temp, mean_DOTS_temp, trialid_DOTS_temp] = spk_corce2(trial_spike_ts,ds_Gaddress_coh,[dots_time1(:)],[pre_dots post_dots],t_array_d,binwidth,idstart);
                
                idstart = trialid_DOTS_temp(end);
                
                if isnan(res_DOTS) % if i==1
                    res_DOTS = res_DOTS_temp; mean_DOTS = mean_DOTS_temp;  trialid_DOTS = trialid_DOTS_temp;
                else
                    res_DOTS = [res_DOTS;res_DOTS_temp];        mean_DOTS = [mean_DOTS; mean_DOTS_temp];        trialid_DOTS = [trialid_DOTS; trialid_DOTS_temp];
                end
                
                clear res_DOTS_temp mean_DOTS_temp res_DOTS_varce mean_DOTS_varce trialid_DOTS_varce
                
            end
            
        end % end of cell-specific part
        
    end
    
    
    f = @(phi) CorCE_opt(phi,binno,res_DOTS,mean_DOTS); % fitting phi based on the theoretical expectation of CorCE for diffusion process
    
    % printing the phi-optimizing process
    errInit = f(phi_init)
    phiFit = fminsearch(f,phi_init);
    errFit = f(phiFit)
    
    % compute CorCE & VarCE with the fitted phi
    [CorCE, theo_CorCE] = compute_CorCE(phiFit,binno,res_DOTS,mean_DOTS);
    VarCE = nanmean(res_DOTS.^2 - phiFit*mean_DOTS);
    
    % linear regression for VarCE 
    timex_dots = pre_dots:binwidth:post_dots;
    [b1,~,st1] = glmfit(timex_dots,VarCE,'normal','link','identity');
    varce_regy = b1(1)+b1(2)*timex_dots;
    if celltype==1
        disp(['VarCE slope for leaders = ',num2str(b1(2)),' (p = ',num2str(st1.p(2)),')'])
    elseif celltype==2
        disp(['VarCE slope for supporters = ',num2str(b1(2)),' (p = ',num2str(st1.p(2)),')'])
    end
    
    % for plotting the actual and theoretical CorCE (Figure S2B)
    for i = 1:binno-1
        CorCE_init(i) = CorCE(1,i+1);
        CorCE_adj(i) = CorCE(i,i+1);
        
        theo_CorCE_init(i) = theo_CorCE(1,i+1);
        theo_CorCE_adj(i) = theo_CorCE(i,i+1);
    end

    
    % plotting 
    figure(27)
    if celltype == 1
        clf;
        subplot(231)
    elseif celltype == 2
        subplot(232)
    end
    imagesc(CorCE)
    box off
    set(gca,'TickDir','out');
    set(gcf,'Color','w');
    
    
    if celltype == 1
        title(['Leaders, phi = ',num2str(phiFit)]);
        subplot(234)
    elseif celltype == 2
        title(['Supporters, phi = ',num2str(phiFit)]);
        subplot(235)
    end
    hold on
    plot(2:binno,CorCE_init,'ko','MarkerEdgeColor','k','MarkerSize',6);
    plot(2:binno,CorCE_adj,'ko','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6);
    plot(2:binno,theo_CorCE_init,'k--','LineWidth',1);
    plot(2:binno,theo_CorCE_adj,'k-','LineWidth',1);
    xlim([2 binno]);
    ylim([0.4 1]);
    set(gca,'TickDir','out');
    
    subplot(2,3,[3 6])
    timex_dots = pre_dots_varce:binwidth:post_dots_varce;
    hold on
    if celltype == 1
        plot(timex_dots,VarCE,'ro','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',6)
        plot(timex_dots,varce_regy,'r--');
    elseif celltype == 2
        plot(timex_dots,VarCE,'bo','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',6)
        plot(timex_dots,varce_regy,'b--');
    end
    
    set(gca,'TickDir','out');
    set(gcf,'Color','w');
    xlim([timex_dots(1) timex_dots(end)])
    hold off
    
    set(gca,'TickDir','out');
    set(gcf,'Color','w');
    box off
    ylim([0.19 0.58])
    
    clear all
end % end of celltype-loop


end

