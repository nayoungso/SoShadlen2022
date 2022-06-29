
%% test 1, comparison of FP4 with monte-carlo 

clear all;

k = 0.5;
B = 20;
C = 0.064;
max_t = 500;

dx = 0.1;
dt = 0.05;
test_bound_change = 1;

    %decide about the change of bound height 
t = dt:dt:max_t;    
if test_bound_change,
    bound_change = [zeros(length(t)/2,1); 5-9*10^-5*(t(1:length(t)/2).^2)'];
    bound_change = [-bound_change bound_change];
else
    bound_change = zeros(length(t), 2);
end;
bound_height_highres = repmat([-B B],[length(t) 1])+bound_change;
bound_height_lowres = repmat([-B B],[max_t 1])+bound_change(1/dt:1/dt:end,:);

    %see the bound height 
figure;
plot(dt:dt:max_t, bound_height_highres);
xlabel('Time (ms)');
ylabel('Bound height');


    %run standard FP4 to find the probability of bound crossing
mu = k*C;
sigma = 1;
Bmax = max(B+bound_change(:,2));
b_margin = [4*sigma+(Bmax-B), 4*sigma+(Bmax-B)];
xmesh = (-Bmax-4*sigma+dx : dx : Bmax+4*sigma-dx)';
uinit = zeros(size(xmesh));
uinit(abs(xmesh)==min(abs(xmesh))) = 1;
b_change = bound_change+1-1;    %this stupid -1+1 is necessary, otherwise Matlab will 
                                %just pass the memory to b_change and FP4 will destroy original bound_change 
tic;
[ufinal, Pt_, Ptb_, Pg0_] = FP4(xmesh, uinit, mu, sigma, b_change, b_margin, dt);        
toc;

pcor_FP4 = sum(Ptb_(1:end-1,2)) + Pg0_(end);

        %change time resolution of Pt, Ptb and Pg0 to 1 ms
len = length(Pt_)-1;        
Pt_ = dt*sum(reshape(Pt_(2:end), [1/dt len*dt]))';
Ptb_ = [sum(reshape(Ptb_(2:end,1), [1/dt len*dt]))
        sum(reshape(Ptb_(2:end,2), [1/dt len*dt]))]';
Pg0_ = dt*sum(reshape(Pg0_(2:end), [1/dt len*dt]))';



    %run monte carlo using FP4_MC
mu = k*C;
sigma = 1;
Bmax = max(B+bound_change(:,2));
b_margin = [4*sigma+(Bmax-B), 4*sigma+(Bmax-B)];
xmesh = (-Bmax-4*sigma+dx : dx : Bmax+4*sigma-dx)';
uinit = zeros(size(xmesh));
uinit(abs(xmesh)==min(abs(xmesh))) = 1;
dt_MC = 1;
b_change = bound_change(dt_MC/dt:dt_MC/dt:end,:);    
num_iter = 5000;
bin_width = 20;
tic;
[ufinal_MC, Pt_MC, Ptb_MC, Pg0_MC] = FP4_MC(xmesh, uinit, mu, sigma, b_change, b_margin, dt_MC, num_iter, bin_width);
toc;



    %look at the correct and wrong rt distribution in the two methods 
figure; 
set(gcf, 'Position', [10 300 600 600]);

subplot(3,2,1);
hold on;
bar(bin_width:bin_width:max_t, Pt_MC);
plot(1:max_t, Pt_, 'r');
xlabel('Time (ms)');
ylabel('Survival');

subplot(3,2,2);
hold on;
bar(xmesh, ufinal_MC{1});
plot(xmesh, ufinal, 'r');
xlabel('X');
ylabel('ufinal');

subplot(3,2,3);
hold on;
bar(1:max_t, Pg0_MC);
plot(1:max_t, Pg0_, 'r');
xlabel('Time (ms)');
ylabel('Pg0');

subplot(3,2,5);
hold on;
bar(bin_width:bin_width:max_t, Ptb_MC(:,1));
plot(1:max_t, Ptb_(:,1)*bin_width, 'r');
xlabel('Time (ms)');
ylabel('Pr(wrgRT)');

subplot(3,2,6);
hold on;
bar(bin_width:bin_width:max_t, Ptb_MC(:,2));
plot(1:max_t, Ptb_(:,2)*bin_width, 'r');
xlabel('Time (ms)');
ylabel('Pr(corRT)');


