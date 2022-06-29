

%function [modelParam, modelLL, predPC, predRT] = fit_FP4 ( trialData , guess , fixed , fitwhat , bc , feedback )
%
% Fits a diffusion model with collapsing bounds to the 2-choice data
%
% Input parameters:
%
% trialData is an array of structs with the following fields
%     Stim    the stimulus presented in the trial, it should contain one of the numbers in the range [1 length(B)]
%     Coh     the coherence of the presented stimulus in the trial, its values should be in the range [0 1]
%     Resp    the response issued by the subject in that trial, again it should be a number in the range [1 length(B)]
%     RT      the reation time of the subject in the trial
% guess     initial guess for the parameters, has 3 elements [k, B, T0], k scales the drift rate, B is the 
%           boundary height and T0 is the non-decision time
% fixed     a vector with three elements, 1 means the corresponding element in guess should remain fixed during
%           fitting, zero means that element is allowed to change. 
%           It provides useful fitting constraints without resorting to constrained fitting methods
% bc        the function and parameters that define the bound change
% fitwhat   a cell with two elements, the first one defines what to be fit ('PC', 'RT', or 'BOTH'), 
%           the second one defines the normalization procedure (use 'TRIALBASED')
% feedback  0 for no feedback, 1 (default) for text feedback at the end of each iteration, 2 for graphical feedback
%
%
% Output parameters:
% 
% modelParam   contains the fitted paramaters
% modelLL      log-likelihood of the model given the fitted parameters
% predPC       predicted probability correct for each coherence level
% predRT       predicted mean RT for each coherence level
%
% _______________________________________________
%
% Example:
% guess = [0.21  20  300];
% fixed = [0     0     0];
% fitwhat = {'BOTH' 'TRIALBASED'};
% bc = struct('func',@boundchange,'ftype','sigmoid','init',[0.001 610],'fixed',[0 0]);
% [modelParam,modelLL,predPC,predRT] = fit_FP4 ( trialData , guess , fixed , fitwhat , bc , 2 );
% _______________________________________________
% 

% 10/2006  RK

function [modelParam, modelLL, predPC, predRT] = fit_FP4 ( trialData , guess , fixed , fitwhat , bc , feedback )

global predPC predRT;
global callNum_ fig_

modelParam = struct ( 'init' , [] , 'fixed' , [] , 'final' , [] , 'se' , [] , 'bc' , [] , 'fitwhat' , [] );
modelLL = [];
predPC = [];
predRT = [];


if nargin < 1 || isempty(trialData),
    error ( 'at least trialData should be provided to the function' );
end;

if nargin < 2 || isempty(guess),
    guess = [0.3 20 300];
end;
if length(guess)~=3,
    error ( 'initial guess does not have correct number of parameters' );
end;

if nargin < 3 || isempty(fixed),      
    fixed = [0 0 0];
elseif length(fixed)~=3 || any(ismember(fixed,[0 1])==0),
    warning ( '''fixed'' is wrongly specified, fixed is automatically set to zeros' );
    fixed = [0 0 0];
end;

if nargin < 4 || isempty(fitwhat),
    fitwhat = {'BOTH' 'TRIALBASED'};
end;

if nargin < 5,
    bc = struct ( 'func' , @boundchange , 'ftype' , 'nochange' , 'init' , [] , 'fixed' , [] , 'final' , [] , 'se' , [] );
else
    if ~isfield(bc,'func') || ~isfield(bc,'ftype'),
        error ( 'bc should be a structure with at least two fields: ''func'' and ''ftype'' ' );
    end;
    if ~isfield(bc,'fixed'),
        bc.fixed = zeros(size(bc.init));
    elseif isempty(bc.fixed),
        bc.fixed = zeros(size(bc.init));
    elseif length(bc.fixed)~=length(bc.init) || any(ismember(bc.fixed,[0 1])==0),
        warning ( '''bc.fixed'' should have values of zero or one. All set to zero' );
        bc.fixed = zeros(size(bc.init));
    end;
    bc.final = [];
    bc.se = [];        
end;

if nargin < 6,
    feedback = 1;
end;

if feedback > 0,
    callNum_ = 1;
    if feedback > 1,
        fig_ = figure;
        hold on;   
        xlabel('Call number' );
        ylabel('-LL');
    end;
end;




        %summerize the data by calculating the percent correct, mean
        %reaction time and standard error of reaction time for each
        %coherence
        %it will considerably speed up the fit because we won't have to
        %calculate these things for each iteration of the loop
summaryData = struct('stimnum',2,'coh_set',unique([trialData.Coh]),'trialnum',[],'PC',[],'RT',[],'RTstd',[],'RTse',[]);
for c = 1 : length(summaryData.coh_set),
    J = find([trialData.Coh]==summaryData.coh_set(c));
    summaryData.trialnum(c)= length(J);
    if summaryData.coh_set(c) ~= 0,
        K = J([trialData(J).Stim]==[trialData(J).Resp]);
        summaryData.PC(c) = length(K)/summaryData.trialnum(c);
    else
        K = J;
        summaryData.PC(c) = 1/summaryData.stimnum;            
    end;
    summaryData.RT(c) = mean([trialData(K).RT]);
    summaryData.RTstd(c) = std([trialData(K).RT]);
    summaryData.RTse(c) = summaryData.RTstd(c)./sqrt(length(K));
end;

gl_guess = [guess bc.init];
gl_fixed = [fixed bc.fixed];

    %save initial settings
modelParam.init = guess;
modelParam.fixed = fixed;
modelParam.bc = bc;
modelParam.fitwhat = fitwhat;

    %check if only model evaluation is requested
if all(fixed==1),
    warning ( 'all parameters are fixed parameters, will return without optimization' ); 
    modelLL = -fitModel_MLEerr([],summaryData,gl_guess,gl_fixed,fitwhat,bc,feedback);
    return;
end;

    %do fitting, apparently there is no need for constrained fitting, our parameter freezing method works fine 
options = optimset ( 'Display' , 'final' , 'TolX' , 1e-2 , 'TolFun' , 1e-2 );
[p,fval] = fminsearch ( @(x) fitModel_MLEerr(x,summaryData,gl_guess,gl_fixed,fitwhat,bc,feedback) , ...
                        gl_guess(gl_fixed==0) , options );
param = getParam(p,gl_guess,gl_fixed);
modelParam.final = param(1:length(guess));
modelParam.bc.final = param(length(guess)+1:end);
modelLL = -fval;

    %to get correct values for predPC and predRT
fitModel_MLEerr([],summaryData,param,ones(size(param)),fitwhat,bc,feedback);



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %returns the -1*log-likelihood of the parameters
            %of the model
            %   if only a subset of parameters are adjustable
            %   the fixed parameters are taken from guess
    
function err = fitModel_MLEerr ( param , summaryData , guess , fixed , fitwhat , bc , feedback )

global predPC predRT;
global callNum_ fig_


param = getParam(param,guess,fixed);

k = param(1);
B = param(2);
T0 = param(3);


stimNum = 2;
M = k*cos(2*pi*(0:stimNum-1)/stimNum)';
V = 0.5*ones(size(M));
Model = [+1 -1];
bcparam = param(4:end);
            
dx = min(0.1,B/100);
dt = 0.1;
            
for c = 1 : length(summaryData.coh_set),
    Pt = cell(1,size(Model,1));
    Ptlow = cell(1,size(Model,1));
    Pthigh = cell(1,size(Model,1));
    for i = 1 : size(Model,1),
            %calc drif rate and standard deviation of the diffusion
        k_ = Model(i,:) * M * summaryData.coh_set(c);
        sigma_ = sqrt( (Model(i,:).^2) * V );
            %use fokker-planck equation to find the probability of termination of the diffusion at each moment
            %save the termination probability for the next stage
        [Pt{i},Ptlow{i},Pthigh{i}] = termpdf(k_,sigma_,B,dx,dt,1e3,bc,bcparam); 
            
            %if coh is 0 and diffusion variances are equal then all diffusion processes will be equivalent
            %and there is no reason to calculate the same probability densities many times
        if i==1 && summaryData.coh_set(c)==0 && all(diff(V)==0),
            for j = 2 : size(Model,1),
                Pt{j} = Pt{i};
                Ptlow{j} = Ptlow{i};
                Pthigh{j} = Pthigh{i};
            end;
            break;
        end;
    end;
            
        %now find the probability that the target diffusion finish but others don't
    if size(Model,1) == 1,
        Pt1 = Pt{1};
    else
        maxsize = 0;    
        for i = 1 : size(Model,1),
            maxsize = max(maxsize,length(Pt{i}));
        end;
        P = zeros(maxsize,size(Model,1));
        P(1:length(Pt{1}),1) = Pt{1};                   %we need the termination probability density function for the target process
        for i = 2 : size(Model,1),
            P(1:length(Pt{i}),i) = 1-cumsum(Pt{i});     %we need the survivor function for the other processes
        end;
                
        Pt1 = prod(P,2);    %the probability that the target process be the one that finishes first is the product of the probability
                            %density of termination of the first process and survivor functions of the other processes
    end;
        %find the probability of getting to the lower and upper bounds separately, conditional to the target process finish before others
    if ~isempty(Ptlow{1}),
        Pt1low = zeros(size(Pt1));                    
        Pt1low(1:length(Ptlow{1})) = Ptlow{1};
    end;
    if ~isempty(Pthigh{1}),
        Pt1high = zeros(size(Pt1));
        Pt1high(1:length(Pthigh{1})) = Pthigh{1};
    end;
            
        %since the area under Pt1 is not necessarily equal to one we need to
        %nomalized it by calculating Pt1/sum(Pt1)
    predRT(c) = (dt*(1:length(Pt1))) * Pt1high/sum(Pt1high);
            
        %the correct probability is the probability of the target process
        %terminate before other processes times the probability of hitting
        %the upper boundary in the process
    L = (Pt1high + Pt1low)~=0;
    predPC(c) = (Pt1high(L)./(Pt1high(L)+Pt1low(L)))' * Pt1(L);                 
end; 
    
LL_PC = 0;
LL_RT = 0;
switch fitwhat{2},
    case 'TRIALBASED',   
        for c = 1 : length(summaryData.coh_set),
            if summaryData.coh_set(c) ~= 0,
                p = predPC(c);
                q = 1-predPC(c);
                if q == 0,      q = 1/summaryData.trialnum(c)/2;        end;
                LL_PC = LL_PC + summaryData.trialnum(c)*(summaryData.PC(c)*log(p)+(1-summaryData.PC(c))*log(q));
            end;
            LL_RT = LL_RT + lognormpdf(predRT(c)+T0,summaryData.RT(c),summaryData.RTse(c));
        end;
    case 'NORMALIZED',
        for c = 1 : length(summaryData.coh_set),
            if summaryData.coh_set(c) ~= 0,
                p = predPC(c);
                q = 1-predPC(c);
                if q == 0,      q = 1/summaryData.trialnum(c)/2;        end;
                LL_PC = LL_PC + (summaryData.PC(c)*log(p)+(1-summaryData.PC(c))*log(q));
            end;
            LL_RT = LL_RT + lognormpdf(predRT(c)+T0,summaryData.RT(c),summaryData.RTstd(c));
        end;
end;
        

switch fitwhat{1},
    case 'RT',      
        err = -sum(LL_RT);
    case 'PC',      
        err = -sum(LL_PC);
    case 'BOTH',    
        err = -sum(LL_PC+LL_RT);
end;


    %print progress report!
if feedback > 0,
    fprintf ( '******************   *******************   *******************\n' );
    fprintf ( 'run %d\n' , callNum_ );
    fprintf ( 'param = ' );     fprintf ( '%f, ' , param(1:3) );        fprintf ( '\n' );
    fprintf ( '        ' );     fprintf ( '%f, ' , param(4:end) );      fprintf ( '\n' );
    fprintf ( 'LL_PC: ' );      fprintf ( '%6.2f, ' , LL_PC );             fprintf ( '\n' );
    fprintf ( 'LL_RT: ' );      fprintf ( '%6.2f, ' , LL_RT );             fprintf ( '\n' );
    fprintf ( 'sum: %f\n' , -err );
    fprintf ( 'predPC:  ' );    fprintf ( '%6.2f ,' , predPC );            fprintf ( '\n' );
    fprintf ( '   goal: ' );    fprintf ( '%6.2f ,' , summaryData.PC );    fprintf ( '\n' );
    fprintf ( 'predRT:  ' );    fprintf ( '%6.2f ,' , predRT+T0 );         fprintf ( '\n' );
    fprintf ( '   goal: ' );    fprintf ( '%6.2f ,' , summaryData.RT );    fprintf ( '\n' );
    if feedback > 1,
        figure(fig_);
        plot(callNum_,err,'.');
        drawnow;
    end;
end;
callNum_ = callNum_ + 1;


    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
            %this function retrieves the full parameter set given the adjustable and
            %fixed parameters
function param2 = getParam ( param1 , guess , fixed )

param2(fixed==0) = param1;              %get adjustable parameters from param1
param2(fixed==1) = guess(fixed==1);     %get fixed parameters from guess



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %returns log of normpdf, this function helps avoinding round off errors for very small probabilities.
            %if you use log(normpdf()) instead of lognormpdf you can get 0 for small probabilities and it's 
            %detrimental to log-likelihood fitting

function l = lognormpdf ( x , mu , sigma )

if sigma==0,
    if all(x==mu),  
        l = 0;
    else
        l = -inf;       
    end;
else
    d = sqrt(2*pi);
    l = -(x-mu).^2/(2*sigma.^2)+log(1./(d*sigma));
end;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
function [Pt, Ptlow, Pthigh] = termpdf ( k , sigma , B , dx , dt , tchunk , bc , bcparam )
       
lb_margin = ceil(4*sigma);                  %size of the lower margin, for calculation of the probability of lower bound crossing
ub_margin = ceil(4*sigma);                  %size of the upper margin, for calculation of the probability of upper bound crossing
xinit = -B-lb_margin+dx;                    %height of the upper bound plus the upper margin    %should be positive
xend = B+ub_margin-dx;                      %height of the lower margin plus the lower margin   %should be negative
lb_change_range = [-lb_margin, B-2*dx];     %range of possible changes in lower boundary position
ub_change_range = [-B+2*dx, ub_margin];     %range of possible changes in upper boundary position
xmesh = (xinit:dx:xend)';                   %spatial grid
Px = zeros(size(xmesh));                    
Px( abs(xmesh) == min(abs(xmesh)) ) = 1.0;    %initial probability density on the spatial grid


Pt = [];     Ptlow = [];     Pthigh = [];
    %continue the diffusion process until the cummulative probability of
    %bund-crossing is greater than 0.999
while abs(sum(Pt)-1.0)>1e-3,
    if ~isempty(bc.func),
        t = length(Pt)*dt+(dt:dt:tchunk);
        [lb_change, ub_change] = feval(bc.func,bc.ftype,bcparam,t,dx,-B,+B,lb_change_range,ub_change_range,1);
        lb_change = lb_change(:);
        ub_change = ub_change(:);
        %figure; plot(t,-B+lb_change,t,B+ub_change);
    end;
    
%%% Mac computers handle the memory incorrectly, make sure that the vectors remain valid on Mac
lb_margin = ceil(4*sigma);                  
ub_margin = ceil(4*sigma);                  
xmesh = (xinit:dx:xend)'; 
%%%
    [Px, Pt_, Ptb_] = FP4_BoundCrossZero ( xmesh, Px, k, sigma, [lb_change, ub_change], [lb_margin; ub_margin], dt );
    
    Pt = [Pt; -diff(Pt_)];
    Ptlow = [Ptlow; Ptb_(2:end,1)];
    Pthigh = [Pthigh; Ptb_(2:end,2)];
end;

    %check if boundary expansion by lb_margin and ub_margin was enough to calculate the correct proabability later
if lb_margin>0 || ub_margin>0,
    err = Pt-(Ptlow+Pthigh);
    if any(abs(err)>1e-6),
        warning('probability correct is inaccurate!\n');
    end;
end;



