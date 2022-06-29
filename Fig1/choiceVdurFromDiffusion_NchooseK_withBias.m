function [err, p] = choiceVdurFromDiffusion_NchooseK_withBias(theta,data,TMAX,pLapse)
% returns minus log likelihood of data given theta=[k B] added TMAX term
% for speed. This version assumes that data has 4 columns: coh, dur,
% nCorrect, nTotal. The error assumes binomial (n-choose-k observations)

% History. MNS wrote it ~Aug 2013 to analyze MIP/LIP data in Lafuente et al
% (2015). I made minor changes and updated the parallel matlab calls in Aug
% 2015.

% NS added bias term

%% preliminaries
if nargin < 3
    TMAX = 2000;
end
if ~exist('pLapse','var')
    pLapse = 0;
end
%% unpack theta
k = theta(1);
B = theta(2);
coh_bias = theta(3);	% bias term

%% unpack data
scoh = data(:,1); 
dur = data(:,2); % in units of ms
nCor = data(:,3);
N = data(:,4);

nPoints = size(data,1); 
cvect = unique(scoh); 


%% run FP
dt = 0.05;
t = (0:dt:TMAX)';
Bup = B*ones(size(t));
Blo = -Bup;
y0 = 0;
uvect = k*(cvect+coh_bias);

[pUpAbs, rtUp, rtLo, upDist, loDist,pLoAbs,dvNotAbs,pGT0,xtDist] = ...
    runFPonSeriesPar(uvect,t,Bup,Blo,y0);

%% likelihoods
p = nan(nPoints,1);
for i = 1:nPoints
    iCoh = find(cvect == scoh(i));
    [d,iT] = min(abs(t-dur(i)));
    if d>dt
        error('this should not happen')
    end
    if isempty(iT)
        sprintf('empty iT')
    end
    if pLapse == 0
        p(i) = binopdf(nCor(i),N(i),pGT0{iCoh}(iT) + sum(upDist{iCoh}(1:iT)));
    else
        p1 = pLapse*0.5 + (1-pLapse)*(pGT0{iCoh}(iT) + sum(upDist{iCoh}(1:iT)));
        p(i) = binopdf(nCor(i),N(i),p1);
    end
    
end

err = -sum(log(p));
fprintf('k=%.3f\tB=%.3f\tC_{bias}=%.3f\terr=%g\n',k,B,coh_bias,err);

