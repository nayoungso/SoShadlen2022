


function [ufinal, Pt, Ptb, Pg0, Pxt] = FP4_MC(xmesh, uinit, k, sigma, b_change, b_margin, dt, num_iter, bin_width)

% Matlab call:
% 
% 	[ufinal, Pt, [Ptlow Pthigh], Pg0, Pxt] = FP4_MC ( xmesh, uinit, k, sigma, [lb_change ub_change], [lb_margin; ub_margin], dt, num_iter);
% 
% Description:
% 
% 	This is a Monte-Carlo version of FP4. See FP4 for further description.
%   Below I describe the differences between FP4 and FP4_MC
%   
%   1. FP4_MC has two extra input arguments; num_iter and bin_width
%           num_iter    defines how many times the monte carlo simulation must run. the
%                       higher this number the more accurate the output. however, memory
%                       can a serious restriction.
%           bin_width   defines how coarseness of the final distribution of bound
%                       crossing times. bin_width is in units of dt
%   2. FP4_MC offers a different format of ufinal compared to FP4
%           ufinal      is a cells with four elements
%                       ufinal{1} similar to ufinal returned by FP4, that is the final distribution on the spatial grid 
%                       ufinal{2} bound crossing times of individual simulated trials, NaN if the bounds were not crossed  
%                       ufinal{3} final value of accumulated evidence at the end of simulated trials, clipped to bound after bound crossing 
%                       ufinal{4} final value of accumulated evidence assuming the value is not clipped to bound after crossing it 
%   3. FP4_MC can accept aribitrary initial distributions but FP4 can accept only delta
%   function
%   4. xmesh can be irregular for FP4_MC but it has to be regular (equal spacing of
%   spatial grid) for FP4
%   5. unlike FP4, FP4_MC does not modify the input arguments
%   6. FP4_MC allows the two bounds to cross zero, but FP4 does not. FP4_MC is actually
%   equivalent to FP4_BoundCrossZero 
%   

%
% Roozbeh Kiani,        Sep 25, 2008
%


if nargin<7 || isempty(dt),
    dt = 1;
end;

if nargin<8 || isempty(num_iter),
    num_iter = 5000;
end;

if size(xmesh,2)~=1 || size(uinit,2)~=1,
    error('FP4_MC:InvalidArg', 'the spatial grid and the initial probability distribution should be column vectors');
end;

if length(xmesh) ~= length(uinit),
    error('FP4_MC:InvalidArg', 'mismatch in the size of the spatial grid and the initial probability distribution');
end;

if any(diff(xmesh)<=0),
    error('FP4_MC:InvalidArg', 'the spatial grid should be monotonically increasing\n\t\te.g., it can''t go from the upper bound to the lower bound');
end;

if size(b_change,2)~=2 || size(b_change,1)<3,
    error('FP4_MC:InvalidArg', 'bound_change should be a matrix with two columns, each of\n\t\twhich corresponding to the change of one of the boundaries across time');
end;
    
if numel(b_margin)~=2,
    error('FP4_MC:InvalidArg', 'bound_margin should be a vector with two elements');
end;

    %how long is each iteration in units of dt
max_t = size(b_change, 1);

if max_t*num_iter > 1e7,
    warning('FP4_MC:ArgsTooBig', 'the number of iteration and the length of time are too big\n\t\tMatlab may run out of memory');
end;

    %retrieve the bound height and its change over time
B = [xmesh(1)+b_margin(1), xmesh(end)-b_margin(2)];     %initial bound
bound_height = repmat(B,[max_t 1]) + b_change;


    %randomly sample from uinit to generate the initial distribution  
Zinit = randsample(xmesh, num_iter, 'true', uinit);

    %generate the sequence of randum numbers, that's is the magnitude of momentary
    %evidence according to a normal distribution with mean k and standard deviation
    %sigma
    %change k and sigma according to dt
k = k*dt;
sigma = sigma*sqrt(dt);
Z = normrnd(k, sigma, [num_iter, max_t]);
    %add the initial values
Z(:,1) = Z(:,1) + Zinit;    
    %calculate the cumulative evidence
Z = cumsum(Z, 2);

    %find the choice and bound crossing time for each trial 
choice = nan(num_iter, 1);
bound_crossing = zeros(num_iter, 1);
bound_crossing_t = nan(num_iter, 1);
    %first find trials in which bound crossing took place
[I, J] = find(Z>=repmat(bound_height(:,2)',[num_iter 1]) | Z<=repmat(bound_height(:,1)',[num_iter 1]));
    %now find the first bound crossing time on these trials, ultimately we'd like I to be
    %an index to the trials with bound crossing and J be the first time of bound crossinf
    %on those trials 
[I,ind] = sort(I, 1, 'ascend');
I = flipud(I);
J = flipud(J(ind));
[I,ind] = unique(I);
J = J(ind);
    %now store the bound crossing times and the choice (1 for the upper bound, 0 for the
    %lower bound)  
bound_crossing(I) = 1;
bound_crossing_t(I) = J;
choice(I) = (Z(sub2ind(size(Z),I,J)) >= bound_height(J,2));

ufinal = cell(4,1);
ufinal{1} = hist(Z(bound_crossing==0,end),xmesh) / num_iter;
ufinal{2} = bound_crossing_t;
ufinal{3} = Z(:,end);
ufinal{3}(choice==1) = bound_height(bound_crossing_t(choice==1),2);
ufinal{3}(choice==0) = bound_height(bound_crossing_t(choice==0),1);
ufinal{4} = Z(:,end);

bin_centers = bin_width:bin_width:max_t;

Ptb = nan(length(bin_centers),2);
Ptb(:,1) = hist(bound_crossing_t(choice==0), bin_centers) / num_iter;
Ptb(:,2) = hist(bound_crossing_t(choice==1), bin_centers) / num_iter;

Pt = 1-cumsum(sum(Ptb,2));

if nargout>3,    
    disp('here!');
    for trial = find(bound_crossing)',
        Z(trial, bound_crossing_t(trial):end) = NaN;
    end;
    Pg0 = sum(Z>=0) / num_iter;
    if nargout>4,
        Pxt = hist(Z, xmesh) / num_iter;
    end;
end;







