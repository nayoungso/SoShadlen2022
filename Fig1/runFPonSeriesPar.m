function [pUpAbs, rtUp, rtLo, upDist, loDist,pLoAbs,dvNotAbs,pGT0,xtDist] = runFPonSeriesPar(uvect,t,Bup,Blo,y0,sigma)
% Parallel version of runFPonSeries which runs the Fokker-Planck to bounds through a series of stimlus difficulties.
% Just a wrapper to run through a coherence series to generate
% probability surface in time under bounary conditions.
%
% Input args
% ~~~~~~~~~~
% uvect vector of the mean drift rates
% t     time axis. Typically [0:dt:tmax] where dt is 0.05 and tmax is
%       whatever you choose
% Bup and Blo are the upper and lower bounds (same dimension as t)
% y0    starting value of the accumulation; zero by default.
%
% dependencies: FP4.mex (from Roozbeh Kiani) or a slightly modified version
% FP4_BoundCrossZero.mex. I think either will work now because roozbeh made
% some adjustments fto FP4.
%
% return args
% ~~~~~~~~~~~
% pUpAbs  for each x member of uvect, prob of absorption at upper bound
% rtUp    for each x member of uvect, the mean time of absorption at upper bound
% rtLo    for each x member of uvect, the mean time of absorption at lower bound
% upDist  for each x member of uvect, the probability distribution of bound
%         crossings at upper bound (cell array of vectors size of t)
% loDist  for each x member of uvect, the probability distribution of bound
%         crossings at upper bound (cell array of vectors size of t)
% pLoAbs  see pUpAbs (see below)
% dvNotAbs  for each x member of uvect, the distribution of the
%           non-absorbed decision variable at t=tmax. This is a column
%           vector that corresponds to the values of xmesh, which is the
%           grid formed within the FP routine. It's not available, but I
%           haven't needed it.
% pGT0     for each x member of uvect, a vector that is the total
%           probability above 0 at each time step. I think this includes
%           absorbed probability. It's useful for calculations in non-RT
%           experiments when a decision is based on the sign of the dv
%           (cell array of vectors size of t)
% xtDist   for each x member of uvect, a matrix containing the
%          probabability of the dv at each time. I can't remember what's
%          stored at the bounds. This is a huge matrix because its
%          dimensions are size of t by the mesh that separates the upper
%          and lower bounds. (cell array of matrices)


% I added the pLoAbs return argument later mainly to help cross check that
% the sum of the total choice probability is 1. 
% 1 - rtUp(i) + rtLo(i) should be equal to the total probabability of
% nonAbsorbed decision variable at tmax, which is sum(dvNotAbs{i}) 

% 10/11/11 mns spawned from runFPonSeries

m = length(uvect);

% Initialization of the probability surface
if nargin<5
    y0 = 0;
end
if nargin<6
    sigma = 1;
end

% Roozbeh wants to see bounds passed as a single scalar from which the mesh
% is constructed. The real upper bound is the 2nd column of b_change added
% to B. The real lower bound is the 1st column of b_change added to -B.
% b_change = [Blo(:) - minBlo, Bup(:)-maxBup]; 
% figure(10), plot(t,b_change(:,1),t,b_change(:,2))


%allocate vectors and cells
pUpAbs = nan(m,1);
pLoAbs = nan(m,1);
rtUp = nan(m,1);
rtLo = nan(m,1);
upDist = cell(m,1);
loDist = cell(m,1);
dvNotAbs = cell(m,1);
if nargout>7
    pGT0 = cell(m,1);
end
if nargout > 8
    xtDist = cell(m,1);
end

% I pulled these out of the parfor loop because of the warnings
DT_  = t(2)-t(1); % so I can reset in loop
MAXBUP_ = max(Bup(:));
MINBLO_ = min(Blo(:));
B_CHANGE_ = [Blo(:) - MINBLO_, Bup(:) - MAXBUP_];

% Best to create the pool
% isOpen = matlabpool('size') > 0;
% if ~isOpen
%     error('open the matlabpool before calling runFPonSeriesPar');
%     %     matlabpool;
% end
% gcp('nocreate')
parfor i = 1:m
    % All arguments to FP4 are corrupted by FP4. Reset them w/in the loop,
    % before each call.
    
    % time
    dt  = DT_;
    %     max_dur = max(t);
    
    % ygrid, what Roozbeh calls the mesh
    maxBup = MAXBUP_;
    minBlo = MINBLO_;
    maxB = max(abs([maxBup minBlo]));
    % B = maxB; % temporary to use RK's nomenclature
    dx = min(0.1,maxB/100);
    b_margin = [4*sigma; 4*sigma];
    xmesh = (minBlo-b_margin(1)+dx : dx : maxBup+b_margin(2)-dx)'; % RK's method
    % adjust the mesh so that it contains zero.
    if ~any(xmesh==0)
        % adjust dx so that the mesh has zero and is therefore
        % symmetric
        [~, I] = min(abs(xmesh));
        I = min(I(:)); % in case there's a tie
        % xmesh(I) is closest to 0.
        b_margin = b_margin+xmesh(I);
        xmesh = (minBlo-b_margin(1)+dx : dx : maxBup+b_margin(2)-dx)'; % recompute the mesh
        xmesh(abs(xmesh)<1e-12) = 0; % now and then you don't get exactly zero
        % check to be sure that 0 is in the mesh
        if sum(xmesh==0) ~= 1
            % save crashfile.mat
            warning('My logic is wrong in the xmesh adjustment');
        end
    end
    
    b_change = B_CHANGE_;  % Roozbeh wants this argument
    
    % Start position
    uinit = zeros(size(xmesh));
    [~, I] = min(abs(xmesh-y0));
    uinit(I) = 1;   % this should be a 1 where the mesh is closest to y0
    if sum(uinit)~=1
        error('defective initialization or mesh')
    end
    
    
    % [ufinal_,Pt_,Ptb_] = FP4(xmesh,uinit,uvect(i),sigma,b_change,b_margin,dt);
    %     if nargout<8
    %         [ufinal,~,Ptb] = FP4_BoundCrossZero(xmesh,uinit,uvect(i),sigma,b_change,b_margin,dt);
    %     elseif nargout<9
    %         [ufinal,~,Ptb,Pg0] = FP4_BoundCrossZero(xmesh,uinit,uvect(i),sigma,b_change,b_margin,dt);
    %     else
    %         [ufinal,~,Ptb,Pg0,Pxt] = FP4_BoundCrossZero(xmesh,uinit,uvect(i),sigma,b_change,b_margin,dt);
    %     end
    [ufinal,~,Ptb,Pg0,Pxt] = FP4_BoundCrossZero(xmesh,uinit,uvect(i),sigma,b_change,b_margin,dt);
    plo = Ptb(:,1);         %probability of crossing the lower bound at each moment
    pup = Ptb(:,2);        %probability of crossing the upper bound
    % survival_ = Pt_;            %survivor function
    % Pt_ = -diff(Pt_);           %the total probability of crossing either of the bounds at each moment
    pUpAbs(i) = sum(pup);
    pLoAbs(i) = sum(plo);
    % These means only make sense if there is no mass in ufinal. It's the
    % mean rt only for absorbed.
    rtUp(i) = sum(pup(2:end) .* t) ./ pUpAbs(i);
    rtLo(i) = sum(plo(2:end).*t(:)) / pLoAbs(i);
    upDist{i} = pup(2:end); % in RK's code, the 1st element of Ptb is absorption at t=0, that is 1 time step before our 1st
    loDist{i} = plo(2:end);
    dvNotAbs{i} = ufinal;
    pGT0{i} = Pg0;
    xtDist{i} = Pxt;
    %     if exist('Pg0','var')
    %         pGT0{i} = Pg0;
    %     end
    %     if exist('Pxt','var')
    %         xtDist{i} = Pxt;
    %     end
        
end

