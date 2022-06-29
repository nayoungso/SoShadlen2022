%
%
% Matlab call:
% 
% 	[ufinal, Pt, Ptlow, Pthigh, Pg0, Pxt] = FP4 ( xmesh, uinit, k, sigma, [lb_change ub_change], [lb_margin; ub_margin], dt);
% 
% Description:
% 
% 	Propagates the probability density <uinit> on the spatial grid <xmesh> across time, using the diffusion
%    parameters k (drift rate) and sigma (diffusion standard deviation). The duration of the propagation is
%    defined by dt and size of [lb_change ub_change] (see below). Don't include Pxt in the output arguments
%    to keep memory management efficient and fast.
% 
% 
% input ->  [xmesh]						//a column vector specifying the spatial grid for probability density propagation in space
% 										//the values in xmesh are monotonically increasing and range from low values (lower bound - lb_margin + dx) to
%                                       //high values (upper bound + ub_margin - dx). xmesh should be regular. define it as:     
%                                       //lower bound - lb_margin + dx : dx : upper bound + ub_margin - dx
%
%           [uinit]						//a column vector specifying the initial probabilities on the grid
%                                       //it is usually a delta function, such as:
%                                       //uinit = zeros(size(xmesh));   uinit(xmesh==0) = 1;
%
%           k							//a scalar specifying the drift rate
%
%           sigma						//a scalar specifying the standard deviation of the diffusion process
%
%           [lb_change ub_change]		//a matrix consisting of two columns, each column shows how the corresponding margin changes across time
%          								//the number of rows in this matrix defines the number of time steps for the propagation of the probability
%                                       //density. stepnum*dt is equal to the total propagation time
%                                       //make sure that the change in the bound height is valid: bound height should change so that the upper 
%                                       //bound stay above zero and the lower bound stay below zero all the time, also neither of them can 
%                                       //exceed the xmesh 
%
%           [lb_margin; ub_margin]		//a column vector specifying the margin of the lower and upper boundaries
%          								//-in Fokker-Planck we normally get the probability of boundary crossing
%            							// but don't distinguish the probability of crossing each of the boundaries. To find the probability
%                         				// of crossing the upper and lower boundaries separately we should reduce xinit and increase xend by 
%                                       // many folds of sigma and calculate how much of the probability falls beyond the actual target boundaries.
%                         				// lb_margin and ub_margin define the distance of the actual boundaries from xinit and xend (the starting
%                                       // and ending point of xmesh. When lb_margin and ub_margin are zero xinit and xend correspond to the boundaries.
%
%           [dt]						//the size of each time step, specifying the resolution of the temporal grid
% 
% 
% output->  [ufinal]					//a column vector, final probability density at the end of the propagation
%
% 			[Pt]						//a column vector, it is the survivor function: the probability of NOT crossing the bounds up to each moment 
%                                       //during the propagation
%
%           [Ptlow, Pthigh]				//a matrix with two columns, the first column corresponds to the probability of crossing the lower bound
%          								//at each moment. the second column corresponds to the probability of crossing the upper bound at each moment
%
%           [Pg0]						//a column vector, the probability of staying above zero at each moment. it is useful when
%          								//the diffusion process terminates before bound crossing
%                                       //I have written the program so that it includes also half of the mass in the spatial grid location that 
%                                       //corresponds to zero 
%
%           [Pxt]						//a matrix with columns corresponding to the probability density at each moment during
%          								//the propagation. its calculation and storage requires a huge amount of memory. don't ask for 
%                                       //it if you are not sure that your computer can handle it or not. 
%
%   

%
% Roozbeh Kiani,        Aug 15, 2006
%
