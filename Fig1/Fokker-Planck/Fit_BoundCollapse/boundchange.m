%  
% function [lb_change,ub_change] = boundchange ( ftype , param , t , dx , lb , ub , lb_change_range , ub_change_range , ignore_dx )
%
% returns the change of lower and upper boundary heights at times defined in t. an increase in the 
% value of the boundary is shown by positive values. note that because lower boundaries are defined 
% by negative values, a positive changes for the lower boundary means getting closer to zero. 
% similarly a negative change for the upper bound means getting closer to zero. 
% 
% Input arguments:
%       ftype -> the function that defines how the bounds change across
%               time. can be 'linear', 'exponential', 'sigmoid', or 'nochange'
%       param -> the parameters for ftype
%              -linear boundary height change: the bounds start to change
%               at time param(2) with the slope param(1)
%              -exponential boundary height change: the bounds start to
%              change at time param(2) with the exponential slope param(1)
%              -sigmoid: the bounds change as a sigmoid function with
%              center param(2) in time and slope param(1)
%       t -> time
%       dx -> resolution of the grid
%       lb -> starting value for the lower bound
%       ub -> starting value for the upper bound
%       lb_change_range -> the range of acceptable changes for the lower bound
%       ub_change_range -> the range of acceptable changes for the upper bound
%       ignore_dx -> if 0 (default value) the changes in bound height will
%               be given in the units of dx, in other words a bound change
%               of dx will be represented as 1 rather than dx
%

% 10/2006   RK

function [lb_change,ub_change] = boundchange ( ftype , param , t , dx , lb , ub , lb_change_range , ub_change_range , ignore_dx )

if nargin<9 || isempty(ignore_dx),
    ignore_dx = 0;
end;

       
switch ftype,
    case 'linear',
        t = t-param(2);
        t(t<0) = 0;
        lb_change = param(1)*t;
        ub_change = -param(1)*t;
    case 'exponential',
        t = t-param(2);
        t(t<0) = 0;
        lb_change = lb*(exp(-param(1)*t)-1);
        ub_change = ub*(exp(-param(1)*t)-1);
    case 'sigmoid',
        t = t-param(2);
        lb_change = lb*(-1./(1+exp(-param(1)*t)));
        ub_change = ub*(-1./(1+exp(-param(1)*t)));
    case 'nochange',
        lb_change = zeros(size(t));
        ub_change = zeros(size(t));
    otherwise,
        error('the requested shape of bound change is not supported!');
end;

if ignore_dx==0,
    lb_change = round(lb_change/dx);
    ub_change = round(ub_change/dx);
end;

    %check that boundary changes do not exceed the expected range for boundary height change
    %boundary should not be changed so that the upper boundary becomes negative or the lower 
    %boundary become positive. They should also not be changed so that 
ub_change(ub_change>ub_change_range(2)) = ub_change_range(2);
ub_change(ub_change<ub_change_range(1)) = ub_change_range(1);
    
lb_change(lb_change<lb_change_range(1)) = lb_change_range(1);
lb_change(lb_change>lb_change_range(2)) = lb_change_range(2);




