function y = Fit_ReLU_gain(x,b, varargin)
% 
% y = zeros(size(x));
persistent a
if nargin>2
    a=varargin{1};
    return;
end

y=  max(x-a,0).*b;


end