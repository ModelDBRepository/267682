function y = Fit_ReLU(x,a,b)
% RELU   A line made of two pieces
% that is not continuous.
% Website https://www.mathworks.com/help/curvefit/fit.html#bto2vuv-1-fitType

% y = zeros(size(x));

y=  max(x-a,0).*b;


end