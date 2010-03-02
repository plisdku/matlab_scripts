function colormap = orangecrush(varargin)

if nargin > 0
    exponent = varargin{1};
else
    exponent = 1;
end

colormap = zeros(256, 3);

topHalf = 129:256;
bottomHalf = 1:128;

%{
color1 = [0.4118, 0.1922, 1];
color32 = [0.012869, 0.006, 0.0312];
color33 = [0.0313, 0.01, 0];
color40 = [
%}

% 1.  Set the top colors to orange

colormap(topHalf,1) = (0:127)/127;
%colormap(topHalf,2) = 0.2720 * ( ((0:31)/31) - 0.4*((0:31)/31).^3 )  ;
colormap(topHalf,2) = 0.3 * ( ((0:127)/127) ).^2 ;
%colormap(topHalf, 3) = 0.3*((0:31)/31).^3;


% 2.  Set the bottom colors to blue

colormap(bottomHalf,1) = 0.4118 * ((127:-1:0)/127);
colormap(bottomHalf,2) = 0.1922 * ((127:-1:0)/127).^2;
colormap(bottomHalf,3) = (127:-1:0)/127;

colormap = colormap .^ exponent;
