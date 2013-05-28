function [colormap, scale] = orangecrush(varargin)

if nargin > 0
    exponent = varargin{1};
else
    exponent = 1;
end

if nargin > 1
    numColors = varargin{2};
else
    numColors = 256;
end

if rem(numColors, 2) == 0
    numColors = numColors + 1;
end

numHalfColors = floor(numColors/2);

colormap = zeros(numColors, 3);

%topHalf = 129:256;
%bottomHalf = 1:128;

vv = linspace(0, 1, numColors);

topHalf = vv(vv > 0.5);
bottomHalf = vv(vv < 0.5);

%{
color1 = [0.4118, 0.1922, 1];
color32 = [0.012869, 0.006, 0.0312];
color33 = [0.0313, 0.01, 0];
color40 = [
%}

% 1.  Set the top colors to orange

vals = linspace(0, 1, numHalfColors+1);

%valTop = 0.5*(vals + vals.^exponent);
valTop = vals.^exponent;

valTop = valTop(2:end);

colormap(vv > 0.5, 1) = valTop;
colormap(vv > 0.5, 2) = 0.3 * valTop.^2;

%colormap(topHalf,1) = (0:127)/127;
%%colormap(topHalf,2) = 0.2720 * ( ((0:31)/31) - 0.4*((0:31)/31).^3 )  ;
%colormap(topHalf,2) = 0.3 * ( ((0:127)/127) ).^2 ;
%%colormap(topHalf, 3) = 0.3*((0:31)/31).^3;


% 2.  Set the bottom colors to blue

valBot = valTop(end:-1:1);

colormap(vv < 0.5,1) = 0.4118 * valBot;
colormap(vv < 0.5,2) = 0.1922 * valBot.^2;
colormap(vv < 0.5,3) = valBot;

%colormap(bottomHalf,1) = 0.4118 * ((127:-1:0)/127);
%colormap(bottomHalf,2) = 0.1922 * ((127:-1:0)/127).^2;
%colormap(bottomHalf,3) = (127:-1:0)/127;

%colormap = colormap .^ exponent;

scale = [valBot 0 valTop];