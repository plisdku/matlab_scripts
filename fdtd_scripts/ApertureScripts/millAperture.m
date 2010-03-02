function heightmap = millAperture(millPattern, millTime, millGaussianSigma,...
    thicknesses, millRates)
% millAperture   Make heightmaps for simulating realistic milled structures
%    heightMap = millAperture(millPattern, millTime, millGaussianSigma, ...
%        thicknesses, millRates)
% This function will return a quasi-realistic heightmap such as you might
% get from milling the specified millPattern with a focused ion beam with
% gaussian beam profile.  The sigma specified is in pixels but can be
% fractional.  Provide an array of layer thicknesses and mill rates to
% simulate milling through structures with different material properties.
% The returned heightmap is not normalized to [0 1].
%
% Example:
% 
% aperture = imread('c_large.bmp');  % black is no mill, white is mill
% rates = [1 8];   % Si3N4 and gold mill rates
% time = 100;  % will penetrate 100 cells of Si3N4 (millRate = 1)
% sigma = 3;    % pixels
% stack = [65 100];  % 65 cells Si3N4, 100 cells gold
%
% heightmap = millAperture(aperture, time, sigma, stack, rates);
% 

% make a nice milled aperture heightmap

stackLevels = [0, cumsum(thicknesses)];

filterWidth = 1 + 2*round(3*millGaussianSigma/2);  % odd numbers best.
gaussFilt = fspecial('gaussian', filterWidth, millGaussianSigma);


%%
millExposures = millTime*imfilter(double(millPattern), gaussFilt, 'replicate');
%%
for (nStack = 1:length(millRates))

    sel = find(millExposures >= stackLevels(nStack));
    overshoot = millExposures(sel) - stackLevels(nStack);
    millExposures(sel) = stackLevels(nStack) + ...
        overshoot*millRates(nStack);
    sel = find(millExposures >= stackLevels(nStack+1));
    overshoot = millExposures(sel) - stackLevels(nStack+1);
    millExposures(sel) = stackLevels(nStack+1) + ...
        overshoot/millRates(nStack);
end
millExposures(find(millExposures > stackLevels(end))) = ...
    stackLevels(end);
%% Color each exposed material differently
materials = zeros(size(millExposures));

%%
for (nStack = 1:length(thicknesses))
    sel = find(millExposures >= stackLevels(nStack) & ...
        millExposures < stackLevels(nStack+1));
    materials(sel) = nStack;
end

heightmap = millExposures;
%heightmap = millExposures/stackLevels(end);

%{
%%
figure(1)
imagesc(millExposures*1e9, ...
    [0, stackLevels(end)*1e9]);
prettify('Times', 18);
colormap gray
colorbar
xlabel('X (nm)');
ylabel('Y (nm)');
title(sprintf('Mill depth for mill time %2.2f (depth in nm)', ...
    millTime*1e9));
savePlot(sprintf('milldepth_%2.2f', millTime*1e9));

%%
figure(2)
imagesc(materials);
prettify('Times', 18);
colormap jet
colorbar
xlabel('X (nm)');
ylabel('Y (nm)');
title(sprintf('Mill layers for mill time %2.2f nm', millTime*1e9));
savePlot(sprintf('materials_%2.2f', millTime*1e9));

%%

xmedian = round(size(millExposures, 2)/2);
ymedian = round(size(millExposures, 1)/2);

yCut = millExposures(ymedian,:);
%yCut(find(yCut > stackLevels(end))) = Inf;
xCut = millExposures(:,xmedian);
%xCut(find(xCut > stackLevels(end))) = Inf;

figure(3)
clf
plot(xCut*1e9);
prettify('Times', 18);
hold on
for (nn = 1:length(stackLevels))
    plot([ys(1), ys(end)]*1e9, [1,1]*stackLevels(nn)*1e9, 'k--', ...
        'LineWidth', 2);
end
xlabel('Y (nm)');
ylabel('Z (nm)');
title(sprintf('X-axis mill cross-section for time %2.2f nm', ...
    millTime*1e9));
%savePlot(sprintf('crossx_%2.2f', millTime*1e9));

figure(4)
clf
plot(yCut*1e9);
prettify('Times', 18);
hold on
for (nn = 1:length(stackLevels))
    plot([xs(1), xs(end)]*1e9, [1,1]*stackLevels(nn)*1e9, 'k--', ...
        'LineWidth', 2);
end
xlabel('X (nm)');
ylabel('Z (nm)');
title(sprintf('Y-axis mill cross-section for time %2.2f nm', ...
    millTime*1e9));
%savePlot(sprintf('crossy_%2.2f', millTime*1e9));

imwrite(millExposures/stackLevels(end), ...
    sprintf('aperture_time%2.2f.bmp', millTime*1e9));
%}
