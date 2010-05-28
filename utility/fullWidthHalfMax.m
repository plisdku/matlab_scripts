function [fwhm, edges] = fullWidthHalfMax(x, y, varargin)

if nargin > 2
    peakIndex = varargin{1};
    peakVal = y(peakIndex);
else
    [peakVal, peakIndex] = max(y);
end

left = find(y(1:peakIndex) < peakVal/2, 1, 'last');
right = peakIndex - 1 + find(y(peakIndex:end) < peakVal/2, 1, 'first');

if isempty(left) || isempty(right)
    error('Cannot find half max positions');
end

xLeft = interp1(y([left, left+1]), x([left, left+1]), peakVal/2);
xRight = interp1(y([right, right-1]), x([right, right-1]), peakVal/2);

fwhm = xRight - xLeft;
edges = [xLeft, xRight];

