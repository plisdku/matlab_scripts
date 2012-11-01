function [rgb, maxVal] = complex2rgb(A, varargin)

X.Mode = 'value';
X.Threshold = inf;
X = parseargs(X, varargin{:});

colVec = @(v) reshape(v, [], 1);

maxVal = max(abs(colVec(A)));
if maxVal > X.Threshold
    maxVal = X.Threshold;
end

if strcmpi(X.Mode, 'value')
    hue = (pi+angle(colVec(A)))/2/pi;
    saturation = ones(size(hue));
    value = min(1, abs(colVec(A)) / maxVal);
elseif strcmpi(X.Mode, 'saturation')
    hue = (pi+angle(colVec(A)))/2/pi;
    saturation = min(1, abs(colVec(A)) / maxVal);
    value = ones(size(saturation));
elseif strcmpi(X.Mode, 'flat')
    hue = (pi+angle(colVec(A)))/2/pi;
    saturation = ones(size(hue));
    value = saturation;
end
    
rgb = reshape(hsv2rgb(cat(2, hue, saturation, value)), [size(A), 3]);
