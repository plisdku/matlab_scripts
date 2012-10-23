function rgb = complex2rgb(A, varargin)

X.Mode = 'value';
X = parseargs(X, varargin{:});

colVec = @(v) reshape(v, [], 1);

if strcmpi(X.Mode, 'value')
    hue = (pi+angle(colVec(A)))/2/pi;
    saturation = ones(size(hue));
    value = abs(colVec(A)) / max(abs(colVec(A)));
elseif strcmpi(X.Mode, 'saturation')
    hue = (pi+angle(colVec(A)))/2/pi;
    saturation = abs(colVec(A)) / max(abs(colVec(A)));
    value = ones(size(saturation));
elseif strcmpi(X.Mode, 'flat')
    hue = (pi+angle(colVec(A)))/2/pi;
    saturation = ones(size(hue));
    value = saturation;
end
    
rgb = reshape(hsv2rgb(cat(2, hue, saturation, value)), [size(A), 3]);
