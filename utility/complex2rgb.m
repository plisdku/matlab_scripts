function rgb = complex2rgb(A, valueFunc)

colVec = @(v) reshape(v, [], 1);

hue = (pi+angle(colVec(A)))/2/pi;
saturation = ones(size(hue));
value = abs(colVec(A)) / max(abs(colVec(A)));
%value = saturation;

%saturation = value;
%value = ones(size(value));

rgb = reshape(hsv2rgb(cat(2, hue, saturation, value)), [size(A), 3]);
