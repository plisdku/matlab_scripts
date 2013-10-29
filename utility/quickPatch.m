function quickPatch(v, f, colorString)

if nargin == 2
    colorString = 'r';
end

patch('Vertices', v, 'Faces', f, 'FaceColor', colorString, ...
    'FaceAlpha', 0.5);