function renderGeometry(varargin)
% renderGeometry('Model', m, 'Parameters', p, 'MeasurePosition', mp, ...
%   'SourcePosition', sp, 'Units', 1, 'Sensitivity', sen);
X.Model = [];
X.Parameters = [];
X.MeasurePosition = [];
X.SourcePosition = [];
X.Units = 1;
X.Sensitivity = [];

X = parseargs(X, varargin{:});

theMesh = Main();

meshes = X.Model.meshes(X.Parameters);

controlVertices = [];
for mm = 1:length(meshes)
    controlVertices = [controlVertices; meshes{mm}.patchVertices];
end

plot3(X.Units*controlVertices(:,1), X.Units*controlVertices(:,2), ...
    X.Units*controlVertices(:,3), 'bo')
hold on

matlColors = {'g', 'k', 'g', 'y'};
%if (numel(theMesh.permittivity) > 1)
    
    for pp = 1:numel(theMesh.permittivity)
        patch('Vertices', X.Units*theMesh.permittivity{pp}.vertices, ...
            'Faces', theMesh.permittivity{pp}.faces,...
            'FaceColor', matlColors{pp}, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.02)
    end
    view(2) % top-down
%end


if ~isempty(X.Sensitivity)
    quiver3(X.Units*controlVertices(:,1), X.Units*controlVertices(:,2),...
        controlVertices(:,3), X.Sensitivity(1,1:3:end), ...
        X.Sensitivity(1,2:3:end), X.Sensitivity(1,3:3:end), ...
        'r', 'LineWidth', 1.5);
end

if ~isempty(X.MeasurePosition)
    plot3(X.Units*X.MeasurePosition(1), X.Units*X.MeasurePosition(2), ...
        X.Units*X.MeasurePosition(3)+0.5, 'go', 'MarkerSize', 10);
end

if ~isempty(X.SourcePosition)
    plot3(X.Units*X.SourcePosition(1), X.Units*X.SourcePosition(2), ...
        X.Units*X.SourcePosition(3)+0.5, 'rx', 'MarkerSize', 10);
end

end

