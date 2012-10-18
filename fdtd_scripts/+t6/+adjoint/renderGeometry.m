function renderGeometry(varargin)
% renderGeometry('Model', m, 'Parameters', p, 'MeasurePosition', mp, ...
%   'SourcePosition', sp, 'Units', 1, 'Sensitivity', sen, ...
%   'Materials', 1:4);
X.Model = [];
X.Parameters = [];
X.MeasurePosition = [];
X.SourcePosition = [];
X.Units = 1;
X.Materials = [];
X.Sensitivity = [];
X.Mesh = [];
X.FaceAlpha = 0.2;
X.EdgeAlpha = 0.02;

X = parseargs(X, varargin{:});

if isempty(X.Mesh)
    theMesh = Main();
else
    func = inline(sprintf('%s()', X.Mesh));
    theMesh = func([]);
end

if isempty(X.Materials)
    X.Materials = 1:numel(theMesh.permittivity);
else
    if any(X.Materials < 1) || any(X.Materials > numel(theMesh.permittivity))
        error('Materials out of bounds, only %i present in simulation', ...
            numel(theMesh.permittivity));
    end
end

if ~isempty(X.Model)
    meshes = X.Model.meshes(X.Parameters);

    controlVertices = [];
    for mm = 1:length(meshes)
        controlVertices = [controlVertices; meshes{mm}.patchVertices];
    end
    
    plot3(X.Units*controlVertices(:,1), X.Units*controlVertices(:,2), ...
        X.Units*controlVertices(:,3), 'bo')
    hold on

end


matlColors = {'g', 'b', 'g', 'y', 'g', 'b', 'y'};
%if (numel(theMesh.permittivity) > 1)
    
    for pp = X.Materials
        flatPatch('Vertices', X.Units*theMesh.permittivity{pp}.vertices, ...
            'Faces', theMesh.permittivity{pp}.faces,...
            'FaceColor', matlColors{pp}, 'FaceAlpha', X.FaceAlpha, ...
            'EdgeAlpha', X.EdgeAlpha)
    end
    view(2) % top-down
    hold on
%end


if ~isempty(X.Sensitivity)
    quiver3(controlVertices(:,1)', controlVertices(:,2)',...
        controlVertices(:,3)', ...
        X.Sensitivity(1,1:3:end), X.Sensitivity(1,2:3:end), ...
        X.Sensitivity(1,3:3:end), ...
        'r', 'LineWidth', 1.5);
    hold on
end

if ~isempty(X.MeasurePosition)
    plot3(X.Units*X.MeasurePosition(1), X.Units*X.MeasurePosition(2), ...
        X.Units*X.MeasurePosition(3)+0.5, 'go', 'MarkerSize', 10);
    hold on
end

if ~isempty(X.SourcePosition)
    plot3(X.Units*X.SourcePosition(1), X.Units*X.SourcePosition(2), ...
        X.Units*X.SourcePosition(3)+0.5, 'rx', 'MarkerSize', 10);
    hold on
end

end

