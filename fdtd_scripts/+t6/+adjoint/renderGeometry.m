function renderGeometry(group, parameters, measurePosition, srcPosition,...
    units, dfdv)
% renderGeometry(group, parameters, measurePosition, srcPosition)

theMesh = Main();

meshes = group.meshes(parameters);

if nargin < 5
    units = 1;
end

controlVertices = [];
for mm = 1:length(meshes)
    controlVertices = [controlVertices; meshes{mm}.patchVertices];
end

plot3(units*controlVertices(:,1), units*controlVertices(:,2), controlVertices(:,3),...
    'bo')
%axis image
hold on

matlColors = {'g', 'k', 'g', 'y'};
if (numel(theMesh.permittivity) > 1)
    
    %for pp = 3:length(theMesh.permittivity)
    for pp = 1
        patch('Vertices', units*theMesh.permittivity{pp}.vertices, ...
            'Faces', theMesh.permittivity{pp}.faces,...
            'FaceColor', matlColors{pp}, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.02)
    end
    view(2) % top-down
%hold on
end


if nargin == 6
    quiver3(units*controlVertices(:,1), units*controlVertices(:,2), controlVertices(:,3),...
        dfdv(1,1:3:end), dfdv(1,2:3:end), dfdv(1,3:3:end), 'r', ...
        'LineWidth', 1.5);
end
%xlabel('x (cells)')
%ylabel('y (cells)')

if exist('measurePosition')
plot3(units*measurePosition(1), units*measurePosition(2), units*measurePosition(3)+0.5, ...
    'go', 'MarkerSize', 10);
end

%if exist('srcPosition')
%plot3(units*srcPosition(1), units*srcPosition(2), units*srcPosition(3)+0.5, 'rx', ...
%    'MarkerSize', 10)
end
%trogRect(totalFieldCells, 'y', 'LineStyle', '-.')


%axis([25 26 0 10 0 1])
%axis vis3d
%view(3)
%axis([25 26 4.5 5.5])
%axis equal
%view(3)
%camlight right