function flatPatch(varargin)
% flatPatch  Replacement for patch() with surface normals instead of vertex
% normals
%
% Usage:
%
% flatPatch('Vertices', v, 'Faces', f);
%
% See the help for patch() for details.  :-D
%

X.Faces = [];
X.Vertices = [];
X.FaceColor = [];
X.FaceAlpha = [];
X.EdgeColor = [];
X.EdgeAlpha = [];
X.FaceVertexCData = [];
X = parseargs(X, varargin{:});

allFields = fieldnames(X);

for ff = 1:length(allFields)
    
    if isempty(X.(allFields{ff}))
        X = rmfield(X, allFields{ff});
    end
end

[v f] = disconnect(X.Vertices, X.Faces);

X.Vertices = v;
X.Faces = f;

argList = {};
allfields = fieldnames(X);
for ff = 1:length(allfields)
    argList{end+1} = allfields{ff};
    argList{end+1} = X.(allfields{ff});
end

patch(argList{:});



function [v f] = disconnect(vertices, faces)

v = zeros(3*size(faces,1), 3);
f = reshape(1:3*size(faces,1), 3, [])';

for ff = 1:size(faces,1)
    v(1+3*(ff-1),:) = vertices(faces(ff,1),:);
    v(2+3*(ff-1),:) = vertices(faces(ff,2),:);
    v(3+3*(ff-1),:) = vertices(faces(ff,3),:);
end