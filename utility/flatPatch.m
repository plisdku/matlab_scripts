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
% Extra options:
%
% FaceFilter    a function @(vert1, vert2, vert3) that must be true for the
%               face to be rendered

% turned off...

% AllVertices   a function @(vertex) that must be true for all vertices to
%               display a face, otherwise it's culled.
%
% AnyVertex     a function @(vertex) that must be true for any vertex in a
%               face for it to be rendered, otherwise it's culled
%

X.Faces = [];
X.Vertices = [];
X.FaceColor = [];
X.FaceAlpha = [];
X.EdgeColor = [];
X.EdgeAlpha = [];
X.FaceVertexCData = [];
%X.AllVertices = [];
%X.AnyVertex = [];
X.FaceFilter = [];
X = parseargs(X, varargin{:});

%allVertices = X.AllVertices; X.AllVertices = [];
%anyVertex = X.AnyVertex; X.AnyVertex = [];
faceFilter = X.FaceFilter; X.FaceFilter = [];

allFields = fieldnames(X);
for ff = 1:length(allFields)
    
    if isempty(X.(allFields{ff}))
        X = rmfield(X, allFields{ff});
    end
end

% Filter the faces:

if isa(faceFilter, 'function_handle')
    
    faceFlags = true(size(X.Faces,1),1);
    for ff = 1:numel(faceFlags)
        faceFlags(ff) = faceFilter(X.Vertices(X.Faces(ff,1),:), ...
            X.Vertices(X.Faces(ff,2),:), ...
            X.Vertices(X.Faces(ff,3),:));
    end
    
    X.Faces = X.Faces(faceFlags,:);
end

%{
if isa(allVertices, 'function_handle')
    
    for ff = 1:numel(faceFlags)
        faceFlags(ff) = faceFlags(ff) & ...
            allVertices(X.Vertices(X.Faces(ff,1),:)) & ...
            allVertices(X.Vertices(X.Faces(ff,2),:)) & ...
            allVertices(X.Vertices(X.Faces(ff,3),:));
    end
    
end

if isa(anyVertex, 'function_handle')
    
    for ff = 1:numel(faceFlags)
        faceFlags(ff) = faceFlags(ff) & ( ...
            anyVertex(X.Vertices(X.Faces(ff,1),:)) | ...
            anyVertex(X.Vertices(X.Faces(ff,2),:)) | ...
            anyVertex(X.Vertices(X.Faces(ff,3),:)) );
    end
end
%}

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
