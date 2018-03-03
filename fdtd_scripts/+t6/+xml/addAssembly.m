function addAssembly(sim, gridXML, documentNode, origin, designParameters)
global TROG_XML_COUNT___;
doc = documentNode;

assemblyXML = doc.createElement('Assembly');

allVertices = [];
allFreeDirections = [];
idxNextControlVert = 0;

for aa = 1:length(sim.Grid.Meshes)
    
    mesh = sim.Grid.Meshes{aa};
    
    vertices = sim.extendIntoPML(mesh.patchVertices);
    faces = mesh.faces-1;
    freeDirections = mesh.freeDirections();
    if isempty(freeDirections)
        freeDirections = zeros(size(vertices));
    end
    
    allVertices = [allVertices; vertices];
    allFreeDirections = [allFreeDirections; freeDirections];
    
    % Create the old-style mesh XML
    elemXML = meshXML(doc, origin, mesh.permittivity, mesh.permeability, ...
        vertices, faces, freeDirections);
    assemblyXML.appendChild(elemXML);
    
    % Add the new-style mesh XML
    elemXML = meshFaceXML(doc, mesh.permittivity, mesh.permeability, ...
        faces + idxNextControlVert);
    idxNextControlVert = idxNextControlVert + size(vertices, 1);
    assemblyXML.appendChild(elemXML);
end

% Add the new-style vertex XML

elemXML = verticesXML(doc, origin, allVertices, allFreeDirections);
assemblyXML.appendChild(elemXML);

% Background!

if ~isempty(sim.Grid.Background)
    elemXML = doc.createElement('Background');
    if isfield(sim.Grid.Background, 'permittivity')
        elemXML.setAttribute('permittivity', ...
            sim.Grid.Background.permittivity);
    end
    if isfield(sim.Grid.Background, 'permeability')
        elemXML.setAttribute('permeability', ...
            sim.Grid.Background.permeability);
    end
    assemblyXML.appendChild(elemXML);
end




gridXML.appendChild(assemblyXML);


%function elemXML = meshXML(doc, mesh)
function elemXML = verticesXML(doc, origin, vertices, freeDirections)

digits = 14;

elemXML = doc.createElement('Vertices');

for vv = 1:length(vertices)
    vertXML = doc.createElement('Vertex');
    vertXML.setAttribute('position', num2str(vertices(vv,:) - origin, digits));
    vertXML.setAttribute('freeDirections', num2str(freeDirections(vv,:)));
    elemXML.appendChild(vertXML);
end



%function elemXML = meshXML(doc, mesh)
function elemXML = meshFaceXML(doc, permittivity, permeability, faces)

digits = 14;

elemXML = doc.createElement('NewMesh');

if ~isempty(permittivity)
    elemXML.setAttribute('permittivity', permittivity);
end

if ~isempty(permeability)
    elemXML.setAttribute('permeability', permeability);
end


for ff = 1:length(faces)
    faceXML = doc.createElement('Face');
    faceXML.setAttribute('vertices', num2str(faces(ff,:), digits));
    elemXML.appendChild(faceXML);
end




%function elemXML = meshXML(doc, mesh)
function elemXML = meshXML(doc, origin, permittivity, permeability, vertices,...
    faces, freeDirections)

digits = 14;

elemXML = doc.createElement('Mesh');

if ~isempty(permittivity)
    elemXML.setAttribute('permittivity', permittivity);
end

if ~isempty(permeability)
    elemXML.setAttribute('permeability', permeability);
end

for vv = 1:length(vertices)
    vertXML = doc.createElement('Vertex');
    vertXML.setAttribute('position', num2str(vertices(vv,:) - origin, digits));
    vertXML.setAttribute('freeDirections', num2str(freeDirections(vv,:)));
    elemXML.appendChild(vertXML);
end

for ff = 1:length(faces)
    faceXML = doc.createElement('Face');
    faceXML.setAttribute('vertices', num2str(faces(ff,:), digits));
    elemXML.appendChild(faceXML);
end
