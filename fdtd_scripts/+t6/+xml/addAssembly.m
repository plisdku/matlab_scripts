function addAssembly(sim, gridXML, documentNode, origin, designParameters)
global TROG_XML_COUNT___;
doc = documentNode;

assemblyXML = doc.createElement('Assembly');

for aa = 1:length(sim.Grid.Meshes)
    
    mesh = sim.Grid.Meshes{aa};
    
    vertices = sim.extendIntoPML(mesh.patchVertices);
    faces = mesh.faces-1;
    freeDirections = mesh.freeDirections();
    if isempty(freeDirections)
        freeDirections = zeros(size(vertices));
    end
    
    elemXML = meshXML(doc, origin, mesh.permittivity, mesh.permeability, ...
        vertices, faces, freeDirections);
    
    assemblyXML.appendChild(elemXML);
end

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
