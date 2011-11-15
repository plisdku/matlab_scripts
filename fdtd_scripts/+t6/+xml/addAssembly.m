function addAssembly(assembly, gridXML, documentNode, origin)
global TROG_XML_COUNT___;
doc = documentNode;

digits = 14;

assemblyXML = doc.createElement('Assembly');

for aa = 1:length(assembly)
    
    switch assembly{aa}.type
        case 'Mesh'
            elemXML = doc.createElement('Mesh');
            
            if isfield(assembly{aa}, 'permittivity')
                elemXML.setAttribute('permittivity', ...
                    assembly{aa}.permittivity);
            end
            if isfield(assembly{aa}, 'permeability')
                elemXML.setAttribute('permeability',...
                    assembly{aa}.permeability);
            end
            
            for vv = 1:length(assembly{aa}.vertices)
                vertXML = doc.createElement('Vertex');
                vertXML.setAttribute('position', ...
                    num2str(assembly{aa}.vertices(vv,:) - origin, digits));
                vertXML.setAttribute('freeDirections', ...
                    num2str(assembly{aa}.vertexFreeDirections(vv,:)));
                elemXML.appendChild(vertXML);
            end
            
            for ff = 1:length(assembly{aa}.faces)
                faceXML = doc.createElement('Face');
                faceXML.setAttribute('vertices', ...
                    num2str(assembly{aa}.faces(ff,:), digits));
                elemXML.appendChild(faceXML);
            end
            
            assemblyXML.appendChild(elemXML);
        
        case 'Background'
            elemXML = doc.createElement('Background');
            if isfield(assembly{aa}, 'permittivity')
                elemXML.setAttribute('permittivity', ...
                    assembly{aa}.permittivity);
            end
            if isfield(assembly{aa}, 'permeability')
                elemXML.setAttribute('permeability', ...
                    assembly{aa}.permeability);
            end
            assemblyXML.appendChild(elemXML);
            
    end % switch assembly instruction type
end % foreach assembly instruction

gridXML.appendChild(assemblyXML);