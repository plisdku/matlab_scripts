function addAssembly(assembly, gridXML, documentNode, originTrogdor)
global TROG_XML_COUNT___;
doc = documentNode;

originTwice = [originTrogdor originTrogdor];

digits = 10;

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
                    num2str(assembly{aa}.vertices(vv,:), digits));
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
            
        case 'Block'
            elemXML = doc.createElement('Block');
            elemXML.setAttribute('yeeBounds', ...
                sprintf('%i ', assembly{aa}.yeeBounds + originTwice));
            if isfield(assembly{aa}, 'permittivity')
                elemXML.setAttribute('permittivity', ...
                    assembly{aa}.permittivity);
            end
            if isfield(assembly{aa}, 'permeability')
                elemXML.setAttribute('permeability',...
                    assembly{aa}.permeability);
            end
            assemblyXML.appendChild(elemXML);
            
        case 'Ellipsoid'
            elemXML = doc.createElement('Ellipsoid');
            elemXML.setAttribute('yeeBounds', ...
                sprintf('%i ', assembly{aa}.yeeBounds + originTwice));
            if isfield(assembly{aa}, 'permittivity')
                elemXML.setAttribute('permittivity', ...
                    assembly{aa}.permittivity);
            end
            if isfield(assembly{aa}, 'permeability')
                elemXML.setAttribute('permeability',...
                    assembly{aa}.permeability);
            end
            assemblyXML.appendChild(elemXML);
        
        case 'HeightMap'
            elemXML = doc.createElement('HeightMap');
            elemXML.setAttribute('yeeBounds', ...
                sprintf('%i ', assembly{aa}.yeeBounds + originTwice));
            if isfield(assembly{aa}, 'permittivity')
                elemXML.setAttribute('permittivity', ...
                    assembly{aa}.permittivity);
            end
            if isfield(assembly{aa}, 'permeability')
                elemXML.setAttribute('permeability',...
                    assembly{aa}.permeability);
            end
            elemXML.setAttribute('row', assembly{aa}.row);
            elemXML.setAttribute('column', assembly{aa}.column);
            elemXML.setAttribute('up', assembly{aa}.up);
            
            %imfilename = t6.xml.randomName('__heightmap_', '.bmp', 8);
            imfilename = sprintf('__heightmap_%i.bmp', ...
                TROG_XML_COUNT___.heightMap);
            TROG_XML_COUNT___.heightMap = TROG_XML_COUNT___.heightMap + 1;
            elemXML.setAttribute('file', imfilename);
            imwrite(assembly{aa}.image, imfilename, 'bmp');
            assemblyXML.appendChild(elemXML);
    end % switch assembly instruction type
end % foreach assembly instruction

gridXML.appendChild(assemblyXML);