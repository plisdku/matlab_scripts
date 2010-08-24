function addAssembly(assembly, gridXML, documentNode, originTrogdor)
global TROG_XML_COUNT___;
doc = documentNode;

originTwice = [originTrogdor originTrogdor];

assemblyXML = doc.createElement('Assembly');

for aa = 1:length(assembly)
    
    switch assembly{aa}.type
        case 'Block'
            elemXML = doc.createElement('Block');
            elemXML.setAttribute('yeeCells', ...
                sprintf('%i ', assembly{aa}.yeeCells + originTwice));
            elemXML.setAttribute('material', assembly{aa}.materialName);
            if isfield(assembly{aa}, 'fillStyle')
                elemXML.setAttribute('fillStyle', assembly{aa}.fillStyle);
            end
            assemblyXML.appendChild(elemXML);
            
        case 'Ellipsoid'
            elemXML = doc.createElement('Ellipsoid');
            elemXML.setAttribute('yeeCells', ...
                sprintf('%i ', assembly{aa}.yeeCells + originTwice));
            elemXML.setAttribute('material', assembly{aa}.materialName);
            if isfield(assembly{aa}, 'fillStyle')
                elemXML.setAttribute('fillStyle', assembly{aa}.fillStyle);
            end
            assemblyXML.appendChild(elemXML);
            
        case 'KeyImage'
            elemXML = doc.createElement('KeyImage');
            elemXML.setAttribute('yeeCells', ...
                sprintf('%i ', assembly{aa}.yeeCells + originTwice));
            elemXML.setAttribute('row', assembly{aa}.row);
            elemXML.setAttribute('column', assembly{aa}.column);
            %imfilename = t5.xml.randomName('__keyimage_', '.bmp', 8);
            imfilename = sprintf('__keyimage_%i.bmp', ...
                TROG_XML_COUNT___.keyImage);
            TROG_XML_COUNT___.keyImage = TROG_XML_COUNT___.keyImage + 1;
            elemXML.setAttribute('file', imfilename);
            
            COLOR = ndims(assembly{aa}.image) == 3;
            GREYSCALE = ndims(assembly{aa}.image) == 2;
            
            if COLOR
                imageToWrite = assembly{aa}.image;
                for nn = 1:length(assembly{aa}.tags)
                    keyTag = doc.createElement('Tag');
                    tag = assembly{aa}.tags{nn};
                    tagPixel = uint8(round(tag.pixel * 255) );
                    assert(length(tagPixel) == 3);
                    hexred = dec2hex(tagPixel(1), 2);
                    hexgreen = dec2hex(tagPixel(2), 2);
                    hexblue = dec2hex(tagPixel(3), 2);
                    keyTag.setAttribute('color', ...
                        ['#', hexred, hexgreen, hexblue]);
                    keyTag.setAttribute('material', tag.material);
                    keyTag.setAttribute('fillStyle', tag.fillStyle);
                    elemXML.appendChild(keyTag);
                end
            elseif GREYSCALE
                imageToWrite = repmat(assembly{aa}.image, [1 1 3]); % RGB
                for nn = 1:length(assembly{aa}.tags)
                    keyTag = doc.createElement('Tag');
                    tag = assembly{aa}.tags{nn};
                    tagPixel = uint8(round(tag.pixel * 255) );
                    tagPixel = [tagPixel tagPixel tagPixel];
                    %assert(length(tagPixel) == 1);
                    hexred = dec2hex(tagPixel(1), 2);
                    hexgreen = dec2hex(tagPixel(2), 2);
                    hexblue = dec2hex(tagPixel(3), 2);
                    keyTag.setAttribute('color', ...
                        ['#', hexred, hexgreen, hexblue]);
                    keyTag.setAttribute('material', tag.material);
                    keyTag.setAttribute('fillStyle', tag.fillStyle);
                    elemXML.appendChild(keyTag);
                end
            else
                error('KeyImage has %i dimensions.', ndims(assembly{aa}.image));
            end
            imwrite(imageToWrite, imfilename, 'bmp');
            assemblyXML.appendChild(elemXML);
            
        case 'HeightMap'
            elemXML = doc.createElement('HeightMap');
            elemXML.setAttribute('yeeCells', ...
                sprintf('%i ', assembly{aa}.yeeCells + originTwice));
            elemXML.setAttribute('material', assembly{aa}.materialName);
            elemXML.setAttribute('row', assembly{aa}.row);
            elemXML.setAttribute('column', assembly{aa}.column);
            elemXML.setAttribute('up', assembly{aa}.up);
            if isfield(assembly{aa}, 'fillStyle')
                elemXML.setAttribute('fillStyle', assembly{aa}.fillStyle);
            end
            %imfilename = t5.xml.randomName('__heightmap_', '.bmp', 8);
            imfilename = sprintf('__heightmap_%i.bmp', ...
                TROG_XML_COUNT___.heightMap);
            TROG_XML_COUNT___.heightMap = TROG_XML_COUNT___.heightMap + 1;
            elemXML.setAttribute('file', imfilename);
            imwrite(assembly{aa}.image, imfilename, 'bmp');
            assemblyXML.appendChild(elemXML);
    end % switch assembly instruction type
end % foreach assembly instruction

gridXML.appendChild(assemblyXML);