function addCustomTFSFSources(grid, gridXML, doc, originTrogdor)

originTwice = [originTrogdor originTrogdor];
directory = t5.TrogdorSimulation.instance().directoryString;

for ll = 1:length(grid.CustomTFSFSources)
    src = grid.CustomTFSFSources{ll};
    
    elemXML = doc.createElement('CustomTFSFSource');
    elemXML.setAttribute('yeeCells', ...
        sprintf('%i ', src.yeeCells + originTwice));
    elemXML.setAttribute('file', [directory, src.spaceTimeFile]);
    elemXML.setAttribute('symmetries', sprintf('%i ', src.symmetries));
    
    for oo = 1:length(src.omitSides)
        omitXML = doc.createElement('OmitSide');
        omitXML.appendChild(doc.createTextNode(...
            sprintf('%i ', src.omitSides{oo})));
        elemXML.appendChild(omitXML);
    end
    
    for dd = 1:size(src.duration, 1)
        durXML = doc.createElement('Duration');
        durXML.setAttribute('firstTimestep', num2str(src.duration(dd,1)));
        durXML.setAttribute('lastTimestep', num2str(src.duration(dd,2)));
        elemXML.appendChild(durXML);
    end
    
    t5.xml.writeSourceSpec(src);
    
    gridXML.appendChild(elemXML);
end
