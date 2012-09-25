function addCustomTFSFSources(sim, gridXML, doc, mode)

for ll = 1:length(sim.Grid.CustomTFSFSources)
    src = sim.Grid.CustomTFSFSources{ll};
    
    if ~isempty(src.mode) && ~strcmpi(src.mode, mode)
        continue;
    end
    
    elemXML = doc.createElement('CustomTFSFSource');
    elemXML.setAttribute('yeeCells', ...
        sprintf('%i ', src.yeeCells));
    elemXML.setAttribute('file', src.spaceTimeFile);
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
    
    t6.xml.writeSourceSpec(sim, src);
    
    gridXML.appendChild(elemXML);
end
