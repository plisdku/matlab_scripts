function addTFSFSources(grid, gridXML, doc, originTrogdor)
global TROG_XML_COUNT___;
originTwice = [originTrogdor originTrogdor];
precisionString = t6.TrogdorSimulation.instance().Precision;

for ll = 1:length(grid.TFSFSources)
    src = grid.TFSFSources{ll};
    
    elemXML = doc.createElement('TFSFSource');
    
    fieldstr = '';
    for ff = 1:length(src.field)
        fieldstr = [fieldstr, src.field{ff}, ' '];
    end
    elemXML.setAttribute('fields', fieldstr);
    
    elemXML.setAttribute('yeeCells', ...
        sprintf('%i ', src.yeeCells + originTwice));
    elemXML.setAttribute('direction', ...
        sprintf('%i ', src.direction));
    
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
    
    %dataFileName = t6.xml.randomName('__tfsfsource_', '', 8);
    dataFileName = sprintf('__tfsfsource_%i', TROG_XML_COUNT___.tfsfTime);
    TROG_XML_COUNT___.tfsfTime = TROG_XML_COUNT___.tfsfTime + 1;
    elemXML.setAttribute('file', dataFileName);
    fh = fopen(dataFileName, 'w');
    try
        % Observe carefully: the array needs to be written in row order,
        % so I need to transpose it.  Matlab writes 2D arrays in column order
        % (column 1, then column 2, etc.).
        if size(src.timeData, 1) == 1
            timeData = repmat(src.timeData, [length(src.field), 1]);
        else
            timeData = src.timeData;
        end
        
        count = fwrite(fh, timeData, precisionString);
    catch
        error('Could not write TFSF source data file.');
    end
    fclose(fh);    
    t6.xml.writeSourceSpec(src, 'AutoTimeFile', dataFileName);
    
    % Count up the fields and decide if this source needs a polarization vector
    % or not.
    fieldTokens = {};
    remainder = src.field;
    while ~strcmp(remainder, '')
        [token, remainder] = strtok(remainder);
        if ~strcmp(token, '')
            fieldTokens = {fieldTokens{:}, token};
        end
    end
    
    if length(fieldTokens) > size(src.timeData, 1)
        elemXML.setAttribute('polarization', '1 1 1');
    end
    
    gridXML.appendChild(elemXML);
end
