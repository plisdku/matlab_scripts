function addTFSFSources(grid, gridXML, doc, origin)
global TROG_XML_COUNT___;

directory = t6.TrogdorSimulation.instance().directoryString;

for ll = 1:length(grid.TFSFSources)
    src = grid.TFSFSources{ll};
    
    elemXML = doc.createElement('TFSFSource');
    
    fieldstr = '';
    for ff = 1:length(src.field)
        fieldstr = [fieldstr, src.field{ff}, ' '];
    end
    elemXML.setAttribute('fields', fieldstr);
    
    elemXML.setAttribute('yeeCells', ...
        sprintf('%i ', src.yeeCells));
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
    
    dataFileName = [directory, sprintf('__tfsfsource_%i', TROG_XML_COUNT___.tfsfTime)];
    TROG_XML_COUNT___.tfsfTime = TROG_XML_COUNT___.tfsfTime + 1;
    elemXML.setAttribute('file', dataFileName);
    
    if ~isempty(src.timeData)
        myWriteTimeData(dataFileName, src.timeData, numel(src.field));
    elseif ~isempty(src.fieldFunction)
        myWriteFunction(dataFileName, src.fieldFunction, src.field, ...
            src.duration) 
    end
    
    t6.xml.writeSourceSpec(src, 'AutoTimeFile', dataFileName);
    
    gridXML.appendChild(elemXML);
end

end % main function


function myWriteTimeData(dataFileName, timeData, numFields)
precisionString = t6.TrogdorSimulation.instance().Precision;
fh = fopen(dataFileName, 'w');

try
    if numFields > 1 && size(timeData, 1) == 1 % if it's a row vector
        timeData = repmat(timeData, [numFields, 1]);
    end
    
    count = fwrite(fh, timeData, precisionString);
catch
    error('Could not write TFSF source data file.');
end
fclose(fh);
end % myWriteTimeData


function tt = fieldTimes(fieldTokens, durationTimesteps)

timesteps = durationTimesteps(1):durationTimesteps(2);

tt = cell(numel(fieldTokens), 1);
for ff = 1:numel(fieldTokens)
    offset = t6.xml.fieldOffset(fieldTokens{ff});
    
    tt{ff} = (offset(4) + timesteps)*t6.simulation().Dt;
end

end % fieldTimes



function myWriteFunction(dataFileName, fieldFunction, fieldTokens, ...
    durationTimesteps)
precisionString = t6.TrogdorSimulation.instance().Precision;
fh = fopen(dataFileName, 'w');

if ~iscell(fieldFunction)
    fieldFunction = {fieldFunction};
end
%validateattributes(fieldFunction, {'cell'}, {}, 'myWriteFunction', 'fieldFunction');

numT = durationTimesteps(2) - durationTimesteps(1) + 1;

tt = fieldTimes(fieldTokens, durationTimesteps);

try
for nn = 1:numT
    for ff = 1:numel(fieldTokens)
        fwrite(fh, fieldFunction{ff}(tt{ff}(nn)), precisionString);
    end
end
catch
    error('Could not write TFSF source data file.');
end

fclose(fh);
end % myWriteFunction























