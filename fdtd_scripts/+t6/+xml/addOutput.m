function addOutput(output, sim, grid, gridXML, doc)

%directory = sim.outputDirectoryString;

elemXML = doc.createElement('FieldOutput');
elemXML.setAttribute('fields', output.fields);
if isfield(output, 'interpolationPoint')
    elemXML.setAttribute('interpolate', ...
        sprintf('%2.2f ', output.interpolationPoint));
end
%elemXML.setAttribute('file', [directory, output.filename]);
elemXML.setAttribute('file', output.filename);

% durations and regions.
for dd = 1:size(output.duration, 1)
    durXML = doc.createElement('Duration');
    durXML.setAttribute('firstTimestep', num2str(output.duration(dd,1)));
    durXML.setAttribute('lastTimestep', num2str(output.duration(dd,2)));
    durXML.setAttribute('period', num2str(output.period(dd)));
    elemXML.appendChild(durXML);
end

for rr = 1:size(output.yeeCells, 1)
    regionXML = doc.createElement('Region');
    regionXML.setAttribute('yeeCells', ...
        sprintf('%i ', output.yeeCells(rr,:)));
    regionXML.setAttribute('stride', sprintf('%i ', output.stride(rr,:)));
    
    if size(output.bounds, 1) == size(output.yeeCells, 1)
        regionXML.setAttribute('bounds', sprintf('%i ', output.bounds(rr,:)));
    end
    
    elemXML.appendChild(regionXML);
end

gridXML.appendChild(elemXML);

