function addCurrentSources(grid, gridXML, doc, originTrogdor)
global TROG_XML_COUNT___;
originTwice = [originTrogdor originTrogdor];

for ss = 1:length(grid.CurrentSources)
    src = grid.CurrentSources{ss};
    
    elemXML = doc.createElement('CurrentSource');
    fieldstr = '';
    for ff = 1:length(src.field)
        fieldstr = [fieldstr, src.field{ff}, ' '];
    end
    elemXML.setAttribute('fields', fieldstr);
    
    if length(src.timeData) ~= 0
        %fname = t6.xml.randomName('__currentsource_time_', '', 8);
        fname = sprintf('__currentsource_time_%i', TROG_XML_COUNT___.currentTime);
        TROG_XML_COUNT___.currentTime = TROG_XML_COUNT___.currentTime + 1;
        elemXML.setAttribute('timeFile', fname);
        fh = fopen(fname, 'w');
        try
            count = fwrite(fh, src.timeData, 'float32');
        catch
            error('Could not write current source data file.');
        end
        fclose(fh);
        t6.xml.writeSourceSpec(src, 'AutoTimeFile', fname);
    end
    
    if length(src.spaceTimeData) ~= 0
        fname = sprintf('__currentsource_spacetime_%i',...
            TROG_XML_COUNT___.currentSpaceTime);
        TROG_XML_COUNT___.currentSpaceTime = TROG_XML_COUNT___.currentSpaceTime + 1;
        elemXML.setAttribute('spaceTimeFile', fname);
        fh = fopen(fname, 'w');
        try
            count = fwrite(fh, src.spaceTimeData, 'float32');
        catch
            error('Could not write current source data file.');
        end
        fclose(fh);
        t6.xml.writeSourceSpec(src, 'AutoSpaceTimeFile', fname);
    end
    
    %if length(src.maskFile) ~= 0
    %    elemXML.setAttribute('maskFile', src.maskFile);
    %    t6.xml.writeSourceSpec(src);
    %end
    
    if length(src.spaceTimeFile) ~= 0
        elemXML.setAttribute('spaceTimeFile', src.spaceTimeFile);
        t6.xml.writeSourceSpec(src);
    end
    
    
    % durations and regions.
    for dd = 1:size(src.duration, 1)
        durXML = doc.createElement('Duration');
        durXML.setAttribute('firstTimestep', num2str(src.duration(dd,1)));
        durXML.setAttribute('lastTimestep', num2str(src.duration(dd,2)));
        elemXML.appendChild(durXML);
    end
    
    for rr = 1:size(src.yeeCells, 1)
        regionXML = doc.createElement('Region');
        regionXML.setAttribute('yeeCells', ...
            sprintf('%i ', src.yeeCells(rr,:) + originTwice));
        elemXML.appendChild(regionXML);
    end
    
    gridXML.appendChild(elemXML);
end