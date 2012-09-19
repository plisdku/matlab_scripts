function addHardSources(sim, grid, gridXML, doc)
global TROG_XML_COUNT___;

precisionString = sim.Precision;
directory = sim.directoryString;

for ss = 1:length(grid.HardSources)
    src = grid.HardSources{ss};
    
    elemXML = doc.createElement('HardSource');
    fieldstr = '';
    for ff = 1:length(src.field)
        fieldstr = [fieldstr, src.field{ff}, ' '];
    end
    elemXML.setAttribute('fields', fieldstr);
    
    % The spec file will include:
    % date
    % dxyz
    % dt
    % field ex (0.5, 0.5, 0.5) 0.0 etc.
    % unitVector0 (1, 0, 0) etc.
    % region [(a, b, c), (d, e, f)] stride (1, 1, 1) etc.
    % duration from A to B period C
    
    if length(src.timeData) ~= 0
        %fname = t6.xml.randomName('__hardsource_time_', '', 8);
        fname = [directory, sprintf('__hardsource_time_%i', TROG_XML_COUNT___.hardTime)];
        TROG_XML_COUNT___.hardTime = TROG_XML_COUNT___.hardTime + 1;
        elemXML.setAttribute('timeFile', fname);
        fh = fopen(fname, 'w');
        try
            count = fwrite(fh, src.timeData, precisionString);
        catch
            error('Could not write TFSF source data file.');
        end
        fclose(fh);
        t6.xml.writeSourceSpec(sim, src, 'AutoTimeFile', fname);
    end
    
    if length(src.maskData) ~= 0
        %fname = t6.xml.randomName('__hardsource_mask_', '', 8);
        fname = [directory, sprintf('__hardsource_mask_%i', TROG_XML_COUNT___.hardMask)];
        TROG_XML_COUNT___.hardMask = TROG_XML_COUNT___.hardMask + 1;
        elemXML.setAttribute('maskFile', fname);
        fh = fopen(fname, 'w');
        try
            % Write the mask in field order.  The mask is a cell array.
            % Inside each field component, write each region.
            % Inside each region, write in xyz order.
            for ff = 1:length(src.field)
            for mm = 1:length(src.maskData)    
                % dims: x y z fields.  Matlab writes in column order!!
                data = src.maskData{mm}(:,:,:,ff);
                %data = permute(src.maskData{mm}, [2 1 3 4]); % wrong
                count = fwrite(fh, data, precisionString);
            end
            end
        catch exception
            rethrow(exception)
            error('Could not write hard source mask data file.');
        end
        fclose(fh);
        t6.xml.writeSourceSpec(sim, src, 'AutoMaskFile', fname);
    end
    
    if length(src.spaceTimeData) ~= 0
        %fname = t6.xml.randomName('__hardsource_spacetime_', '', 8);
        fname = [directory, sprintf('__hardsource_spacetime_%i', ...
            TROG_XML_COUNT___.hardSpaceTime)];
        TROG_XML_COUNT___.hardSpaceTime = TROG_XML_COUNT___.hardSpaceTime + 1;
        elemXML.setAttribute('spaceTimeFile', fname);
        fh = fopen(fname, 'w');
        try
            % Write the per-timestep data in field order.  Regions make a cell
            % array.
            % Inside each field component, write each region.
            % Inside each region, write in xyz order.
            numT = size(src.spaceTimeData{1}, 5);
            % To speed up fwrite, we'll buffer all this data in-place.
            alldat = [];
            for tt = 1:numT
            for ff = 1:length(src.field)
            for mm = 1:length(src.spaceTimeData)    
                % dims: x y z fields.  Matlab writes in column order!!
                data = src.spaceTimeData{mm}(:,:,:,ff);
                %data = permute(src.maskData{mm}, [2 1 3]); % wrong
                alldat = [alldat, data(:)];
            end
            end
            end
            
            count = fwrite(fh, alldat, precisionString);
        catch
            error('Could not write hard source space-time data file.');
        end
        fclose(fh);
        t6.xml.writeSourceSpec(sim, src, 'AutoSpaceTimeFile', fname);
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
            sprintf('%i ', src.yeeCells(rr,:)));
        elemXML.appendChild(regionXML);
    end
    
    gridXML.appendChild(elemXML);
end