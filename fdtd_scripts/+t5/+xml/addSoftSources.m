function addSoftSources(grid, gridXML, doc, originTrogdor)
global TROG_XML_COUNT___;
originTwice = [originTrogdor originTrogdor];

for ss = 1:length(grid.SoftSources)
    src = grid.SoftSources{ss};
    
    elemXML = doc.createElement('AdditiveSource');
    fieldstr = '';
    for ff = 1:length(src.field)
        fieldstr = [fieldstr, src.field{ff}, ' '];
    end
    elemXML.setAttribute('fields', fieldstr);
    
    if length(src.timeData) ~= 0
        %fname = t5.xml.randomName('__softsource_time_', '', 8);
        fname = sprintf('__softsource_time_%i', TROG_XML_COUNT___.softTime);
        TROG_XML_COUNT___.softTime = TROG_XML_COUNT___.softTime + 1;
        elemXML.setAttribute('timeFile', fname);
        fh = fopen(fname, 'w');
        try
            count = fwrite(fh, src.timeData, 'float32');
        catch
            error('Could not write source data file.');
        end
        fclose(fh);
        t5.xml.writeSourceSpec(src, 'AutoTimeFile', fname);
    end
    
    if length(src.maskData) ~= 0
        %fname = t5.xml.randomName('__softsource_mask_', '', 8);
        fname = sprintf('__softsource_mask_%i', TROG_XML_COUNT___.softMask);
        TROG_XML_COUNT___.softMask = TROG_XML_COUNT___.softMask + 1;
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
                count = fwrite(fh, data, 'float32');
            end
            end
        catch exception
            rethrow(exception)
            error('Could not write soft source mask data file.');
        end
        fclose(fh);
        t5.xml.writeSourceSpec(src, 'AutoMaskFile', fname);
    end
    
    if length(src.spaceTimeData) ~= 0
        %fname = t5.xml.randomName('__softsource_spacetime_', '', 8);
        fname = sprintf('__softsource_spacetime_%i', ...
            TROG_XML_COUNT___.softSpaceTime);
        TROG_XML_COUNT___.softSpaceTime = TROG_XML_COUNT___.softSpaceTime + 1;
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
                data = src.spaceTimeData{mm}(:,:,:,ff,tt);
                %data = permute(src.maskData{mm}, [2 1 3]); % wrong!
                alldat = [alldat, data(:)];
            end
            end
            end
            count = fwrite(fh, alldat, 'float32');
        catch
            error('Could not write soft source space-time data file.');
        end
        fclose(fh);
        
        t5.xml.writeSourceSpec(src, 'AutoSpaceTimeFile', fname);
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