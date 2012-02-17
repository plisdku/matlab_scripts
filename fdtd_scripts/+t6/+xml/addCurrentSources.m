function addCurrentSources(grid, gridXML, doc, origin)
global TROG_XML_COUNT___;

directory = t6.TrogdorSimulation.instance().directoryString;

for ss = 1:length(grid.CurrentSources)
    src = grid.CurrentSources{ss};
    
    elemXML = doc.createElement('CurrentSource');
    fieldstr = '';
    for ff = 1:length(src.field)
        fieldstr = [fieldstr, src.field{ff}, ' '];
    end
    elemXML.setAttribute('fields', fieldstr);
    
    if ~isempty(src.timeData)
        fname = [directory, sprintf('__currentsource_time_%i', TROG_XML_COUNT___.currentTime)];
        TROG_XML_COUNT___.currentTime = TROG_XML_COUNT___.currentTime + 1;
        elemXML.setAttribute('timeFile', fname);
        
        myWriteTimeFile(src.timeData);
        
        t6.xml.writeSourceSpec(src, 'AutoTimeFile', fname);
    end
    
    if ~isempty(src.spaceTimeData)
        fname = [directory, sprintf('__currentsource_spacetime_%i',...
            TROG_XML_COUNT___.currentSpaceTime)];
        TROG_XML_COUNT___.currentSpaceTime = TROG_XML_COUNT___.currentSpaceTime + 1;
        elemXML.setAttribute('spaceTimeFile', fname);
        
        myWriteSpaceTimeData(fname, src.spaceTimeData);
        
        t6.xml.writeSourceSpec(src, 'AutoSpaceTimeFile', fname);
    end
    
    if ~isempty(src.fieldFunction)
        fname = [directory, sprintf('__currentsource_spacetime_%i',...
            TROG_XML_COUNT___.currentSpaceTime)];
        TROG_XML_COUNT___.currentSpaceTime = TROG_XML_COUNT___.currentSpaceTime + 1;
        elemXML.setAttribute('spaceTimeFile', fname);
        
        if isempty(src.bounds)
            myWriteCurrent_YeeCells(fname, src.yeeCells, src.field, ...
                src.duration, src.fieldFunction);
        else
            myWriteCurrent_Bounds(fname, src.yeeCells, src.bounds, src.field, ...
                src.duration, src.fieldFunction);
        end
        
        t6.xml.writeSourceSpec(src, 'AutoSpaceTimeFile', fname);
    end
    
    if isfield(src, 'spaceTimeFile')
        error('SpaceTimeFile should be gone now');
    end
    %if ~isempty(src.spaceTimeFile)
    %    elemXML.setAttribute('spaceTimeFile', src.spaceTimeFile);
    %    t6.xml.writeSourceSpec(src);
    %end
    
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


function myWriteTimeFile(fname, timeData)
precisionString = t6.TrogdorSimulation.instance().Precision;
fh = fopen(fname, 'w');
try
    count = fwrite(fh, timeData, precisionString);
catch
    error('Could not write current source data file.');
end
fclose(fh);

return


function myWriteSpaceTimeData(fname, spaceTimeData)
precisionString = t6.TrogdorSimulation.instance().Precision;
fh = fopen(fname, 'w');
try
    count = fwrite(fh, spaceTimeData, precisionString);
catch
    error('Could not write current source data file.');
end
fclose(fh);

return


% Find the correct (x,y,z,t) coordinates to evaluate the source functions
% at.
function myWriteCurrent_YeeCells(fname, yeeRegion, fieldTokens, duration, fieldFunction)
%function src = myCurrent_YeeCells(yeeRegion, duration, fieldTokens,...
%    fieldFunction)
precisionString = t6.TrogdorSimulation.instance().Precision;
fh = fopen(fname, 'w');

for ff = 1:numel(fieldTokens)
    offset = t6.xml.fieldOffset(fieldTokens{ff});

    yeeX = yeeRegion(1):yeeRegion(4);
    yeeY = yeeRegion(2):yeeRegion(5);
    yeeZ = yeeRegion(3):yeeRegion(6);
    timesteps = duration(1):duration(2);

    fieldArgs(ff).x = t6.grid().Origin(1) + offset(1) + yeeX*t6.sim().Dxyz(1);
    fieldArgs(ff).y = t6.grid().Origin(2) + offset(2) + yeeY*t6.sim().Dxyz(2);
    fieldArgs(ff).z = t6.grid().Origin(3) + offset(3) + yeeZ*t6.sim().Dxyz(3);
    fieldArgs(ff).t = offset(4) + timesteps*t6.sim().Dt;
    
    [fieldArgs(ff).xx fieldArgs(ff).yy fieldArgs(ff).zz] = ndgrid(x,y,z);
end

numT = duration(2)-duration(1)+1;
for nn = 1:numT
    for ff = 1:numel(fieldTokens)
        fwrite(fh, fieldFunction{ff}(...
            fieldArgs(ff).xx, fieldArgs(ff).yy, fieldArgs(ff).zz, ...
            fieldArgs(ff).t(nn)), ...
            precisionString);
    end
end

return








function myWriteCurrent_Bounds(fname, yeeRegion, bounds, fieldTokens, duration, fieldFunction)
%function src = myCurrent_Bounds(yeeRegion, bounds, duration, fieldTokens,...
%    fieldFunction)
precisionString = t6.TrogdorSimulation.instance().Precision;
fh = fopen(fname, 'w');

dxyz = t6.sim().Dxyz;
dt = t6.sim().Dt;

%src = zeros([yeeRegion(4:6)-yeeRegion(1:3)+1, numel(fieldTokens), ...
%    duration(2) - duration(1) + 1]);

for ff = 1:numel(fieldTokens)
    offset = t6.xml.fieldOffset(fieldTokens{ff}) .* [dxyz dt];
    
    support = t6.boundsToYee(bounds(1,:), fieldTokens{ff});
    
    % Yee cells over which field ff is nonzero.  This is a subset of
    % yeeRegion.
    supportYee{1} = support(1):support(4);
    supportYee{2} = support(2):support(5);
    supportYee{3} = support(3):support(6);
    
    % Indices into the source array for each timestep
    fieldArgs(ff).indicesX = 1 + supportYee{1} - yeeRegion(1);
    fieldArgs(ff).indicesY = 1 + supportYee{2} - yeeRegion(2);
    fieldArgs(ff).indicesZ = 1 + supportYee{3} - yeeRegion(3);
    
    % Physical points in space at which the current will be provided
    currCoords = cell(3,1);
    currCoords{1} = t6.grid().Origin(1) + offset(1) + supportYee{1}*dxyz(1);
    currCoords{2} = t6.grid().Origin(2) + offset(2) + supportYee{2}*dxyz(2);
    currCoords{3} = t6.grid().Origin(3) + offset(3) + supportYee{3}*dxyz(3);
    
    evalCoords = cell(3,1); % points in real space at which to evaluate function
    fieldArgs(ff).weights = cell(3,1); % interpolation weights
    for xyz = 1:3
        
        if t6.grid().numCells(xyz) == 1
            if numel(supportYee{xyz}) > 1
                error(['Current source bounds are too large in the ', ...
                    '%s direction'], char('w'+xyz));
            end
            
            evalCoords{xyz} = bounds(xyz);
            fieldArgs(ff).weights{xyz} = 1.0;
            
        elseif bounds(xyz) == bounds(xyz+3) % zero-width in this dim
            if numel(supportYee{xyz}) > 1
                % if the grid itself is not lower-dimensioned, then this is
                % really a delta-function distribution in one dimension.
                % Take the specified current to be a total current and
                % determine the appropriate density by dividing by dx.
                evalCoords{xyz} = [bounds(xyz), bounds(xyz)];
                fieldArgs(ff).weights{xyz} = (1 - abs(evalCoords{xyz}-currCoords{xyz})/dxyz(xyz))/dxyz(xyz);
                
            else % the grid is low-dimensioned.
                evalCoords{xyz} = bounds(xyz);
                fieldArgs(ff).weights{xyz} = 1.0/dxyz(xyz);
            end
        else
            % Physical bounds of cells centered at current samples
            cellLeft = max(bounds(xyz), currCoords{xyz} - dxyz(xyz));
            cellRight = min(bounds(xyz+3), currCoords{xyz} + dxyz(xyz));
            
            % This is a midpoint rule for a boxcar smoothing.
            evalCoords{xyz} = 0.5*(cellLeft + cellRight);
            fieldArgs(ff).weights{xyz} = 0.5*(cellRight-cellLeft)/dxyz(xyz);
        end
        
    end
    
    % Things get a little gross here.  At timestep 0,
    % D(1) - D(0) = J(0.5)
    % B(1.5) - B(0.5) = M(1.0)
    %
    % so I have to manually adjust the offset of the M field to move
    % forward one timestep.
    
    actualTimeOffset = offset(4);
    if offset(4) == 0
        actualTimeOffset = offset(4) + dt;
    end
    
    timesteps = duration(1):duration(2);
    fieldArgs(ff).t = actualTimeOffset + timesteps*t6.sim().Dt;
    
    [fieldArgs(ff).xx fieldArgs(ff).yy fieldArgs(ff).zz] = ndgrid(evalCoords{:});
end

% The chunkwise approach does not make this faster.
%frameSize = [yeeRegion(4:6)-yeeRegion(1:3)+1, numel(fieldTokens)];
%frameBytes = prod(frameSize)*8;
%chunkBytes = 1e4;
%[chunkStarts, chunkEnds, chunkLengths] = t6.OutputFile.chunkTimesteps(...
%    duration(1), duration(2), frameBytes, chunkBytes);

numT = duration(2)-duration(1)+1;
%for cc = 1:length(chunkStarts)
for nn = 1:numT
    src = zeros([yeeRegion(4:6)-yeeRegion(1:3)+1, numel(fieldTokens)]); %, chunkLengths(cc)]);
    for ff = 1:numel(fieldTokens)
        % Chunkwise: slow.
        %t = fieldArgs(ff).t( (chunkStarts(cc):chunkEnds(cc))-duration(1)+1);
        %xs = repmat(fieldArgs(ff).xx, [1 1 1 length(t)]);
        %ys = repmat(fieldArgs(ff).yy, [1 1 1 length(t)]);
        %zs = repmat(fieldArgs(ff).zz, [1 1 1 length(t)]);
        %ts = repmat(reshape(t, 1, 1, 1, []), [size(fieldArgs(ff).xx) 1]);
        %rawCurrent = fieldFunction{ff}(xs, ys, zs, ts);
        
        rawCurrent = fieldFunction{ff}(...
            fieldArgs(ff).xx,fieldArgs(ff).yy,fieldArgs(ff).zz,...
            fieldArgs(ff).t(nn));
        
        scaledCurrent = bsxfun(@times, reshape(fieldArgs(ff).weights{1}, [], 1, 1), ...
            bsxfun(@times, reshape(fieldArgs(ff).weights{2}, 1, [], 1), ...
            bsxfun(@times, reshape(fieldArgs(ff).weights{3}, 1, 1, []), rawCurrent)));
        
        % The scaled current may not be the full size of the source region.
        % Fill in the appropriate elements.  This amounts to finding an index 
        % offset.
        
        src(fieldArgs(ff).indicesX, fieldArgs(ff).indicesY, ...
            fieldArgs(ff).indicesZ,ff,:) = scaledCurrent;
    end
    
    fwrite(fh, src(:), precisionString);
end



fclose(fh);

return



