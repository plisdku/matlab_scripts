function addCurrentSource(src, sim, gridXML, doc)
global TROG_XML_COUNT___;

directory = sim.directoryString;

fname = [directory, sprintf('__current_%i', TROG_XML_COUNT___.current)];
TROG_XML_COUNT___.current = TROG_XML_COUNT___.current + 1;

elemXML = doc.createElement('CurrentSource');
fieldstr = '';
for ff = 1:length(src.field)
    fieldstr = [fieldstr, src.field{ff}, ' '];
end
elemXML.setAttribute('fields', fieldstr);

% durations and regions.
for dd = 1:size(src.timesteps, 1)
    durXML = doc.createElement('Duration');
    durXML.setAttribute('firstTimestep', num2str(src.timesteps(dd,1)));
    durXML.setAttribute('lastTimestep', num2str(src.timesteps(dd,2)));
    elemXML.appendChild(durXML);
end

for rr = 1:size(src.yeeCells, 1)
    regionXML = doc.createElement('Region');
    regionXML.setAttribute('yeeCells', ...
        sprintf('%i ', src.yeeCells(rr,:)));
    elemXML.appendChild(regionXML);
end

gridXML.appendChild(elemXML);

writeFile = true;
if exist(fname, 'file') && ~src.overwrite
    writeFile = false;
    fprintf('Not overwriting current source %s.\n', fname);
end

% Time data
if ~isempty(src.timeData)
    elemXML.setAttribute('timeFile', fname);
    
    if writeFile
        myWriteTimeFile(fname, src.timeData, sim.Precision);
    end
    
    t6.xml.writeSourceSpec(sim, src, 'AutoTimeFile', fname);
end

% Space-time data
if ~isempty(src.spaceTimeData)
    elemXML.setAttribute('spaceTimeFile', fname);
    
    if writeFile
        myWriteSpaceTimeData(fname, src.spaceTimeData, sim.Precision);
    end
    
    t6.xml.writeSourceSpec(sim, src, 'AutoSpaceTimeFile', fname);
end

% Field function
if ~isempty(src.fieldFunction)
    elemXML.setAttribute('spaceTimeFile', fname);
    
    if writeFile
        if isempty(src.bounds)
            writeFunctionCurrent_Yee(fname, src.yeeCells, src.field, ...
                src.timesteps, src.fieldFunction, sim.Precision);
        else
            writeFunctionCurrent_Bounds(sim, fname, src.yeeCells, src.bounds,...
                src.field, src.timesteps, src.fieldFunction, sim.Precision);
        end
    end
    
    t6.xml.writeSourceSpec(sim, src, 'AutoSpaceTimeFile', fname);
end

% Field functor
if ~isempty(src.fieldFunctor)
    elemXML.setAttribute('spaceTimeFile', fname);
    
    if writeFile
        if isempty(src.bounds)
            writeFunctorCurrent_Yee(fname, src.yeeCells, src.field, ...
                src.timesteps, src.fieldFunctor, sim.Precision);
        else
            writeFunctorCurrent_Bounds(sim, fname, src.yeeCells, src.bounds, src.field, ...
                src.timesteps, src.fieldFunctor, sim.Precision);
        end
    end
    
    %warning('Not writing current')
    
    t6.xml.writeSourceSpec(sim, src, 'AutoSpaceTimeFile', fname);
end

if isfield(src, 'spaceTimeFile')
    error('SpaceTimeFile should be gone now');
end



function myWriteTimeFile(fname, timeData, precisionString)

fh = fopen(fname, 'w');
try
    count = fwrite(fh, timeData, precisionString);
catch
    error('Could not write current source data file.');
end
fclose(fh);

return


function myWriteSpaceTimeData(fname, spaceTimeData, precisionString)
fh = fopen(fname, 'w');
try
    count = fwrite(fh, spaceTimeData, precisionString);
catch
    error('Could not write current source data file.');
end
fclose(fh);

return


function fieldArgs = yeeCellArguments(yeeRegion, fieldTokens, duration)

%grid = sim.CurrentGrid;

fieldArgs = struct('x', cell(size(fieldTokens)), ...
    'y', cell(size(fieldTokens)), ...
    'z', cell(size(fieldTokens)), ...
    't', cell(size(fieldTokens)));
for ff = 1:numel(fieldTokens)
    offset = t6.fieldOffset(fieldTokens{ff});

    yeeX = yeeRegion(1):yeeRegion(4);
    yeeY = yeeRegion(2):yeeRegion(5);
    yeeZ = yeeRegion(3):yeeRegion(6);
    timesteps = duration(1):duration(2);

    fieldArgs(ff).x = sim.Grid.Origin(1) + (offset(1) + yeeX)*sim.Dxyz(1);
    fieldArgs(ff).y = sim.Grid.Origin(2) + (offset(2) + yeeY)*sim.Dxyz(2);
    fieldArgs(ff).z = sim.Grid.Origin(3) + (offset(3) + yeeZ)*sim.Dxyz(3);
    fieldArgs(ff).t = (offset(4) + timesteps)*sim.Dt;
    
    [fieldArgs(ff).xx, fieldArgs(ff).yy, fieldArgs(ff).zz] = ndgrid(x,y,z);
end

return


function writeFunctionCurrent_Yee(fname, yeeRegion, fieldTokens, duration, fieldFunction, precisionString)

fh = fopen(fname, 'w');

fieldArgs = yeeCellArguments(yeeRegion, fieldTokens, duration);

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


function writeFunctorCurrent_Yee(fname, yeeRegion, fieldTokens, duration, fieldFunctor, precisionString)

fh = fopen(fname, 'w');

warning('This function has not been tested.')

fieldArgs = yeeCellArguments(yeeRegion, fieldTokens, duration);

fieldFunction = cell(size(fieldTokens));
for ff = 1:numel(fieldTokens)
    fieldFunction{ff} = fieldFunctor{ff}(...
        fieldArgs(ff).xx, fieldArgs(ff).yy, fieldArgs(ff).zz);
end

numT = duration(2)-duration(1)+1;
for nn = 1:numT
    for ff = 1:numel(fieldTokens)
        fwrite(fh, fieldFunction{ff}(fieldArgs(ff).t(nn)), precisionString);
    end
end

return




function fieldArgs = boundsArguments(sim, yeeRegion, bounds, fieldTokens, duration)

empty = cell(size(fieldTokens));
fieldArgs = struct('xx', empty, 'yy', empty, 'zz', empty, 't', empty, ...
    'indicesX', empty, 'indicesY', empty, 'indicesZ', empty, ...
    'weights', empty);

dxyz = sim.Dxyz;
dt = sim.Dt;

supportYee = cell(3,1);
for ff = 1:numel(fieldTokens)
    offset = t6.fieldOffset(fieldTokens{ff}) .* [dxyz dt];
    
    support = sim.boundsToYee(bounds(1,:), fieldTokens{ff});
    
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
    currCoords{1} = sim.Grid.Origin(1) + offset(1) + supportYee{1}*dxyz(1);
    currCoords{2} = sim.Grid.Origin(2) + offset(2) + supportYee{2}*dxyz(2);
    currCoords{3} = sim.Grid.Origin(3) + offset(3) + supportYee{3}*dxyz(3);
    
    evalCoords = cell(3,1); % points in real space at which to evaluate function
    fieldArgs(ff).weights = cell(3,1); % interpolation weights
    for xyz = 1:3
        
        if sim.Grid.numCells(xyz) == 1
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
                fieldArgs(ff).weights{xyz} = ...
                    (1 - abs(evalCoords{xyz}-currCoords{xyz})/dxyz(xyz)) ...
                    / dxyz(xyz);
                
            else % the grid is low-dimensioned.
                evalCoords{xyz} = bounds(xyz);
                fieldArgs(ff).weights{xyz} = 1.0/dxyz(xyz);
            end
        else
            % Physical bounds of cells centered at current samples.
            % This segment is the intersection of the current source extent
            % and the segment [x0-dx, x0+dx] centered around x0 =
            % currCoords{xyz}.
            cellCenter = currCoords{xyz};
            cellLeft = max(bounds(xyz), cellCenter - dxyz(xyz));
            cellRight = min(bounds(xyz+3), cellCenter + dxyz(xyz));
            
            % Evaluate the current at the center of the "cell".  We will,
            % however, assume that the current is constant.
            evalCoords{xyz} = 0.5*(cellLeft + cellRight);
            
            % The weight function comes from the triangular weight function
            %   w(x) = 1+x, -1 <= x <= 0
            %        = 1-x, 0 <= x <= 1
            % To perform the integral int_cell J(x) w(x) dx, we assume the
            % current J(x) = J(x_center) is constant:
            %   integral = J(x_center) * int_cell w(x) dx
            % and then figure out the integral of the triangle function
            % over the cell.  Do the integral piecewise in two parts:
            
            % Left-hand integral:
            xl = cellLeft;
            xr = min(cellRight, cellCenter);            
            yl = (xl - cellCenter + dxyz(xyz))/dxyz(xyz);
            yr = (xr - cellCenter + dxyz(xyz))/dxyz(xyz);
            wLeft = max( (xr-xl)*0.5.*(yl + yr),  0);
            
            xr = cellRight;
            xl = max(cellLeft, cellCenter);
            yl = (cellCenter + dxyz(xyz) - xl)/dxyz(xyz);
            yr = (cellCenter + dxyz(xyz) - xr)/dxyz(xyz);
            wRight = max( (xr-xl)*0.5.*(yl + yr),  0);
            
            fieldArgs(ff).weights{xyz} = (wLeft + wRight)/dxyz(xyz);
            
            %fieldArgs(ff).weights{xyz} = 0.5*(cellRight-cellLeft)/dxyz(xyz);
        end
        
    end
    
    % Things get a little gross here.  At timestep 0,
    % D(1) - D(0) = J(0.5)
    % B(1.5) - B(0.5) = M(1.0)
    %
    % so I have to manually adjust the offset of the M field to move
    % forward one timestep.
    
    %fieldArgs(1).weights{1}
    %fieldArgs(1).weights{2}
    
    actualTimeOffset = offset(4);
    if offset(4) == 0
        actualTimeOffset = offset(4) + dt;
    end
    
    timesteps = duration(1):duration(2);
    fieldArgs(ff).t = actualTimeOffset + timesteps*dt;
    
    [fieldArgs(ff).xx fieldArgs(ff).yy fieldArgs(ff).zz] = ndgrid(evalCoords{:});
end

return


function writeFunctionCurrent_Bounds(sim, fname, yeeRegion, bounds, fieldTokens, duration, fieldFunction, precisionString)

fh = fopen(fname, 'w');

fieldArgs = boundsArguments(sim, yeeRegion, bounds, fieldTokens, duration);

numT = duration(2)-duration(1)+1;

for nn = 1:numT
    src = zeros([yeeRegion(4:6)-yeeRegion(1:3)+1, numel(fieldTokens)]); %, chunkLengths(cc)]);
    for ff = 1:numel(fieldTokens)
        
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



function writeFunctorCurrent_Bounds(sim, fname, yeeRegion, bounds, ...
    fieldTokens, duration, fieldFunctor, precisionString)

fh = fopen(fname, 'w');

fieldArgs = boundsArguments(sim, yeeRegion, bounds, fieldTokens, duration);

fieldFunction = cell(size(fieldTokens));
for ff = 1:numel(fieldTokens)
    fieldFunction{ff} = fieldFunctor{ff}(...
        fieldArgs(ff).xx, fieldArgs(ff).yy, fieldArgs(ff).zz);
end

numT = duration(2)-duration(1)+1;

for nn = 1:numT
    src = zeros([yeeRegion(4:6)-yeeRegion(1:3)+1, numel(fieldTokens)]); %, chunkLengths(cc)]);
    for ff = 1:numel(fieldTokens)
        
        rawCurrent = fieldFunction{ff}(fieldArgs(ff).t(nn));
        
        scaledCurrent = bsxfun(@times, reshape(fieldArgs(ff).weights{1}, [], 1, 1), ...
            bsxfun(@times, reshape(fieldArgs(ff).weights{2}, 1, [], 1), ...
            bsxfun(@times, reshape(fieldArgs(ff).weights{3}, 1, 1, []), rawCurrent)));
        
        % The scaled current may not be the full size of the source region.
        % Fill in the appropriate elements.  This amounts to finding an index 
        % offset.
        
        src(fieldArgs(ff).indicesX, fieldArgs(ff).indicesY, ...
            fieldArgs(ff).indicesZ,ff,:) = scaledCurrent;
    end
    
    %fprintf('%i of %i\n', nn, numT);
    fwrite(fh, src(:), precisionString);
end

fclose(fh);

return


% This function is faster than the non-chunked version UNTIL the copy into
% a subarray of src.  That stage is preposterously slow.
function writeFunctorCurrent_Bounds_Chunked(sim, fname, yeeRegion, bounds, fieldTokens, duration, fieldFunctor, precisionString)

fh = fopen(fname, 'w');

fieldArgs = boundsArguments(sim, yeeRegion, bounds, fieldTokens, duration);

fieldFunction = cell(size(fieldTokens));
for ff = 1:numel(fieldTokens)
    fieldFunction{ff} = fieldFunctor{ff}(...
        fieldArgs(ff).xx, fieldArgs(ff).yy, fieldArgs(ff).zz);
end

numT = duration(2)-duration(1)+1;

timestepSize = [yeeRegion(4:6)-yeeRegion(1:3)+1, numel(fieldTokens)];
maxChunkSize = 1e5;
[chunkStarts, chunkEnds, chunkLengths] = t6.OutputFile.chunkTimesteps(...
    duration(1)+1, duration(2)+1, prod(timestepSize), maxChunkSize);
numChunks = numel(chunkStarts);
fprintf('Chunk size %i or %i\n', chunkLengths(1), chunkLengths(end));

for cc = 1:numChunks
    src = zeros([yeeRegion(4:6)-yeeRegion(1:3)+1, numel(fieldTokens), chunkLengths(cc)]);
    for ff = 1:numel(fieldTokens)
        
        rawCurrent = fieldFunction{ff}(fieldArgs(ff).t(chunkStarts(cc):chunkEnds(cc)));
        
        scaledCurrent = bsxfun(@times, reshape(fieldArgs(ff).weights{1}, [], 1, 1), ...
            bsxfun(@times, reshape(fieldArgs(ff).weights{2}, 1, [], 1), ...
            bsxfun(@times, reshape(fieldArgs(ff).weights{3}, 1, 1, []), rawCurrent)));
        
        % The scaled current may not be the full size of the source region.
        % Fill in the appropriate elements.  This amounts to finding an index 
        % offset.
        
        %numel(src(fieldArgs(ff).indicesX, fieldArgs(ff).indicesY, ...
        %    fieldArgs(ff).indicesZ, ff, :))
        %numel(scaledCurrent)
        
        src(fieldArgs(ff).indicesX, fieldArgs(ff).indicesY, ...
            fieldArgs(ff).indicesZ,ff,:) = scaledCurrent;
    end
    
    fwrite(fh, src(:), precisionString);
end

fclose(fh);

return




