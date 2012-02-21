function addCurrentSource(varargin)
%addCurrentSource Add a source of electric and/or magnetic current density to the grid.
%   addCurrentSource('Field', 'my', 'YeeCells', [0 0 0 0 100 0], 'TimeData', ...
%       sin(1:numT)) will drive the grid with a sinusoidal magnetic current
%       density my along a line from (0,0,0) to (0,100,0) on every
%       timestep.
%
%   Usage: addCurrentSource(named parameters)
%
%   Named parameters:
%       Field       The electromagnetic fields to source.  Valid fields:
%                   'mx', 'my', 'mz'    magnetic current densities
%                   'jx', 'jy', 'jz'    electric current densities
%                   Any combination of fields may be used, but they should be
%                   specified in order: mx before my, electric before
%                   magnetic.
%                   (required)
%       YeeCells    The region of the grid in which to add electromagnetic
%                   current; [x0 y0 z0 x1 y1 z1] will source all cells (x, y, z)
%                   where x0 <= x <= x1, y0 <= y <= y1, z0 <= z <= z1.  Multiple
%                   rows may be used to source in multiple regions.
%                   (YeeCells or Bounds required)
%       Bounds      The region of the simulation space in which to add
%                   electromagnetic current; [x0 y0 z0 x1 y1 z1] in real units
%                   will be used to choose the YeeCells [m0 n0 p0 m1 n1 p1] in
%                   which to source fields, suitably for the grid resolution
%                   (YeeCells or Bounds required)
%       Duration    The range of timesteps on which to source currents; [t0 t1]
%                   will source on timesteps t such that t0 <= t <= t1.  Using
%                   multiple rows specifies multiple ranges of timesteps.
%                   (default: all timesteps)
%       FieldFunction   A function of space and time, e.g.
%                   @(x, y, z, t) sin(x).*cos(y).*tan(z).*exp(-t).
%                   If more than one current component is specified, then this
%                   argument must be a cell array with one function per field
%                   component, e.g. {jxFunc, jyFunc, jzFunc}.
%                   (FieldFunction, TimeData, or SpaceTimeData required)
%       TimeData    An array of size [nFields nTimesteps].  If the Duration
%                   is specified as [0 10] then TimeData needs 11 columns, one
%                   for each sourced timestep.
%                   (FieldFunction, TimeData or SpaceTimeFile required)
%       SpaceTimeData An array of size [nFields nTimesteps*nCells].
%       SpaceTimeFile   Name of a file that will contain the current for each
%                   cell at each timestep.  When using SpaceTimeFile, Trogdor
%                   will run once to generate a data request file named
%                   [SpaceTimeFile, '.m'] containing an ordered list of points
%                   and a set of timesteps at which to provide electromagnetic
%                   current values.  The space-time data file must be provided
%                   according to this data request before the simulation may
%                   run.  (SpaceTimeFile or TimeData required.)
%
%   Example:
%
%   addCurrentSource('Field', 'jz', 'YeeCells', [50 50 50 50 50 50; ...
%       60 60 60 60 60 60], 'Duration', [0 10; 100 110], 'SpaceTimeFile', 
%       'currentInput');
%
%   Specify a source of jz electric current in cells (50,50,50) and (60,60,60).
%   The first time Trogdor is run it will generate a data request file called
%   'currentInput.m' containing the ordered list of fields, cells and timesteps
%   at which data is required.  After this data is written in the correct order
%   to a binary file called 'currentInput', Trogdor can be run again to
%   completion.

import t6.*

grid = t6.TrogdorSimulation.instance().currentGrid();

X.Field = '';
X.YeeCells = [];
X.Bounds = [];
X.Duration = [];
X.FieldFunction = [];
X.TimeData = [];
X.SpaceTimeData = [];
%X.SpaceTimeFile = [];
X = parseargs(X, varargin{:});

t6.validateDataRequestParameters(X);
t6.validateYeeCellsAndBounds(X);

% Validate fields; should be a single string with some tokens in it
fieldTokens = tokenizeFields(X.Field, 'j m je mh');

% If we obtained Bounds and not YeeCells, set the YeeCells appropriately
if ~isempty(X.Bounds)
    X.YeeCells = boundsToYee(X.Bounds, fieldTokens);
end

% Validate duration
if isempty(X.Duration)
    X.Duration = [0, t6.TrogdorSimulation.instance().NumT-1];
elseif size(X.Duration, 2) ~= 2
    error('Duration must have two columns (first and last timestep).');
end

% Evaluate the field function at all the right places and times
if ~isempty(X.FieldFunction)
    
    if ~iscell(X.FieldFunction)
        X.FieldFunction = {X.FieldFunction};
    end
    
    %{
    if isempty(X.Bounds)
        X.SpaceTimeData = myCurrent_YeeCells(X.YeeCells(1,:),...
            X.Duration(1,:), fieldTokens, X.FieldFunction);
    else
        X.SpaceTimeData = myCurrent_Bounds(X.YeeCells(1,:), ...
            X.Bounds(1,:), X.Duration(1,:), fieldTokens, X.FieldFunction);
    end
    %}
    
end

obj = struct;
obj.type = 'CurrentSource';
obj.field = fieldTokens;
obj.yeeCells = X.YeeCells;
obj.bounds = X.Bounds;
obj.duration = X.Duration;
obj.fieldFunction = X.FieldFunction;
obj.timeData = X.TimeData;
obj.spaceTimeData = X.SpaceTimeData;

%obj.spaceTimeFile = X.SpaceTimeFile;

grid.CurrentSources = {grid.CurrentSources{:}, obj};



%{
% Find the correct (x,y,z,t) coordinates to evaluate the source functions
% at.
function src = myCurrent_YeeCells(yeeRegion, duration, fieldTokens,...
    fieldFunction)

src = zeros([yeeRegion(4:6)-yeeRegion(1:3)+1, numel(fieldTokens), ...
    duration(2) - duration(1,1) + 1]);

for ff = 1:numel(fieldTokens)
    offset = t6.xml.fieldOffset(fieldTokens{ff});

    yeeX = yeeRegion(1):yeeRegion(4);
    yeeY = yeeRegion(2):yeeRegion(5);
    yeeZ = yeeRegion(3):yeeRegion(6);
    timesteps = duration(1):duration(2);

    x = t6.grid().Origin(1) + offset(1) + yeeX*t6.simulation().Dxyz(1);
    y = t6.grid().Origin(2) + offset(2) + yeeY*t6.simulation().Dxyz(2);
    z = t6.grid().Origin(3) + offset(3) + yeeZ*t6.simulation().Dxyz(3);
    t = offset(4) + timesteps*t6.simulation().Dt;

    [xx yy zz tt] = ndgrid(x,y,z,t);
    
    src(:,:,:,ff,:) = fieldFunction{ff}(xx,yy,zz,tt);
    
    %src(:,:,:,ff,:) = reshape(X.FieldFunction{ff}(xx(:), yy(:), zz(:), tt(:)),...
    %    size(xx));
end

return



function src = myCurrent_Bounds(yeeRegion, bounds, duration, fieldTokens,...
    fieldFunction)

dxyz = t6.simulation().Dxyz;
dt = t6.simulation().Dt;

src = zeros([yeeRegion(4:6)-yeeRegion(1:3)+1, numel(fieldTokens), ...
    duration(2) - duration(1) + 1]);

for ff = 1:numel(fieldTokens)
    offset = t6.xml.fieldOffset(fieldTokens{ff}) .* [dxyz dt];
    
    support = t6.boundsToYee(bounds(1,:), fieldTokens{ff});
    
    % Yee cells over which field ff is nonzero.  This is a subset of
    % yeeRegion.
    supportYee{1} = support(1):support(4);
    supportYee{2} = support(2):support(5);
    supportYee{3} = support(3):support(6);
    
    % Physical points in space at which the current will be provided
    currCoords = cell(3,1);
    currCoords{1} = t6.grid().Origin(1) + offset(1) + supportYee{1}*dxyz(1);
    currCoords{2} = t6.grid().Origin(2) + offset(2) + supportYee{2}*dxyz(2);
    currCoords{3} = t6.grid().Origin(3) + offset(3) + supportYee{3}*dxyz(3);
    
    evalCoords = cell(3,1); % points in real space at which to evaluate function
    weights = cell(3,1); % interpolation weights
    for xyz = 1:3
        
        if t6.grid().numCells(xyz) == 1
            if numel(supportYee{xyz}) > 1
                error(['Current source bounds are too large in the ', ...
                    '%s direction'], char('w'+xyz));
            end
            
            evalCoords{xyz} = bounds(xyz);
            weights{xyz} = 1.0;
            
        elseif bounds(xyz) == bounds(xyz+3) % zero-width in this dim
            if numel(supportYee{xyz}) > 1
                % if the grid itself is not lower-dimensioned, then this is
                % really a delta-function distribution in one dimension.
                % Take the specified current to be a total current and
                % determine the appropriate density by dividing by dx.
                evalCoords{xyz} = [bounds(xyz), bounds(xyz)];
                weights{xyz} = (1 - abs(evalCoords{xyz}-currCoords{xyz})/dxyz(xyz))/dxyz(xyz);
                
            else % the grid is low-dimensioned.
                evalCoords{xyz} = bounds(xyz);
                weights{xyz} = 1.0/dxyz(xyz);
            end
        else
            % Physical bounds of cells centered at current samples
            cellLeft = max(bounds(xyz), currCoords{xyz} - dxyz(xyz));
            cellRight = min(bounds(xyz+3), currCoords{xyz} + dxyz(xyz));
            
            % This is a midpoint rule for a boxcar smoothing.
            evalCoords{xyz} = 0.5*(cellLeft + cellRight);
            weights{xyz} = 0.5*(cellRight-cellLeft)/dxyz(xyz);
            
            % This is a midpoint rule for a triangle smoothing.
            % Also, this is a crappy, crappy approximation, riddled with errors.
            % Whenever evalCoords == currCoords, w(x) == w(0) == 2, when it
            % really ought to be 1... (that is, 2*(cellRight-cellLeft)).
            %w = @(x) 1-abs(x);
            %
            %evalCoords{xyz} = 0.5*(cellLeft + cellRight);
            %weights{xyz} = w( (evalCoords{xyz} - currCoords{xyz})/dxyz(xyz) ) .* ...
            %    (cellRight-cellLeft)/dxyz(xyz);
        end
        
    end
    
    %fprintf('Field %i weights:\n', ff);
    %weights{1}
    %fprintf('Sum of weights: %2.4g\n', sum(weights{1}))
    %keyboard
    
    %fprintf('Source distrib:\n');
    %for xx = 1:numel(evalCoords{1})
    %    fprintf('Weight %2.2f at %2.2f\n', weights{1}(xx), ...
    %        currCoords{1}(xx));
    %end
    
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
    t = actualTimeOffset + timesteps*t6.simulation().Dt;
    
    [xx yy zz tt] = ndgrid(evalCoords{:},t);
    
    rawCurrent = fieldFunction{ff}(xx,yy,zz,tt);
    
    %figure(1000)
    %plot(rawCurrent(:));
    %pause
    
    scaledCurrent = bsxfun(@times, reshape(weights{1}, [], 1, 1), ...
        bsxfun(@times, reshape(weights{2}, 1, [], 1), ...
        bsxfun(@times, reshape(weights{3}, 1, 1, []), rawCurrent)));
    
    %keyboard
    
    % The scaled current may not be the full size of the source region.
    % Fill in the appropriate elements.  This amounts to finding an index 
    % offset.
    
    indicesX = 1 + supportYee{1} - yeeRegion(1);
    indicesY = 1 + supportYee{2} - yeeRegion(2);
    indicesZ = 1 + supportYee{3} - yeeRegion(3);
    
    src(indicesX, indicesY, indicesZ,ff,:) = scaledCurrent;
end

return

%}






