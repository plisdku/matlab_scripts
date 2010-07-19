function validateSourceParameters(X)
% Validate the various data parameters to addHardSource and addSoftSource.
% This function will:
% - make sure Yee cells are Nx6
% - make sure Duration is Mx2
% - make sure TimeData is as long as numT (total from Durations with period 1)
%   and has one value per input field
% - make sure MaskData has dimensions [x y z nFields] or, for multiple
%   regions, is a cell array where each cell has dimensions [x y z nFields] for
%   its region
% - make sure SpaceTimeData has dimensions [x y z nFields numT] or is right-
%   sized cell array as for MaskData
% Furthermore this makes sure that the pair (TimeData, MaskData) does not occur
% along with SpaceTimeData.
%
% As soon as a problem is found this function will throw an error.
% It has no return value.
grid = t6.TrogdorSimulation.instance().currentGrid();

% Validate fields; should be a single string with some tokens in it
fieldTokens = {};
remainder = X.Field;
while ~strcmp(remainder, '')
    [token, remainder] = strtok(remainder);
    if ~strcmp(token, '')
        fieldTokens = {fieldTokens{:}, token};
    end
end

if size(X.YeeCells, 2) ~= 6
    error('YeeCells must have six columns.');
end

if length(X.Duration) == 0
    X.Duration = [0, t6.TrogdorSimulation.instance().NumT-1];
elseif size(X.Duration, 2) ~= 2
    error('Duration must have two columns (first and last timestep).');
end

numT = sum(X.Duration(:,2) - X.Duration(:,1) + 1);

% 1/3  Validate time data.
if length(X.TimeData) ~= 0
    if length(X.TimeData) ~= numT
        error('TimeData must have the same length as the Duration');
    elseif size(X.TimeData, 1) ~= length(fieldTokens)
        error('TimeData must have size [nFields, timesteps]');
    end
end

% 2/3  Validate mask data.
if length(X.MaskData) ~= 0
    if size(X.YeeCells, 1) == 1
        % One mask region only; dims [x y z field]
        dim = X.YeeCells(4:6) - X.YeeCells(1:3) + 1;
        if size(X.MaskData,1) ~= dim(1) || size(X.MaskData,2) ~= dim(2) ...
            || size(X.MaskData,3) ~= dim(3)
            error('MaskData must have same size as YeeCells.');
        elseif size(X.MaskData, 4) ~= length(fieldTokens)
            error('Must provide MaskData for every field.');
        end
    else
        % Many mask regions; must pack into cell array
        if ~iscell(X.MaskData)
            error('MaskData for multiple-rectangle sources must be cell array.');
        end
        if length(X.MaskData) ~= size(X.YeeCells, 1)
            error('Provide MaskData for each rectangle.');
        end
        
        for mm = 1:length(X.MaskData)
            dim = X.YeeCells(mm, 4:6) - X.YeeCells(mm, 1:3) + 1;
            if size(X.MaskData{mm},1) ~= dim(1) || ...
                size(X.MaskData{mm},2) ~= dim(2) || ...
                size(X.MaskData{mm},3) ~= dim(3)
                error('MaskData #%i must have same size as YeeCells.', mm);
            elseif size(X.MaskData{mm}, 4) ~= length(fieldTokens)
                error('Must provide MaskData for every field.');
            end
        end
    end
end

% 3/3  Validate space-time data
if length(X.SpaceTimeData) ~= 0
    if size(X.YeeCells, 1) == 1
        dim = X.YeeCells(4:6) - X.YeeCells(1:3) + 1;
%        size(X.SpaceTimeData)
        if size(X.SpaceTimeData) ~= [dim, length(fieldTokens), numT]
            error(['Dimensions of SpaceTimeData do not match ', ...
                '[YeeCells, nFields, nTimesteps].']);
        end
    else
        % Many yee cell regions; must pack into cell array
        if ~iscell(X.SpaceTimeData)
            error(['SpaceTimeData for multiple-rectangle sources ',...
                'must be cell array.']);
        end
        if length(X.SpaceTimeData) ~= size(X.YeeCells, 1)
            error('Provide SpaceTimeData for each rectangle.');
        end
        
        for mm = 1:length(X.SpaceTimeData)
            dim = X.YeeCells(mm, 4:6) - X.YeeCells(mm, 1:3) + 1;
            if size(X.SpaceTimeData{mm},1) ~= dim(1) || ...
                size(X.SpaceTimeData{mm},2) ~= dim(2) || ...
                size(X.SpaceTimeData{mm},3) ~= dim(3)
                error('SpaceTimeData #%i must have same size as YeeCells.', mm);
            elseif size(X.SpaceTimeData{mm}, 4) ~= length(fieldTokens)
                error('Must provide SpaceTimeData for every field.');
            elseif size(X.SpaceTimeData{mm}, 5) ~= numT
                error(['SpaceTimeData for each rectangle must be long ',...
                    'enough for all timesteps.']);
            end
        end
    end
end

% Last: cross-validate.  Some of the options are exclusive:
if length(X.TimeData) ~= 0
    if length(X.SpaceTimeData) ~= 0
        error('Only one of TimeData and SpaceTimeData can be specified.');
    end
end

if length(X.SpaceTimeData) ~= 0
    if length(X.MaskData) ~= 0
        error('Space-Time sources cannot use MaskData.');
    end
end



