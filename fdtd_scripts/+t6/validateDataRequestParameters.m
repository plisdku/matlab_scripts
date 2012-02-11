function validateDataRequestParameters(X)
% Throw an error if:
%   X.YeeCells doesn't have six columns, [x0 y0 z0 x1 y1 z1]
%   X.Duration doesn't have two columns, [t0 t1]
%   X.TimeData doesn't have one value per field per timestep, size 
%       [numFields, numTimesteps]
%   X.SpaceTimeFile exists along with X.TimeFile or X.MaskFile
% Doesn't check if the fields are ordered right (E before H, K before J).

% Validate fields; should be a single string with some tokens in it
fieldTokens = {};
remainder = X.Field;
while ~strcmp(remainder, '')
    [token, remainder] = strtok(remainder);
    if ~strcmp(token, '')
        fieldTokens = {fieldTokens{:}, token};
    end
end

%if size(X.YeeCells, 2) ~= 6
%    error('YeeCells must have six columns.');
%end

if length(X.Duration) == 0
    X.Duration = [0, t6.TrogdorSimulation.instance().NumT-1];
elseif size(X.Duration, 2) ~= 2
    error('Duration must have two columns (first and last timestep).');
end

numT = sum(X.Duration(:,2) - X.Duration(:,1) + 1);

%  Validate time data.
if length(X.TimeData) ~= 0
    if length(X.TimeData) ~= numT
        error('TimeData must have the same length as the Duration');
    elseif size(X.TimeData, 1) ~= length(fieldTokens)
        error('TimeData must have size [nFields, timesteps] or [timesteps].');
    end
end

% Check for contradictory specifications
% space time file is deprecated
%{
if length(X.SpaceTimeFile) ~= 0
    if length(X.TimeData) ~= 0
        error('Cannot use both TimeData and SpaceTimeFile.');
    elseif length(X.MaskFile) ~= 0
        error('Cannot use both MaskFile and SpaceTimeFile.');
    end
end
%}