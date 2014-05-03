function UU = applyCoords_Multiply(U, whichCoords, xs, ys, zs, ts)
% UU = applyCoords_Multiply(U, whichCoords, xs, ys, zs, ts)
%
% Evaluate U at the points in xs, ys, zs, ts.
% whichCoords selects arguments for U from xs, ys, zs ts.  For instance:
%  if U = U(y), then whichCoords = [2].
%  if U = U(x,y), then whichCoords = [1 2].
% UU will still be the appropriate size from xs, ys, zs, ts.

% How does the user signal that U is separately evaluated for each field?
% Well, if U returns four fields, that's a signal.
% If U depends on t, that's a signal too.

if ~isa(U, 'function_handle')
    UU = U;
    return;
end

validateArguments(U, whichCoords, xs, ys, zs, ts);

% We sort of have to guess about what to do here.
%
% We separately evaluate per field if:
%   U depends on a coordinate which depends on field (e.g. t, always)
%
% and other no other circumstances... dang!

% Figure it out!
areCells = [iscell(xs) iscell(ys) iscell(zs) false iscell(ts)];
useCoord = ismember(1:5, whichCoords);

if any(areCells & useCoord)
    UU = oneArrayPerField(U, whichCoords, xs, ys, zs, ts);
else
    UU = justOneArray(U, whichCoords, xs, ys, zs, ts);
end


function aa = cellOrNot(A, index)
% cellOrNot   Return ith element of cell array, else entire array.

if iscell(A)
    aa = A{index};
else
    aa = A;
end




function UU = oneArrayPerField(U, whichCoords, xs, ys, zs, ts)

% If U returns a sparse array then we need to make a cell array.
% If U returns a full array then we just concatenate arrays.
% So do the first array and see whether it's sparse or not.

if iscell(xs)
    numFields = numel(xs);
else
    numFields = numel(ts);
end


fieldCoords = @(ff) {cellOrNot(xs,ff), cellOrNot(ys, ff), ...
    cellOrNot(zs, ff), ff, cellOrNot(ts, ff)};
cell2list = @(A, ii) A{ii};  % freaking Matlab... you suck.

UU = cell(1, numFields);

for ff = 1:numFields
    
    xyzt = fieldCoords{ff};
    nxyzt = cellfun(@numel, xyzt);
    
    xyztBig = cell(numel(whichCoords), 1);
    [xyztBig{:}] = ndgrid(xyzt{whichCoords});
    
    Uvals = U(xyztBig{:});  % this doesn't interweave singular dims yet...
    
    szOut = [1 1 1 1 1];
    szOut(whichCoords) = nxyzt(whichCoords);
    
    UU{ff} = reshape(Uvals, szOut);
    
    %UU{ff} = U(cell2list(fieldCoords(ff), whichCoords));
end

% I think I won't actually want to concatenate these things... eh?

%if issparse(UU{1}) % Cell array is what we need for txca()
%    % ... we're fine, do nothing.
%else
%   UU = cat(4, UU{:}); 
%end




function UU = justOneArray(U, whichCoords, xs, ys, zs, ts)

fieldCoords = @(ff) {cellOrNot(xs,ff), cellOrNot(ys, ff), ...
    cellOrNot(zs, ff), ff, cellOrNot(ts, ff)};

%args = fieldCoords(1); % any field would do here, 1 is arbitrary

xyzft = fieldCoords(1); % cell array containing xs, ys, zs, fields, ts
nxyzft = cellfun(@numel, xyzft); % the sizes of xs, ys, zs, fields, ts

% This little kludge just allows me to use two incompatible behaviors to
% solve problems for my thesis.  :-D
doUseNDGRID = 0;

if doUseNDGRID
    xyzftBig = cell(numel(whichCoords),1); % container for ndgrid(xs,ys...)
    [xyzftBig{:}] = ndgrid(xyzft{whichCoords});
    Uvals = U(xyzftBig{:});
else
    Uvals = U(xyzft{whichCoords});
end

% The only complication now is that U may return several fields at once...
% we need to be field-agnostic!

szOut = [1 1 1 1 1];
szOut(whichCoords) = nxyzft(whichCoords);

% here we let the field dimension be whatever size it has to be.
UU = reshape(Uvals, szOut(1), szOut(2), szOut(3), [], szOut(5));

%UU = U(args{whichCoords});






%{
if iscell(fieldTimes)
    fprintf('Times are per-field.\n');
    
    if ismember(4, whichCoords)
        fprintf('Creating a cell array of filter data');
    end
    
else
    fprintf('Time is not per-field.\n');
    
    
    
end

%}

function validateArguments(U, whichCoords, xs, ys, zs, ts)

% 1.  Check consistency of use of cell arrays.
if ~isequal(iscell(xs), iscell(ys), iscell(zs))
    error('xs, ys and zs must either all be arrays or all be cell arrays');
end

if iscell(xs)
    if ~isequal(size(xs), size(ys), size(zs))
        error('As cell arrays, xs, ys and zs must all be the same size');
    end
    
    if iscell(ts) && ~isequal(size(xs), size(ts))
        error('If xs and ts are both cell arrays, they must have the same size');
    end
    
    sizeX = cellfun(@size, xs);
    sizeY = cellfun(@size, ys);
    sizeZ = cellfun(@size, zs);
    
    if ~isequal(sizeX{:})
        error('Each member of xs must be the same size array');
    end
    
    if ~isequal(sizeY{:})
        error('Each member of ys must be the same size array');
    end
    
    if ~isequal(sizeZ{:})
        error('Each member of zs must be the same size array');
    end
end

% 2.  Make sure whichCoords works.

if numel(whichCoords) > 5
    error('whichCoords has too many elements (%i > 4)', ...
        numel(whichCoords));
end

if any(whichCoords > 5) || any(whichCoords < 1)
    error('whichCoords may contain only the values 1, 2, 3 and 4');
end

if numel(unique(whichCoords)) ~= numel(whichCoords)
    error('whichCoords has repeated elements (%s)', num2str(whichCoords)');
end






