function UU = applyCoordinates(U, whichCoords, xs, ys, zs, ts)
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
    % old way: just this one line!
    UU{ff} = U(cell2list(fieldCoords(ff), whichCoords));
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
cell2list = @(A, ii) A{ii};  % freaking Matlab... you suck.

args = fieldCoords(1); % any field would do here, 1 is arbitrary

UU = U(args{whichCoords});






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






