function B = multilinearInterp(varargin)
%multilinearInterp Faster, memory-efficient interpn for regular cartesian grids
%   v = gridInterp(X1,X2,X3,...,A,Y) interpolates to find v, the
%   values of the underlying N-D function A at all points Y(1,:), Y(2,:)...
%   Arrays X1, X2, etc. are 1D arrays of coordinates at which the data in A
%   are given.  Y is a [N Ndims] array of N point coordinates at which to
%   look up A.
%
%   If Xn for some n are empty arrays, then no interpolation
%   will be carried out in those dimensions.  This provides a way to select
%   which dimensions to interpolate.

numdim = nargin - 2;
%numdim = floor(nargin/2);

x = cell(1, numdim);

x0 = zeros(1,numdim);
dx = x0;
interpDims = x0;
keepDims = x0;

nInterp = 0;
nKeep = 0;
for dd = 1:numdim
    xdd = reshape(varargin{dd}, [], 1);
    
    x{dd} = xdd;
    
    if ~isempty(xdd)
        nInterp = nInterp+1;
        dx(nInterp) = xdd(2) - xdd(1);
        x0(nInterp) = xdd(1);
        interpDims(nInterp) = dd;
    else
        nKeep = nKeep+1;
        keepDims(nKeep) = dd;
    end
end

A = varargin{numdim+1};
y = varargin{numdim+2};

numSamples = size(y,1);
assert(size(y, 2) == nInterp);

%% Try using griddedInterpolant.
%{
% Convert distances to subscripts.  (Divide by dx.)

distance = max(bsxfun(@plus, y, -x0(1:nInterp)), 0);
subscripts = bsxfun(@times, distance, 1./dx(1:nInterp)) + 1;

%yCell = num2cell(y, 1);
subsCell = num2cell(subscripts,1);

if nKeep == 0 % interpolate over all coords
    gridMonster = griddedInterpolant(A);
    B = gridMonster(subsCell{:});
else
    szA = size(A);
    numKeepElements = prod(szA(keepDims(1:nKeep)));
    numInterpElements = prod(szA(interpDims(1:nInterp)));
    
    %{
    for kk = 1:numKeepElements
        subsA = indsA;
        subsA
        gridMonster.Values = squish(A(subsA{:}));
        B(kk,subsB{:}) = gridMonster(subsCell{:}
    %}
    % This was slow.
    
    A = permute(A, [keepDims(1:nKeep) interpDims(1:nInterp)]);
    A = reshape(A, [numKeepElements szA(interpDims(1:nInterp))]);
    
    szSlice = szA(interpDims(1:nInterp));
    
    % First make an interpolant that has the right dimensions.
    % Fill in the data later...
    gridMonster = griddedInterpolant(zeros(szSlice));
    
    B = zeros(numKeepElements, numSamples);
    
    for kk = 1:numKeepElements
        gridMonster.Values = reshape(A(kk, :), szSlice);
        B(kk,:) = reshape(gridMonster(subsCell{:}), 1, numSamples);
    end
    
    B = reshape(B, [szA(keepDims(1:nKeep)) numSamples]);
    
end
%}

%% My slow way... can't beat interpn... dangit.

%% Convert distances to subscripts.  (Divide by dx.)

distance = max(bsxfun(@plus, y, -x0(1:nInterp)), 0);
subscripts = bsxfun(@times, distance, 1./dx(1:nInterp)) + 1;

%% Method differs depending on how many dimensions to interpolate over!

if nKeep > 0
    B = subsInterp(A, subscripts, interpDims(1:nInterp));
else
    % I tried using a persistent variable here but it didn't seem to help.
    %persistent gridMonster;
    %if isempty(gridMonster) || ~isequal(size(gridMonster.Values), size(A))
        gridMonster = griddedInterpolant(A);
    %else
    %    gridMonster.Values = A;
    %end
    
    subsCell = num2cell(subscripts,1);
    B = gridMonster(subsCell{:});
end

%% Get x0 and dx for all directions.  (the OLD way.)

%szA = size(A);
%interpDims = find(cellfun(@(A) ~isempty(A), x));

%x0 = cellfun(@(A) A(1), x(interpDims));
%x1 = cellfun(@(A) A(2), x(interpDims));
%Nx = cellfun(@(A) length(A), x(interpDims));
%dx = x1 - x0;

end


function keepDims = otherDimensions(numDims, interpDims)

    keepDims = true(1, numDims);
    keepDims(interpDims) = false;
    keepDims = find(keepDims);
    
end