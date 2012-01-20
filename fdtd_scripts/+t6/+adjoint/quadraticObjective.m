function [f, Df, f_w, freqs] = quadraticObjective(data, varargin)
% quadraticObjective  Evaluate objective function and its gradient
%
% [f Df] = quadraticObjective(data) provides the sum of squares of the
% field values in data and the corresponding adjoint source Df.
%
% Named parameters:
%   Kernel
%   FilterT
%   FilterX
%   FilterY
%   FilterZ
%   Mask
%   Dt = 1
%
% The order of operations is first to do spatial and temporal filtering,
% then do operations on the fields (e.g. ExH and such).  The spatial
% operations are:
%
% mask * filterT * filterZ * filterY * filterX * fields
%
% Math background: let the objective function be x'*A*x.  Then the adjoint
% source should be (A' + A)*x.  I think.  I'm trying to do that.

import t6.*

X.Kernel = [];
X.FilterT = [];
X.FilterX = [];
X.FilterY = [];
X.FilterZ = [];
X.Mask = [];
X.Dt = 1;

X = parseargs(X, varargin{:});
numFields = size(data, 4);

if isempty(X.Kernel)
    X.Kernel = eye(numFields);
end

rightSize = @(A) (iscell(A) && numel(A) == numFields) || numel(A) == 0;

assert(rightSize(X.FilterT));
assert(rightSize(X.FilterX));
assert(rightSize(X.FilterY));
assert(rightSize(X.FilterZ));

for ff = 1:numFields
    currDat = data(:,:,:,ff,:);
    if ~isempty(X.FilterX) && ~isempty(X.FilterX{ff})
        currDat = multTensor(currDat, X.FilterX{ff}, 1);
    end
    if ~isempty(X.FilterY) && ~isempty(X.FilterY{ff})
        currDat = multTensor(currDat, X.FilterY{ff}, 2);
    end
    if ~isempty(X.FilterZ) && ~isempty(X.FilterZ{ff})
        currDat = multTensor(currDat, X.FilterZ{ff}, 3);
    end
    if ~isempty(X.FilterT) && ~isempty(X.FilterT{ff})
        currDat = multTensor(currDat, X.FilterT{ff}, 5);
    end
    if ~isempty(X.Mask) && ~isempty(X.Mask{ff})
        currDat = bsxfun(@times, currDat, X.Mask{ff});
    end
    
    if ff == 1
        filteredData = currDat;
    else
        filteredData(:,:,:,ff,:) = currDat;
    end
end

% Evaluate the objective function itself
% This should be f = x'*A*x.

quadraticAddend = conj(filteredData) .* multTensor(filteredData, X.Kernel, 4);

f = sum(quadraticAddend(:));

% Obtain Df
% it should be Df = (A' + A)*x, which we evaluate carefully.
% If K is the kernel, the other operations can be U.  Then
%
% A = U'*K*U
% A' = U'*K'*U
%
% so Df = U'*(K + K')*U*x.  The U*x part has been done already (it's data).

Df = multTensor(filteredData, X.Kernel + X.Kernel', 4);

for ff = 1:numFields
    currDf = Df(:,:,:,ff,:);
    if ~isempty(X.Mask) && ~isempty(X.Mask{ff})
        currDf = bsxfun(@times, currDf, X.Mask{ff});
    end
    if ~isempty(X.FilterT) && ~isempty(X.FilterT{ff})
        currDf = multTensor(currDf, X.FilterT{ff}', 5);
    end
    if ~isempty(X.FilterZ) && ~isempty(X.FilterZ{ff})
        currDf = multTensor(currDf, X.FilterZ{ff}', 3);
    end
    if ~isempty(X.FilterY) && ~isempty(X.FilterY{ff})
        currDf = multTensor(currDf, X.FilterY{ff}', 2);
    end
    if ~isempty(X.FilterX) && ~isempty(X.FilterX{ff})
        currDf = multTensor(currDf, X.FilterX{ff}', 1);
    end
    
    if ff == 1
        Df = currDf;
    else
        Df(:,:,:,ff,:) = currDf;
    end
end

% Obtain f_w if desired

[data_w, freqs] = t6.analysis.spectrum(data, 'Dt', X.Dt);

for ff = 1:numFields
    currData_w = data_w(:,:,:,ff,:);
    if ~isempty(X.FilterX) && ~isempty(X.FilterX{ff})
        currData_w = multTensor(currData_w, X.FilterX{ff}, 1);
    end
    if ~isempty(X.FilterY) && ~isempty(X.FilterY{ff})
        currData_w = multTensor(currData_w, X.FilterY{ff}, 2);
    end
    if ~isempty(X.FilterZ) && ~isempty(X.FilterZ{ff})
        currData_w = multTensor(currData_w, X.FilterZ{ff}, 3);
    end
    if ~isempty(X.FilterT) && ~isempty(X.FilterT{ff})
        currData_w = multTensor(currData_w, X.FilterT{ff}, 5);
    end
    if ~isempty(X.Mask) && ~isempty(X.Mask{ff})
        currData_w = bsxfun(@times, currData_w, X.Mask{ff});
    end
    
    if ff == 1
        filteredData_w = currData_w;
    else
        filteredData_w(:,:,:,ff,:) = currData_w;
    end
end

f_w = squeeze(sum(sum(sum(sum(conj(filteredData_w).*multTensor(filteredData_w, X.Kernel, ...
    4), 4), 3), 2), 1));

%% Multiply Tensors

function B = multTensor(T, A, dim)
% mult(T, A, dim) multiplies tensor T by matrix A over dimension dim
% (A*T)

%if isscalar(T) && dim ~= 1
%    warning('Tensor T is a scalar so multiplication is being shifted to first dimension');
%    dim = 1;
%end

if dim < 1 || dim > ndims(T)
    if size(A, 2) > 1
        error('dim must be valid dimension of T');
    else
        dim = 1;
    end
end

permutedDims = [dim, 1:dim-1, dim+1:ndims(T)];
sz = size(T);
szPermuted = sz(permutedDims);


tPerm = permute(T, permutedDims); % Put the desired dimension first
prodTensor = reshape(A*tPerm(:,:), [size(A,1), szPermuted(2:end)]);

B = ipermute(prodTensor, permutedDims);





