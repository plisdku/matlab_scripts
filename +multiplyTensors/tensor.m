function A = tensor(szA, varargin)
% tensor  Index-value constructor for tensors
%
% tensor(sz, i1, i2, ..., iN, vals) creates an array representation of a tensor
% of dimension sz where vals(1) is in element (i1(1), i2(1), ..., iN(1)), etc.
% It is analogous to Matlab's sparse() function.

indices = cell(1, numel(szA));

if numel(varargin) ~= numel(indices) + 1
    error(['For %i-dimensional array, please provide %i index arrays ',...
        'and one value array'], numel(szA), numel(szA));
end
[indices{:}] = deal(varargin{1:end-1});
vals = varargin{end};

ii = sub2ind(szA, indices{:});

A = zeros([szA 1]);
A(ii) = vals;