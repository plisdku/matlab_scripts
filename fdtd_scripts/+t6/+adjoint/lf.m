function filter = lf(varargin)
% filter = lf(...)   Create a linear filter operation
%
% Usage:
%
%   filter = lf('Matrix', mat, 'dim', n);
%
%   filter = lf('Matrix', @(coords) mat, 'dim', n)
%
%   filter = lf('Multiply', array)
%
%   filter = lf('Multiply', @(x1,x2,x3) array, 'dim', [n m p])
%
%   filter = lf('Integrate', dim)
%
%   filter = lf('MatrixArray', cells, 'dim', n) % one matrix per field
%
% For any of these you may add the 'Field' argument to restrict operation
% to one or several fields (the fourth dimension of the data array, e.g.
% Ex and Hy only).
%

X.Matrix = [];
X.Multiply = [];
X.MatrixArray = [];
X.Integrate = [];
X.Dim = [];
X.Field = [];

X = parseargs(X, varargin{:});

if ~isempty(X.Matrix)
    filter = struct('Operation', 'Matrix', 'Data', X.Matrix, ...
        'dim', X.Dim, 'field', X.Field);
elseif ~isempty(X.MatrixArray)
    % mind the horrible syntax!
    filter = struct('Operation', 'MatrixArray', 'Data', {X.MatrixArray}, ...
        'dim', X.Dim, 'field', X.Field);
elseif ~isempty(X.Multiply)
    filter = struct('Operation', 'Pointwise', 'Data', X.Multiply, ...
        'dim', X.Dim, 'field', X.Field);
elseif ~isempty(X.Integrate)
     filter = struct('Operation', 'Matrix', ...
         'Data', @(x) 0.5*([0 diff(x)] + [diff(x) 0]), 'dim', X.Integrate);
end



