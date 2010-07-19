function addGridReport(filename, varargin)
%addGridReport Tell Trogdor to write the material contents of the grid to a file
%   addGridReport('gridReportFile', 'YeeCells', [20 20 20 30 30 30]) will tell
%   Trogdor to save all the materials in the cube from (20,20,20) to (30,30,30)
%   to a file named 'gridReport'.  A spec file 'gridReport.txt' will also be
%   written, and the data may be opened as an OutputFile.
%
%   Usage: addGridReport(fileName, named parameters)
%
%   Named parameters:
%       YeeCells    The region of the grid to report;  [x0 y0 z0 x1 y1 z1] will
%                   save all cells (x, y, z) where x0 <= x <= x1, y0 <= y <= y1,
%                   z0 <= z <= z1.  (default: entire grid)
grid = t6.TrogdorSimulation.instance().currentGrid();

X.YeeCells = [];
X = parseargs(X, varargin{:});

if length(X.YeeCells) ~= 0 && ~t6.validateRect(X.YeeCells)
    error('Invalid YeeCells rectangle.');
end

obj = struct;
obj.type = 'GridReport';
obj.filename = filename;
obj.yeeCells = X.YeeCells;

grid.GridReports = {grid.GridReports{:}, obj};