function xyz = outputPositions(fileName, varargin)
%outputPositions Return the x, y and z coordinates of Trogdor output data
%   xyz = outputPositions('out') will return the positions at which
%   fields are measured in the output file called "out".
%
%   This is a convenience function wrapping OutputFile.positions().


%   Example: 'out' contains the hz field in a rectangle in 2d.
%   
%   [x y] = outputPositions('out');
%   hz = spectrum('out', 'Frequency', 2*pi*3e8/500e-9); % extract 500 nm
%   figure;
%   imagesc(x, y, abs(hz)');
%   axis xy;
%

fi = t6.OutputFile(fileName);
xyz = fi.positions(varargin{:});
