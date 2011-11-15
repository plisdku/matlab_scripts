function [x y z] = outputPositions(fileName)
%outputPositions Return the x, y and z coordinates of Trogdor output data
%   [x y z] = outputPositions('out') will return the positions at which
%   fields are measured in the output file called "out".  If there are
%   multiple regions and fields in the file, x, y and z will be cell arrays
%   indexed as x{region, field}.
%
%   Example: 'out' contains the hz field in a rectangle in 2d.
%
%   [x y] = outputPositions('out');
%   hz = spectrum('out', 'Frequency', 2*pi*3e8/500e-9); % extract 500 nm
%   figure;
%   imagesc(x, y, abs(hz)');
%   axis xy;
%
%   This is a convenience function wrapping OutputFile.positions().

fi = t6.OutputFile(fileName);
[x y z] = fi.positions();