function data = readOutputFile(fileName)
%readOutputFile Grab all the data from the named output file
%   alldat = readOutputFile(filename) is a shortcut for using OutputFile to
%   read the file.  (Provided for back-compatibility with Trogdor 4 scripts.)
%
%   Usage:
%
%       fieldData = readOutputFile('electricFields');
%
%   fieldData will be an array of size [Nx Ny Nz Nf Nt] where
%       Nx is the number of output YeeCells in x
%       Ny is the number of output YeeCells in y
%       Nz is the number of output YeeCells in z
%       Nf is the number of fields (e.g. ex, ey and ez)
%       Nt is the number of saved timesteps
%
%   If there are multiple rows in the 'YeeCells' argument when addOutput is
%   called, readOutputFile will produce a cell array with one cell per row
%   of 'YeeCells'.  Each cell will contain an array of size [Nx Ny Nz Nf
%   Nt] corresponding to its YeeCell bounds.
%

import t5.*
file = OutputFile(fileName);
data = file.read;
