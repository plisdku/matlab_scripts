function data = readOutputFile(fileName, varargin)
%readOutputFile Grab all the data from the named output file
%   alldat = readOutputFile(filename) is a shortcut for using OutputFile to
%   read the file.  (Provided for back-compatibility with Trogdor 4 scripts.)

import t6.*

X.Regions = 'Separate';
X = parseargs(X, varargin{:});

file = OutputFile(fileName);
data = file.read('Regions', X.Regions);
