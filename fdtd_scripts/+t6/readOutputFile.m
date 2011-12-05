function data = readOutputFile(fileName, varargin)
%readOutputFile Grab all the data from the named output file
%   alldat = readOutputFile(filename) is a shortcut for using OutputFile to
%   read the file.  (Provided for back-compatibility with Trogdor 4 scripts.)

X.Regions = 'Separate';
X = parseargs(X, varargin{:});

file = t6.OutputFile(fileName);
data = file.read('Regions', X.Regions);
