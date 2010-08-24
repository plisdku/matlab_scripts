function data = readOutputFile(fileName)
%readOutputFile Grab all the data from the named output file
%   alldat = readOutputFile(filename) is a shortcut for using OutputFile to
%   read the file.  (Provided for back-compatibility with Trogdor 4 scripts.)

import t5.*
file = OutputFile(fileName);
data = file.read;
