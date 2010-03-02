function fftDat = outputFFT(fileName)
%outputFFT Return the FFT of FDTD output data in a file.
%	fftDat = outputFFT(filePrefix) will have the same effect as
%	opening the fdtd output file and calling fft on the data
%	manually.  This is a nice shortcut, presuming that all the 
%	data can be fit in memory at once.
import t5.*
file = OutputFile(fileName);
data = file.read;

if iscell(data)
    error('This function does not work on multi-region outputs yet.');
end

fftDat = fft(data, [], ndims(data));