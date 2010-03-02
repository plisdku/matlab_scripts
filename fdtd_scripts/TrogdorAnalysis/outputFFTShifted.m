function fftDat = outputFFTShifted(filePrefix)
%
%   version 4.5
%   July 29, 2008


[fid, dim] = openOutputFile(filePrefix);
data = getOutputFrames(fid, dim);
fclose(fid);

datasize = size(data);

fftDat = zeros(size(data));

if (dim(1) == 1 & dim(2) == 1 & dim(3) == 1)
    %disp('0D data');
    %   single-point data in sequence
    fftDat = fftshift(fft(data));
    
elseif (dim(1)+dim(2) == 2 | dim(2)+dim(3) == 2 | dim(1)+dim(3) == 2)
    %disp('1D data');
    %   1D spatial data in sequence
    % Shift left to get the time index to the first dimension.
    % Then take the FFT and shift back to put frequency in the time
    % dimension again.
    fftDat = shiftdim(fftshift(fft(shiftdim(data, 3))), 1);
    
elseif (dim(1) == 1 | dim(2) == 1 | dim(3) == 1)
    %disp('2D data');
    %   2D spatial data in sequence
    %   Shift left to get the time index to the first dimension for fft...
    %   same as above.
    fftDat = shiftdim(fftshift(fft(shiftdim(data, 3))), 1);
    
else
    %disp('3D data');
    %   3D spatial data in sequence.
    % see above.
    fftDat = shiftdim(fftshift(fft(shiftdim(data, 3))), 1);
    
end

% Singleton dimensions have been stripped by the shiftdim calls, so let's
% put them back in.
fftDat = reshape(fftDat, datasize);




