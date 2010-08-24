function harmonicComponent = outputHarmonic(filename, periods)
% outputHarmonic    Extract one frequency of a Fourier transform from data.
%   harmonicComponent = outputHarmonic(filename, period)
%   
%   filename    prefix for .dat and .txt FDTD output files
%
%   period      number of frames per oscillation at the frequency of
%               interest; Matlab's built-in FFT measures periods of
%               Inf, N, N/2, N/3, and so on where N is the number
%               of frames of the data.  Can be an array.
%   Provided for back-compatibility with Trogdor 4.

%filename = 'outYZE';
%periods = 1:2;
import t6.*

file = OutputFile(filename);
if length(file.Regions) > 1
    error('Cannot handle multi-region output files yet.');
end

if length(file.Fields) > 1
    dim = [file.Regions{1}.Size, length(file.Fields)];
else
    dim = file.Regions{1}.Size;
end

numFrames = file.numFramesAvailable;

endOfFile = 0;

harmonicComponent = zeros([dim, length(periods)]);

frameNumber = 1;

file.open;
if (length(periods) > 1)
    harmonicComponent = reshape(harmonicComponent, [], length(periods));
    
    harmonic = cell(1, length(periods));
    for pp = 1:length(periods)
        harmonic{pp} = zeros(prod(dim), 1);
    end
    
    while frameNumber <= numFrames
        phaseFactor = exp(-1i*2*pi*(frameNumber-1)./periods);
        framedat = file.readFrames('NumFrames', 1);
        for pp = 1:length(periods)
            harmonic{pp} = harmonic{pp} + phaseFactor(pp)*framedat(:);
        end
        
        frameNumber = frameNumber+1;
    end
    
    for pp = 1:length(periods)
        harmonicComponent(:,pp) = harmonic{pp};
    end
    
    harmonicComponent = reshape(harmonicComponent, [dim, length(periods)]);
else
    while frameNumber <= numFrames
        phaseFactor = exp(-1i*2*pi*(frameNumber-1)./periods);
        framedat = file.readFrames('NumFrames', 1);
        harmonicComponent = harmonicComponent + phaseFactor*framedat;
        frameNumber = frameNumber+1;
    end
end

file.close