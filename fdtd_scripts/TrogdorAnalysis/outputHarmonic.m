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
%
%   version 4.5 revised
%   March 3, 2009


%filename = 'outYZE';
%periods = 1:2;

[fid, dim] = openOutputFile(filename);

endOfFile = 0;

harmonicComponent = zeros([dim, length(periods)]);
harmBuffer = zeros([length(periods), dim]);

frameNumber = 1;

if (length(periods) > 1)

    while (~endOfFile)

        phaseFactor = exp(-i*2*pi*(frameNumber-1)./periods);

        [framedat, count] = getOutputFrame(fid, dim);
        
        if (count ~= 0)
            for pp = 1:length(periods)
                phaseFrame = phaseFactor(pp)*framedat;
                harmBuffer(pp,:) = phaseFrame(:);
            end
            harmonicComponent = harmonicComponent + shiftdim(harmBuffer,1);
            frameNumber = frameNumber+1;
        else
            endOfFile = 1;
        end
    end
else
    
    while (~endOfFile)

        phaseFactor = exp(-i*2*pi*(frameNumber-1)./periods);

        [framedat, count] = getOutputFrame(fid, dim);

        if (count ~= 0)
            harmonicComponent = harmonicComponent + phaseFactor*framedat;
            %harmonicComponent = harmonicComponent + phaseBlock.*framedatBlock;
            frameNumber = frameNumber+1;
        else
            endOfFile = 1;
        end
    end
    
    
    
end
fclose(fid);
