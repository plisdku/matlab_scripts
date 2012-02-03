function [f, freqs, resid, residFrac] = spectrum(fileName, varargin)
%spectrum Calculate discrete-time Fourier components of time-sampled data
%   THE NORMALIZATION OF THIS FUNCTION DIFFERS FROM FFT().  Read the note
%   at the end for details.
%
%   [F, freqs] = spectrum(f, 'Dt', dt) returns the discrete-time Fourier
%   transform of vector f (length N) and its corresponding frequencies
%   
%       freqs(k) = 2*pi*(k-1)/N*dt.
%   
%   The values of f are separated by a time interval of dt.  The Fourier
%   transform is always taken along the last dimension of f.
%
%   [F, freqs] = spectrum(f, 'Dt', dt, 'Frequency', frequencies) extracts
%   several frequency components of f without taking the full fft.  For
%   each given frequency,
%   
%                       N
%       F(k) = (1/N) * sum  f(n)*exp(-j*frequency*n*dt)
%                      n=1
%  
%   If the frequencies correspond to the frequencies of the fft, F(k) will
%   agree frequency-by-frequency with the result of spectrum(f, 'Dt', dt).
%
%   [F, freqs, residual, relativeResidual] = spectrum(f, 'Time', timeVals,
%   'SteadyStateFrequency', frequency) is used to obtain the amplitude and
%   phase of a signal that is known to have one dominant frequency
%   component: f(n) = F*exp(j*frequency*n*dt).  As few as two samples of f
%   may be used to obtain F.
%   
%   Example:
%       freq = 0.1;
%       t = 0:10;
%       f = cos(freq*t);
%       F = t5.analysis.spectrum(f(9:10), 'Time', t(9:10), ...
%           'SteadyStateFrequency', freq)
%   
%   Output:
%
%   F =
%
%       0.5000 - 0.0000i
%   
%   Providing more than two samples of f permits the calculation of a
%   residual, representing the part of f that is not well-represented by
%   the given sinusoid.  For convenience, the relative magnitude of the
%   residual relativeResidual = norm(residual)/norm(f) is given as well.
%   
%   [F, freqs] = spectrum(filename) returns the discrete-time Fourier
%   transform of the FDTD data from the given Trogdor output file.  The
%   time of each sample is obtained from the output file, and a phase
%   correction is applied to E and H as necessary to compensate for FDTD
%   interleaving.  (This function is equivalent to the old outputFFT() with
%   the additional phase correction and different amplitude normalization.
%   Read the note at the end of the help for details.)
%   
%   [F, freqs] = spectrum(filename, 'Frequency', frequencies) returns the
%   requested frequency components of the output file without using fft().
%   Large data files may be impossible to load into memory for fft().
%   Using this function the data are read one timestep at a time, so
%   frequency decompositions for large files may still be performed.  The
%   time of each sample is obtained from the output file, and a phase
%   correction is applied to E and H as necessary to compensate for FDTD
%   interleaving. (This function is equivalent to the old outputHarmonic()
%   with the additional phase correction and different amplitude
%   normalization.  Read the note at the end of the help for details.)
%
%   [F, freqs, residual, relativeResidual] = spectrum(filename,
%   'SteadyStateFrequency', frequency) is used to obtain the amplitude and
%   phase of a signal that is known to have one dominant frequency
%   component: f(n) = F*exp(j*frequency*n*dt).  As few as two samples of f
%   may be used to obtain F.  In practice this is useful for obtaining the
%   phasor components of an electromagnetic field in a volume.  After
%   running a simulation to steady-state, the last two (or N) timesteps of
%   the simulation may be saved to a file and used to extract the phasors.
%   The amount of disk space required is miniscule compared to saving every
%   timestep in the whole volume.  The time of each sample is obtained from
%   the output file, and a phase correction is applied to E and H as
%   necessary to compensate for FDTD interleaving.  The residual is also
%   provided (see above).
%
%   DEFINITION OF FREQUENCY DECOMPOSITION
%
%   A length N input vector f with time spacing dt can be decomposed into
%   its frequency components as
%
%               N
%       f(n) = sum  F(k)*exp(j*freq(k)*n*dt), 1 <= n <= N
%              k=1
%
%   The frequencies freq(k) = 2*pi*(k-1)/(N*dt).  This decomposition
%   differs from Matlab's ifft() by the omission of the prefactor (1/N).
%   By omitting this factor, the decomposition F(k) of a truncated signal
%   will have roughly the same amplitude regardless of the number of
%   samples in the signal and the sampling interval.  This facilitates
%   comparison of spectra of signals taken with different sampling
%   intervals or from FDTD simulations of different total durations.
%
%   DIMENSIONS:
%
%       fieldData = t5.analysis.spectrum('electricFields');
%
%   fieldData will be an array of size [Nx Ny Nz Nf Nfreq] where
%       Nx is the number of output YeeCells in x
%       Ny is the number of output YeeCells in y
%       Nz is the number of output YeeCells in z
%       Nf is the number of fields (e.g. ex, ey and ez)
%       Nfreq is the number of frequencies
%
%   If the output file contained multiple blocks of Yee cells, as would
%   happen when addOutput is called with a multi-row 'YeeCells' argument,
%   analysis.spectrum will produce a cell array with one cell per row
%   of 'YeeCells'.  Each cell will contain an array of size [Nx Ny Nz Nf
%   Nt] corresponding to its YeeCell bounds.

import t5.*;

X.Frequency = [];
X.SteadyStateFrequency = [];
X.Dt = [];
X.Time = [];
X = parseargs(X, varargin{:});

if ~isempty(X.Frequency) && ~isempty(X.SteadyStateFrequency)
    error('Frequency and SteadyStateFrequency may not both be specified');
end

if ischar(fileName)
    if nargout < 3
        [f, freqs] = doFileSpectrum(fileName, X);
    else
        [f, freqs, resid, residFrac] = doFileSpectrum(fileName, X);
    end
else
    if nargout < 3
        [f, freqs] = doDataSpectrum(fileName, X);
    else
        [f, freqs, resid, residFrac] = doDataSpectrum(fileName, X);
    end
end


function [f, freqs, resid, residFrac] = doDataSpectrum(data, X)
f = 0;
freqs = 0;

if iscell(data)
    error(['Cannot take spectrum of cell array/multi region data',...
        ' except from output file']);
end


if ndims(data) == 2
    if size(data, 2) > 1
        timeAxis = 2;
    else
        timeAxis = 1;
    end
else
    timeAxis = ndims(data);
end 
    
numT = size(data, timeAxis);

if ~isempty(X.Time) && ~isempty(X.Dt)
    error('Please provide either Dt or Times');
end

if ~isempty(X.Time)
    if numel(X.Time) ~= numT
        error('Number of time positions must equal last dimension of data');
    end
    
    X.Dt = X.Time(2) - X.Time(1);
else
    if isempty(X.Dt)
        error('Must provide either Dt or Time');
    end
    
    X.Time = (0:numT-1)*X.Dt;
end

if isempty(X.Frequency) && isempty(X.SteadyStateFrequency)
    
    freqs = 2*pi*(0:numT-1) / (numT * X.Dt);
    
    phaseFactor = exp(-1i*freqs*X.Time(1));
    phaseFactorSize = [ones(1, timeAxis-1), numel(phaseFactor)];
    
    f = bsxfun(@times, fft(data, [], timeAxis) / numT, ...
    	reshape(phaseFactor, phaseFactorSize));
    
elseif ~isempty(X.SteadyStateFrequency)
    
    if numel(X.SteadyStateFrequency) > 1
        error('Can only do one steady state frequency');
    end
    
    freqs = X.SteadyStateFrequency;
    
    calcResid = (nargout > 2);
    nd = ndims(data);
    timeFirstIndices = [timeAxis, 1:timeAxis-1, timeAxis+1:nd];
    timeFirstData = permute(data, timeFirstIndices);
    szTFD = size(timeFirstData);
    
    time = reshape(X.Time, 1, []);
    
    A = 2*[cos(freqs*time'), -sin(freqs*time')];
    pinvA = pinv(A);
    
    if calcResid
        timeFirstDecomp = pinvA*timeFirstData(:,:);
        
        timeFirstResid = reshape(A*timeFirstDecomp - timeFirstData(:,:),...
            szTFD);
        decomp = reshape(timeFirstDecomp, [2 szTFD(2:end)]);
        
        resid = ipermute(timeFirstResid, timeFirstIndices);
        residFrac = norm(resid(:)) / norm(data(:));
    else
        decomp = reshape(pinvA*timeFirstData(:,:), [2 szTFD(2:end)]);
    end
    
    szDecomp = size(decomp);
    
    complexDecomp = [1 1i]*decomp(:,:);
    f = ipermute(reshape(complexDecomp, [1 szDecomp(2:end)]), timeFirstIndices);
    
else
    % Output harmonic...
    
    nd = ndims(data);
    timeFirstIndices = [timeAxis, 1:timeAxis-1, timeAxis+1:nd];
    timeFirstData = permute(data, timeFirstIndices);
    szTFD = size(timeFirstData);
    szHarm = [numel(X.Frequency), szTFD(2:end)];
    
    freqs = reshape(X.Frequency, [], 1);
    timestamps = reshape(X.Time, 1, []);
    
    phaseFactors = exp(-1i*freqs*timestamps);
    
    harmonics = reshape(phaseFactors*timeFirstData(:,:), szHarm);
    f = ipermute(harmonics, timeFirstIndices) / numT;
    
end


function [f, freqs, resid, residFrac] = doFileSpectrum(fileName, X)
import t5.*;
file = OutputFile(fileName);
numT = file.numFramesAvailable;
numFields = numel(file.Fields);
numRegions = numel(file.Regions);

if (isempty(X.Frequency) && isempty(X.SteadyStateFrequency)) || ...
    strcmpi(X.Frequency, 'fft')

    % This is a beefed-up equivalent of outputFFT()
    
    if length(file.Durations) > 1
        error('Cannot do FFT on multi-duration outputs');
    end
    
    file.open();
    timeData = file.readFrames(numT);
    file.close();
    
    if numRegions > 1
        f = cell(size(timeData));
        for cc = 1:length(timeData)
            f{cc} = fft(timeData{cc}, [], ndims(timeData{cc})) / numT;
        end
    else
        f = {fft(timeData, [], ndims(timeData)) / numT};
    end
    
    deltaT = file.Durations{1}.Period * file.Dt;
    
    freqs = 2*pi*(0:numT-1) / (numT*deltaT);
    
    % Correct for field offset time
    if numFields > 1
        for rr = 1:numRegions
            for ff = 1:numFields
                offset = file.Fields{ff}.Offset(4);
                phaseFactors = exp(-1i*offset*file.Dt*freqs);

                for pp = 1:numel(freqs)
                    f{rr}(:,:,:,ff,pp) = f{rr}(:,:,:,ff,pp) * phaseFactors(pp);
                end
            end
        end
    else
        for rr = 1:numRegions
            offset = file.Fields{1}.Offset(4);
            phaseFactors = exp(-1i*offset*file.Dt*freqs);
            
            for pp = 1:numel(freqs)
                f{rr}(:,:,:,pp) = f{rr}(:,:,:,pp) * phaseFactors(pp);
            end
        end
    end
    
    if numel(f) == 1
        f = f{1};
    end
    
elseif ~isempty(X.SteadyStateFrequency)
    freq = X.SteadyStateFrequency;
    % do a single-frequency decomposition and provide a residual if
    % there are more than two timesteps present.
    
    if (numel(freq) > 1)
        error('Presently can only handle one steady-state frequency');
    end
    
    calcResid = (nargout > 2);
    
    file.open();
    timeData = file.readFrames(numT);
    file.close();
    if ~iscell(timeData)
        timeData = {timeData};
    end
    
    time = file.times(); % will be a cell array for multiple fields
    if ~iscell(time)
        time = {time};
    end
    
    % One region at a time:
    
    for rr = 1:numRegions
        st = size(timeData{rr});
        f{rr} = zeros(st(1:end-1));
        nd = ndims(timeData{rr});
        
        % Putting the time index first makes matrix multiplication simple.
        timeFirstData = permute(timeData{rr}, [nd, 1:nd-1]);
        
        if calcResid
            resid{rr} = 0*timeData{rr};
            
            sz = size(timeFirstData);
            timeFirstResids = 0*timeFirstData;
        end
        
        for ff = 1:numFields
            A = 2*[cos(freq*time{ff}'), -sin(freq*time{ff}')];
            pinvA = pinv(A);

            fieldData = timeFirstData(:,:,:,:,ff);
            sfd = size(fieldData);
            
            if calcResid
                timeFirstDecomp = pinvA*fieldData(:,:);
                
                timeFirstResids(:,:,:,:,ff) = reshape(A*timeFirstDecomp - fieldData(:,:),...
                    sfd);
                
                decomp = reshape(timeFirstDecomp, [2 sfd(2:end)]);
            else
                decomp = reshape(pinvA*fieldData(:,:), [2 sfd(2:end)]);
            end
            
            cosAndSin = ipermute(decomp, [nd, 1:nd-1]);

            f{rr}(:,:,:,ff) = cosAndSin(:,:,:,1) + ...
                1i*cosAndSin(:,:,:,2);
        end
        
        if calcResid
            resid{rr} = ipermute(timeFirstResids, [nd, 1:nd-1]);
            residFrac{rr} = norm(resid{rr}(:))/norm(timeData{rr}(:));
        end
    end
    
    if numel(f) == 1
        f = f{1};
        
        if calcResid
            resid = resid{1};
            residFrac = residFrac{1};
        end
    end
    
    if nargout > 1
        freqs = freq;
    end
    
    
else % OutputHarmonic-type functionality
    freqs = reshape(X.Frequency, 1, []); % a row vector!!
    deltaT = file.Durations{1}.Period * file.Dt;
    
    if length(file.Durations) > 1
        error('Cannot extract single frequencies from multi-duration output.');
    end
    
    harmonic = cell(numRegions,1);
    for rr = 1:numRegions
        dim{rr} = [file.Regions{rr}.Size, numFields];
        %if numFields > 1; dim{rr} = [file.Regions{rr}.Size, numFields];
        %else dim{rr} = file.Regions{rr}.Size;
        %end
        harmonic{rr} = zeros([prod(dim{rr}), numel(freqs)]);
    end
    
    frameNumber = 1;
    file.open;
    if (length(freqs) > 1)
        while frameNumber <= numT
            
            phaseFactors = exp(-1i*(frameNumber-1)*deltaT.*freqs);
            
            if numRegions == 1; frameData = {file.readFrames(1)};
            else frameData = file.readFrames(1);
            end
            
            for rr = 1:numRegions;
                % outer product: column vector * row vector
                dat = reshape(frameData{rr}, [], 1); % data is a column vector
                harmonic{rr} = harmonic{rr} + dat*phaseFactors/numT;
            end
            
            frameNumber = frameNumber + 1;
        end
    else
        while frameNumber <= numT
            phaseFactor = exp(-1i*(frameNumber-1)*deltaT*freqs);
            
            if numRegions == 1; frameData = {file.readFrames(1)};
            else frameData = file.readFrames(1);
            end
            
            for rr = 1:numRegions
                harmonic{rr} = harmonic{rr} + frameData{rr}(:)*phaseFactor/numT;
            end
            
            frameNumber = frameNumber+1;
        end
    end
    file.close
    
    for rr = 1:numRegions
        f{rr} = reshape(harmonic{rr}, [dim{rr} numel(freqs)]);
    end
    
    % Correct for field offset time
    if numFields > 1
        for ff = 1:numFields
            offset = file.Fields{ff}.Offset(4);
            phaseFactors = exp(-1i*offset*file.Dt*freqs);

            for pp = 1:numel(freqs)
                f{rr}(:,:,:,ff,pp) = f{rr}(:,:,:,ff,pp) * phaseFactors(pp);
            end
        end
    else
        offset = file.Fields{1}.Offset(4);
        phaseFactors = exp(-1i*offset*file.Dt*freqs);
        for pp = 1:numel(freqs)
            f{rr}(:,:,:,pp) = f{rr}(:,:,:,pp) * phaseFactors(pp);
        end
    end
    
    if numel(f) == 1
        f = f{1};
    end
    
end



