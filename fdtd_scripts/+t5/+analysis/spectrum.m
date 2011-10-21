function [f, freqs, resid, residFrac] = spectrum(fileName, varargin)
%
% [f, freqs] = spectrum(timeData, 'Dt', dt)
% [f, freqs] = spectrum(timeData, 'Dt', dt, 'Frequency', freq)
% [f, freqs] = spectrum(fileName)
% [f, freqs] = spectrum(fileName, 'Frequency', freq)
% [f, freqs] = spectrum(fileName, 'SteadyStateFrequency', freq)
% [f, freqs, residual, residualFraction] = spectrum(fileName, 'SteadyStateFrequency', freq)
%
%
% To replace outputFFT() or fft(), use spectrum(fileName) or spectrum(data)
% To replace outputHarmonic(), use spectrum(fileName, 'Frequency', freqs)
%

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

sz = size(data);
numT = sz(end);

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
    f = fft(data, [], ndims(data)) / numT;
    
elseif ~isempty(X.SteadyStateFrequency)
    
    if numel(X.SteadyStateFrequency) > 1
        error('Can only do one steady state frequency');
    end
    
    freqs = X.SteadyStateFrequency;
    
    calcResid = (nargout > 2);
    nd = ndims(data);
    timeFirstData = permute(data, [nd, 1:nd-1]);
    szTFD = size(timeFirstData);
    
    time = reshape(X.Time, 1, []);
    
    A = 2*[cos(freqs*time'), -sin(freqs*time')];
    pinvA = pinv(A);
    
    if calcResid
        timeFirstDecomp = pinvA*timeFirstData(:,:);
        
        timeFirstResid = reshape(A*timeFirstDecomp - timeFirstData(:,:),...
            szTFD);
        decomp = reshape(timeFirstDecomp, [2 szTFD(2:end)]);
        
        resid = ipermute(timeFirstResid, [nd, 1:nd-1]);
        residFrac = norm(resid(:)) / norm(data(:));
    else
        decomp = reshape(pinvA*timeFirstData(:,:), [2 szTFD(2:end)]);
    end
    
    szDecomp = size(decomp);
    
    complexDecomp = [1 1i]*decomp(:,:);
    f = ipermute(reshape(complexDecomp, [1 szDecomp(2:end)]), [nd, 1:nd-1]);
    
else
    % Output harmonic...
    
    nd = ndims(data);
    timeFirstData = permute(data, [nd, 1:nd-1]);
    szTFD = size(timeFirstData);
    szHarm = [numel(X.Frequency), szTFD(2:end)];
    
    freqs = reshape(X.Frequency, [], 1);
    times = reshape(X.Time, 1, []);
    
    phaseFactors = exp(-1i*freqs*times);
    
    harmonics = reshape(phaseFactors*timeFirstData(:,:), szHarm);
    f = ipermute(harmonics, [nd, 1:nd-1]) / numT;
    
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



