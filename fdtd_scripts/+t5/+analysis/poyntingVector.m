function [pv, timeOrFreq, positions] = poyntingVector(fileName, varargin)

X.SteadyStateFrequency = [];
X.Frequency = [];
X = parseargs(X, varargin{:});


file = t5.OutputFile(fileName);
if numel(file.Fields) ~= 6
    error('File must store ex, ey, ez, hx, hy and hz');
end

neededFields = {'ex', 'ey', 'ez', 'hx', 'hy', 'hz'};

for nn = 1:6
    if ~strcmp(file.Fields{nn}.Name, neededFields{nn})
        error('File must store ex, ey, ez, hx, hy and hz');
    end
end

if isempty(X.Frequency)
    % time-domain PV.
    % timeOrFreq = timesteps of E, and will be understood to be a little
    % wrong
    
    data = t5.readOutputFile(fileName);
    
    if iscell(data)
        pv = cell(size(file.Regions));
        for rr = 1:length(file.Regions)
            pv{rr} = cross(data{rr}(:,:,:,1:3,:), data{rr}(:,:,:,4:6,:), 4);
        end
    else
        pv = cross(data(:,:,:,1:3,:), data(:,:,:,4:6,:), 4);
    end
    
else
    % Frequency domain!
    
    [data, timeOrFreq] = analysis.spectrum(fileName, varargin{:});
    
    if iscell(data)
        pv = cell(size(file.Regions));
        for rr = 1:length(file.Regions)
            pv{rr} = 0.5*cross(data{rr}(:,:,:,1:3,:), ...
                conj(data{rr}(:,:,:,4:6,:)), 4);
        end
    else
        pv = 0.5*cross(data(:,:,:,1:3,:), conj(data(:,:,:,4:6,:)), 4);
    end
end

if nargout > 2
    [posX, posY, posZ] = file.positions();
    
    aField = 1;
    if iscell(data)
        positions = cell(numel(file.Regions),3);
        size(file.Regions)
        
        for rr = 1:numel(file.Regions)
            positions{rr,1} = posX{rr,aField};
            positions{rr,2} = posY{rr,aField};
            positions{rr,3} = posZ{rr,aField};
        end
        
        size(positions)
        
    else
        positions{1} = posX{aField};
        positions{2} = posY{aField};
        positions{3} = posZ{aField};
    end
end

