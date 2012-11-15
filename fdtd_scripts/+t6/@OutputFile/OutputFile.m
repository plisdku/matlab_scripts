% OutputFile
% This is my output file

classdef OutputFile < handle
    % Output file!
    properties (Hidden)
        FrameSize = 0; %  = nFields * total Yee cells, for field outputs
        FieldOffsetsInFrames = [0];  % measured in # values, not bytes
        RegionOffsetsInFields = [0]; % measured in # values, not bytes
        BytesPerValue = 4;
        FileHandle = -1;
        SavedFrame = {};
        SavedFrameNumber = -1;
        NextFrameNumber = -1;
        SampleTimes = [];
    end
    
    % The following properties can be set only by class methods
    properties (SetAccess = private)
        FileName = '';
        Regions = struct('YeeCells', [], 'Size', [], 'Stride', [], ...
            'NumYeeCells', [], 'Bounds', []);
        Durations = {}; % struct: d.First, d.Last, d.Period, d.NumTimesteps
        Materials = {};  % used in grid reports
        HalfCells = {};  % used in grid reports
        UnitVectors = [1 0 0; 0 1 0; 0 0 1];
        RunlineDirection = 0;
        Fields = {};
        Dxyz = [0 0 0];
        Dt = 0;
        Origin = [0 0 0];
        DateString = '';
        SpecFileName = '';   % set in constructor
        TrogdorVersionString = '';
        Precision = 'float32';
    end
    
    methods
        % Constructor
        function obj = OutputFile(fileName)
            obj.FileName = fileName;
            
            % Check for existence of the specfile (required)
            obj.SpecFileName = [fileName, '.txt'];
            if ~exist(obj.SpecFileName)
                error(['Could not open spec file ', obj.SpecFileName, '.']);
            end
            
            % Check for existence of the data file (optional)
            if ~exist(fileName, 'file')
                warning(['Warning: data file ', fileName, ' does not exist.']);
            end
            
            obj.readSpecFile();
        end
        
        function nR = numRegions(obj)
            nR = size(obj.Regions.YeeCells, 1);
        end
        
        function nF = numFields(obj)
            nF = numel(obj.Fields);
        end
        
        function yesNo = hasBounds(obj)
            yesNo = ~any(isnan(obj.Regions.Bounds(:)));
        end
        
        % obtain # of samples per field per timestep
        function fv = fieldValues(obj)
            fv = 0;
            for rr = 1:obj.numRegions
                fv = fv + prod(obj.Regions.Size(rr,:));
            end
        end
        
        n = numFramesAvailable(obj);
        ijk = yeeCells(obj, varargin); % ('Regions', [1 3], 'Fields', {'ex', 'ez'})
        xyz = positions(obj, varargin);
        nn = timesteps(obj);
        tt = times(obj, varargin);
        
        open(obj)
        [data, positions] = readFrames(obj, varargin)
        close(obj)
        
        seekFrame(obj, frame)
    end
    
    methods (Access = private)
        readSpecFile(obj)
        data = readFrames_SeparateRegions(obj, numFrames);
        data = readFrames_RegionsTogether(obj, numFrames);
        pos = naturalSamplingPositions(obj, bounds, region, numSamples);
        tt = naturalSamplingTimes(obj, duration, numSamples);
    end
    
    methods (Static, Access = private)
        n = getIndexOfStringInCellArray(cellArray, stringToFind);
    end
    
    methods (Static)
        [starts, ends, lengths] = chunkTimesteps(nFirst, nLast, ...
            frameSize, bytes);
    end
end