% OutputFile

classdef OutputFile < handle
    properties (Hidden)
        FrameSize = 0; %  = nFields * total Yee cells, for field outputs
        FieldOffsetsInFrames = [0];  % measured in # values, not bytes
        RegionOffsetsInFields = [0]; % measured in # values, not bytes
        BytesPerValue = 4;
        FileHandle = -1;
    end
    
    % The following properties can be set only by class methods
    properties (SetAccess = private)
        FileName = '';
        Regions = {};  % struct: r.YeeCells, r.Size, r.Stride, r.NumYeeCells
        Durations = {}; % struct: d.First, d.Last, d.Period, d.NumTimesteps
        Materials = {};  % used in grid reports
        HalfCells = {};  % used in grid reports
        UnitVectors = [1 0 0; 0 1 0; 0 0 1];
        RunlineDirection = 0;
        Fields = {};
        Dxyz = [0 0 0];
        Dt = 0;
        DateString = '';
        SpecFileName = '';   % set in constructor
        %DataFileName = '';   % set in constructor
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
        
        n = numFramesAvailable(obj);
        [ii, jj, kk] = yeeCells(obj, varargin); % ('Regions', [1 3], 'Fields', {'ex', 'ez'})
        [xx, yy, zz] = positions(obj);
        nn = timesteps(obj);
        tt = times(obj);
        
        data = read(obj); % return fields by region
        
        open(obj)
        data = readFrames(obj, varargin)
        close(obj)
    end
    
    methods (Access = private)
        readSpecFile(obj)
    end
    
    methods (Static, Access = private)
        n = getIndexOfStringInCellArray(cellArray, stringToFind);
    end
end
