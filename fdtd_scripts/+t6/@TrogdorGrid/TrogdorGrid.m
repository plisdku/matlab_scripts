classdef TrogdorGrid < handle
    
    properties
        Name = '';
        Meshes = {};
        ParameterizedMeshes = {};
        Background = [];
        PML = [0 0 0 0 0 0];
        PMLParams = struct('alpha', '', 'kappa', '', 'sigma', '');
        
        Outputs = {};
        HardSources = {};
        SoftSources = {};
        CurrentSources = {};
        TFSFSources = {};
        CustomTFSFSources = {};
        Links = {};
        ModelReports = {};
        GridReports = {};
        YeeCells = [0 0 0 0 0 0]; % one cell
        Origin = [0 0 0];
    end
    
    methods
        function nc = numCells(obj, varargin)
            if numel(varargin) == 0
                nc = obj.YeeCells([4:6]) - obj.YeeCells([1:3]) + 1;
            else
                nc = obj.numCells();
                nc = nc(varargin{1});
            end
        end
    end
end