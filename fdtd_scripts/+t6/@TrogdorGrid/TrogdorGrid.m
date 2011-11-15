classdef TrogdorGrid < handle
    
    properties
        Name = '';
        Assembly = {}; % cell array of various assembly structs
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
        yeeCells = extent(obj);
    end
end