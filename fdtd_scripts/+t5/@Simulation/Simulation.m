% Simulation

classdef Simulation < handle
    
    properties
        Grids = [];  % array of Grid instances
        Dxyz = [0 0 0];
        Dt = 0;
        NumT = 0;
        Directory = ''; % where to put all the auxiliary setup files
    end
    
    methods
        function obj = Simulation(dxyz, dt, numTimesteps)
            if length(dxyz) ~= 3
                error('dxyz needs to be a length-3 vector.');
            end
            obj.Dxyz = dxyz;
            obj.Dt = dt;
            obj.NumT = numTimesteps;
        end
        
        function str = directoryString(obj)
            if isempty(obj.Directory)
                str = '';
            else
                str = [obj.Directory, filesep()];
            end
        end
    end
end