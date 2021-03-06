classdef TrogdorSimulation < handle
    
    properties
        Materials = {};
        Grid = t6.TrogdorGrid;
        %Grids = {}
        %CurrentGrid % handle to, well, the current grid
        ElectromagneticMode = '3d'
        Dxyz = 0
        Dt = 0
        NumT = 0
        NumCells = [1 1 1];
        NonPMLBounds = [];
        OuterBounds = [];
        Precision = 'float64';
        Directory = '';
        OutputDirectory = '';
    end
    
    methods
        function obj = TrogdorSimulation()
            obj.Grid = t6.TrogdorGrid;
        end
        
        function str = directoryString(obj)
            if isempty(obj.Directory)
                str = '';
            else
                str = [obj.Directory, filesep()];
            end
        end
        
        function str = outputDirectoryString(obj)
            if isempty(obj.OutputDirectory)
                str = '';
            else
                str = [obj.OutputDirectory, filesep()];
            end
        end
        
        timesteps = timeToTimesteps(obj, timeBounds, fieldTokens);
        yeeCells = boundsToYee(obj, bounds, fieldTokens); 
        y = yeeCells(obj, bounds);
        v = extendIntoPML(obj, verts);
    end
    
    %methods (Access = private)
    %    function obj = TrogdorSimulation(Dxyz, Dt, NumT)
    %    end
    %end 
    
    methods (Static)
        %{
        function clear()
            global TROGDOR_SIMULATION
            %if exist('TROGDOR_SIMULATION', 'var')
            %    delete(TROGDOR_SIMULATION);
            %end
            TROGDOR_SIMULATION = t6.TrogdorSimulation();
        end
        
        function obj = instance()
            global TROGDOR_SIMULATION
            if ~isa(TROGDOR_SIMULATION, 't6.TrogdorSimulation')
                TROGDOR_SIMULATION = t6.TrogdorSimulation();
            end
            obj = TROGDOR_SIMULATION;
        end
        %}
        %write(fileName);
    end
end
