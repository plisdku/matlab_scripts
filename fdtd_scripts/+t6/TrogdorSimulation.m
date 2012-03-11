classdef TrogdorSimulation < handle
    
    properties
        Materials = {};
        Grids = {}
        CurrentGrid % handle to, well, the current grid
        Dxyz = 0
        Dt = 0
        NumT = 0
        NumCells = [1 1 1];
        NonPMLBounds = [];
        OuterBounds = [];
        Precision = 'float32';
        Directory = '';
        OutputDirectory = '';
    end
    
    methods
        function grid = currentGrid(obj)
            if ~isa(obj.CurrentGrid, 't6.TrogdorGrid')
                error('Please call addGrid to create a new grid.');
            end
            grid = obj.CurrentGrid;
        end
        
        function setCurrentGrid(obj, grid)
            if ~isa(grid, 't6.TrogdorGrid')
                error('Not a valid grid handle.');
            end
            obj.CurrentGrid = grid;
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
    end
    
    methods (Access = private)
        function obj = TrogdorSimulation(Dxyz, Dt, NumT)
        end
    end 
    
    methods (Static)
        function clear()
            global TROGDOR_SIMULATION
            TROGDOR_SIMULATION = t6.TrogdorSimulation();
        end
        
        function obj = instance()
            global TROGDOR_SIMULATION
            if ~isa(TROGDOR_SIMULATION, 't6.TrogdorSimulation')
                TROGDOR_SIMULATION = t6.TrogdorSimulation();
            end
            obj = TROGDOR_SIMULATION;
        end
        
        write(fileName);
    end
end
