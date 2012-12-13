classdef DesignObject < handle

    properties
        Sim
        XML
    end
    
    
    methods
        function obj = DesignObject(sim)
            obj.Sim = sim;
            obj.XML = 'params.xml';
            obj.Sim.Directory = 'sim';
            obj.Sim.OutputDirectory = 'output';
        end
        
        function writeSimulation(obj, designParameters, mode)
            import com.mathworks.xml.XMLUtils.*;
            
            if ~strcmpi(mode, 'forward') && ~strcmpi(mode, 'adjoint')
                error('Mode must be forward or adjoint');
            end
            
%            if numel(obj.Sim.Grids) > 1
%                error('Trogdor 6 does not support multiple grids');
%            end
            
            if ~isempty(obj.Sim.Directory) && ~exist(obj.Sim.Directory, 'dir')
                try mkdir(obj.Sim.Directory)
                catch exception
                    error('Could not create helper directory!');
                end
            end

            if ~isempty(obj.Sim.OutputDirectory) && ~exist(obj.Sim.OutputDirectory, 'dir')
                try mkdir(obj.Sim.OutputDirectory)
                catch exception
                    error('Could not create output directory!');
                end
            end
            
%            assert(numel(obj.Sim.Grids) == 1);
%            if ~isempty(obj.Sim.Grids{1}.Links)
%                warning('Ignoring Link objects');
%            end
            
            doc = t6.xml.generateXML(obj.Sim, designParameters, mode);
            xmlwrite([obj.Sim.Directory, filesep, obj.XML], doc);
        end
        
        function applyParameters(obj, designParameters)
            
            %for gg = 1:length(obj.Sim.Grids)
            %    grid = obj.Sim.Grids{gg};
                
                obj.Sim.Grid.Meshes = obj.Sim.Grid.NodeGroup.meshes(designParameters);
            %end
            
        end
        
        function [f dfdp dfdv dvdp] = evaluate(obj, designParameters)
            
            meas = obj.Sim.Grid.Measurement;
            
            fname = [obj.Sim.OutputDirectory filesep meas.filename];
            
            %f = t6.adjoint.evalQuadFormFile(fname, meas.function);
            f = evalQuadraticFormFile(fname, meas.function);
            
            [dfdp dfdv dvdp] = t6.adjoint.calculateAdjointSensitivity(...
                obj.Sim.Grid.NodeGroup, designParameters);
            
        end
        
    end
    
end