classdef DesignObject < handle

    properties
        Sim
        XML
    end
    
    
    methods
        function obj = DesignObject(sim, tag)
            
            if ~exist('tag', 'var')
                tag = '';
            end
            obj.Sim = sim;
            
            if ~isempty(tag)
                obj.Sim.Directory = ['sim_', tag];
                obj.Sim.OutputDirectory = ['output_', tag];
            else
                obj.Sim.Directory = 'sim';
                obj.Sim.OutputDirectory = 'output';
            end
            
            obj.XML = 'params.xml';
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
            obj.Sim.Grid.Meshes = obj.Sim.Grid.NodeGroup.meshes(designParameters);
        end
        
        function f = evaluateForward(obj)
            
            f = 0;
            
            for mm = 1:numel(obj.Sim.Grid.Measurements)
                meas = obj.Sim.Grid.Measurements{mm};
            
                fname = [obj.Sim.OutputDirectory filesep meas.filename];
            
                f = f + t6.adjoint.evalQuadraticFormFile(fname, meas.filters, meas.kernel);
            end
        end
        
        function [f, dfdp, dfdv, dvdp] = evaluate(obj, designParameters)
            
            [dfdp, dfdv, dvdp] = t6.adjoint.calculateAdjointSensitivity(...
                obj.Sim.Grid.NodeGroup, designParameters, ...
                obj.Sim.OutputDirectory);
                
            for mm = 1:numel(obj.Sim.Grid.Measurements)
                meas = obj.Sim.Grid.Measurements{mm};
                
                fname = [obj.Sim.OutputDirectory filesep meas.filename];
                
                f_ = t6.adjoint.evalQuadraticFormFile(fname, meas.filters, meas.kernel);
                
                if mm == 1
                    f = f_;
                else
                    f = f + f_;
                end
            end
        end
        
    end
    
end