classdef AFPMultilayer < handle
    
    properties (SetAccess = private)
        boundariesE = [];
        boundariesH = [];
        layerOrigins = 0;
        dxyz = [0 0 0];
        dt = 0;
        permittivityFuncs = {@(z) 8.85e-12 * ones(size(z))};
        permeabilityFuncs = {@(z) 4e-7*pi * ones(size(z))};
        
        sourceIndices = [];
        nonSourceIndices = [];
        %sourcePhasors = [0 0 0 0];
        source = [0 0 0 0];  % forward Ex and Ey, backward Ex and Ey
        sourceDirection = [0 0 1];
        sourceOrigin = [0 0 0];
        
        % these anonymous functions execute many times faster than
        % actual class methods.
        eForward = @(xyz,layer) 12*(layer-1) + (xyz-1)*2 + 1;
        eBackward = @(xyz, layer) 12*(layer-1) + (xyz-1)*2 + 2;
        hForward = @(xyz, layer) 12*(layer-1) + (xyz-1)*2 + 7;
        hBackward = @(xyz, layer) 12*(layer-1) + (xyz-1)*2 + 8;
    end
    
    methods
        % Constructor
        function obj = AFPMultilayer(dxyz, dt, boundariesE, boundariesH,...
                permittivities, permeabilities)
            
            assert(length(dxyz) == 3);
            assert(length(dt) == 1);
            assert(length(boundariesE) == length(boundariesH));
            assert(length(permittivities) == length(permeabilities));
            assert(length(permittivities) == length(boundariesE) + 1);
            
            obj.dxyz = dxyz;
            obj.dt = dt;
            obj.boundariesE = reshape(boundariesE, 1, []);
            obj.boundariesH = reshape(boundariesH, 1, []);
            if ~isempty(boundariesE)
                obj.layerOrigins = [obj.boundariesE(1), obj.boundariesE];
            else
                obj.layerOrigins = 0;
            end
            obj.permittivityFuncs = permittivities;
            obj.permeabilityFuncs = permeabilities;
            
            % Source indices represent forward Ex and Ey, and
            % backward Ex and Ey.
            obj.sourceIndices = ...
                [obj.eForward(1,1), ...
                obj.eForward(2,1), ...
            	obj.eBackward(1, obj.numLayers()), ...
            	obj.eBackward(2, obj.numLayers())];
            obj.nonSourceIndices = sort(...
                setdiff(1:12*obj.numLayers(), obj.sourceIndices));
        end
        
        function n = numLayers(obj)
            n = length(obj.permittivityFuncs);
        end
        
        % addSource(o, dir, polX, srcX, polY, srcY)
        addSource(obj, origin, direction, polarization, source, varargin);
        
        % phasors(frequencies, cutoffFrequency)
        [fieldPhasors, z, zParallel, zNormal] = phasors(omegas, varargin);
        
        % solve({posEx posEy posEz}, {posHx posHy posHz},
        %   cutoffFrequency, coarsening)
        [E, H, timesteps] = solve(posE, posH, varargin);
        
    end
    
    methods(Access = private)
        phasors = solveStack(obj, permittivities, permeabilities, ...
            sourcePhasors, z, zParallel, zNormal);
        
        sourcePhasors = makeSourcePhasors(obj, coarsening);
    end
end