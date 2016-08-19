function meshes = endSim(varargin)
    
    global LL_MODEL;
    
    if strcmpi(LL_MODEL.physics, 'maxwell')
        ll.endSim_maxwell(varargin{:});
    elseif strcmpi(LL_MODEL.physics, 'poisson')
        meshes = ll.endSim_poisson(varargin{:});
    end
end


