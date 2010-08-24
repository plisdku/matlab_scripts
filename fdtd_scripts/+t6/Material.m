% Material

classdef Material < handle
    
    properties
        name = '';
        model = '';
        PMLParams = struct('kappa', '', 'alpha', '', 'sigma', '');
        numerator = [];
        denominator = [];
    end
    
    methods
        
    end
    
end