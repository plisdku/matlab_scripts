function [f Df] = evalQuadraticFormFile(file, userFilters)

import t6.*

data = readOutputFile(file, 'InterpolateSpace', false);
of = OutputFile(file);

pos = cell(of.numFields(), 1);
posInterp = of.positions();
timeVals = pos;

filters = {};

for xyz = 1:3
    matrices = cell(1, of.numFields());
    
    for ff = 1:of.numFields()
    
        pos{ff} = of.positions('Field', ff, 'InterpolateSpace', false);
        timeVals{ff} = of.times('Field', ff);
        
        matrices{ff} = adjoint.interpolationMatrix(pos{ff}{xyz}, ...
            posInterp{xyz});
    end
    
    filters = [filters, {adjoint.lf('MatrixArray', matrices, 'Dim', xyz)}];
end

filters = [filters, userFilters];

%%

[f Df] = adjoint.evalQF(data, filters, posInterp{1}, posInterp{2}, ...
    posInterp{3}, timeVals);


