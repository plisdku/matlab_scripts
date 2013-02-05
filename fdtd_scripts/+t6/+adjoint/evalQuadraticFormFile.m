function [f Df] = evalQuadraticFormFile(file, userFilters, userKernel)

import t6.*

data = readOutputFile(file, 'InterpolateSpace', false);
of = OutputFile(file);

pos = cell(of.numFields(), 1);
posInterp = of.positions();
timeVals = pos;

[~, tInterp] = timeExtent(of);

filters = {};

% space filters
for xyz = 1:3
    spaceMatrices = cell(1, of.numFields());
    timeMatrices = cell(1, of.numFields());
    
    for ff = 1:of.numFields()
    
        pos{ff} = of.positions('Field', ff, 'InterpolateSpace', false);
        
        spaceMatrices{ff} = adjoint.interpolationMatrix(pos{ff}{xyz}, ...
            posInterp{xyz});
    end
    
    filters = [filters, {adjoint.lf('MatrixArray', spaceMatrices, 'Dim', xyz)}];
end

% time filters
for ff = 1:of.numFields()
    timeVals{ff} = of.times('Field', ff);
    timeMatrices{ff} = adjoint.interpolationMatrix(timeVals{ff}, ...
        tInterp);
end
filters = [filters, {adjoint.lf('MatrixArray', timeMatrices, 'Dim', 5)}];

% tack on the user filters now
filters = [filters, userFilters];

%%

[f Df] = adjoint.evalQF(data, filters, userKernel, posInterp{1}, posInterp{2}, ...
    posInterp{3}, tInterp);


function [tSpan, tSamples] = timeExtent(of)
% timeExtent  Return a span of times covered by all field measurements

tSpan = [-Inf, Inf];

for ff = 1:of.numFields()
    tt = of.times('Field', ff);
    tSpan(1) = max(tSpan(1), tt(1));
    tSpan(2) = min(tSpan(2), tt(end));
    
    dt = tt(2)-tt(1);
end

numSamples = ceil( (tSpan(2)-tSpan(1))/dt );
tSamples = linspace(tSpan(1), tSpan(2), numSamples);

