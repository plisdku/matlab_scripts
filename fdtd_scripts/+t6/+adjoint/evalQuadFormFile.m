function [f, Df, f_w, freqs, filteredData] = evalQuadFormFile(file, objFunStruct)
% evalQuadFormFile  Evaluate quadratic form on EM fields from file
%
% [f Df] = evalQuadForm(file, objectiveFunctionStruct)
%
% This provides the sum of squares of field values in the data file.
% Really it just calls evalQuadForm.  The reason this boilerplate is not in
% evalQuadForm is for ease of unit-testing, basically.

import t6.*

%% Read the file to get the goodies out
data = t6.readOutputFile(file, 'InterpolateSpace', false);
of = OutputFile(file);

pos = cell(of.numFields(), 1);
posInterp = of.positions(); % should just obtain natural sampling positions
timeVals = pos;

for ff = 1:of.numFields()
    pos{ff} = of.positions('Field', ff, 'InterpolateSpace', false);
    timeVals{ff} = of.times('Field', ff);
    
    for xyz = 1:3
        interpField = ['Interp', char('W'+xyz)];
        if isempty(objFunStruct.(interpField))
            objFunStruct.(interpField){ff} = ...
                adjoint.interpolationMatrix(pos{ff}{xyz}, posInterp{xyz});
        end
    end
end

%% Call the worker function

if nargout > 2
    [f, Df, f_w, freqs, filteredData] = adjoint.evalQuadForm(data, posInterp, timeVals, objFunStruct);
else
    [f, Df] = adjoint.evalQuadForm(data, posInterp, timeVals, objFunStruct);
end

%{
if nargout > 2
    [f, Df, f_w, freqs] = adjoint.evalQuadForm(data, pos, timeVals, objFunStruct);
else
    [f, Df] = adjoint.evalQuadForm(data, pos, timeVals, objFunStruct);
end
%}


