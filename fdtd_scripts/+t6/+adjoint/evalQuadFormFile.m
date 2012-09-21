function [f, Df, f_w, freqs] = evalQuadFormFile(file, objFunStruct)
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
posInterp = pos;
timeVals = pos;

for ff = 1:of.numFields()
    pos{ff} = of.positions('Field', ff, 'InterpolateSpace', false);
    posInterp{ff} = of.positions('Field', ff);
    timeVals{ff} = of.times('Field', ff);
    
    xyzChars = 'XYZ';
    for xyz = 1:3
        if isempty(objFunStruct.(xyzChars(xyz)))
            objFunStruct.(xyzChars(xyz)){ff} = ...
                adjoint.interpolationMatrix(pos{ff}{xyz}, posInterp{ff}{xyz});
        end
    end
end

%% Call the worker function

if nargout > 2
    [f, Df, f_w, freqs] = adjoint.evalQuadForm(data, pos, timeVals, objFunStruct);
else
    [f, Df] = adjoint.evalQuadForm(data, pos, timeVals, objFunStruct);
end


