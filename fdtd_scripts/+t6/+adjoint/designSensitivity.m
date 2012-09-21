function [f dfdp dfdv] = designSensitivity(designObject, parameters, callback)

%global fs;
%global ps;

designObject.applyParameters(parameters);
designObject.writeSimulation(parameters, 'forward');
exitval = unix('snapdragon --geometry --sensitivity --outputDirectory output sim/params.xml > out.txt');

if exitval
    error('Unix is unhappy with forward sim');
end

designObject.writeSimulation(parameters, 'adjoint');
exitval = unix('snapdragon --adjoint --sensitivity --outputDirectory output sim/params.xml > out.txt');

if exitval
    error('Unix is unhappy with adjoint sim');
end

[f dfdp dfdv] = designObject.evaluate(parameters);

%fs(end+1) = f;
%ps(end+1) = parameters;

if nargin > 2
    callback(designObject, parameters, f, dfdp, dfdv);
end

% Change the signs here to increase rather than decrease the objective function
dfdp = -dfdp;
f = -f;