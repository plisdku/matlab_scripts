function [f, dfdp, dfdv] = designSensitivity(designObject, parameters, callback)
% [f dfdp dfdv] = designSensitivity(designObject, parameters, callback)

%global fs;
%global ps;

designObject.applyParameters(parameters);
designObject.writeSimulation(parameters, 'forward');

if strcmpi(computer, 'GLNXA64')
    snapdragon = 'env -u LD_LIBRARY_PATH snapdragon';
else
    snapdragon = 'snapdragon';
end

unixCall = sprintf('%s --geometry --sensitivity --outputDirectory %s %s/params.xml > out.txt', ...
    snapdragon, designObject.Sim.OutputDirectory, ...
    designObject.Sim.Directory);

exitval = unix(unixCall);

if exitval
    warning('Unix is unhappy with forward sim');
    vertSize = numel(designObject.Sim.Grid.NodeGroup.vertices(parameters));
    f = -realmax;
    dfdp = zeros(1, length(parameters));
    dfdv = zeros(vertSize, length(parameters));
    return;
end

designObject.writeSimulation(parameters, 'adjoint');

unixCall = sprintf('%s --adjoint --sensitivity --outputDirectory %s %s/params.xml > out.txt', ...
    snapdragon, designObject.Sim.OutputDirectory, ...
    designObject.Sim.Directory);
exitval = unix(unixCall);

if exitval
    warning('Unix is unhappy with adjoint sim');
    vertSize = numel(designObject.Sim.Grid.NodeGroup.vertices(parameters));
    f = -realmax;
    dfdp = zeros(1, length(parameters));
    dfdv = zeros(vertSize, length(parameters));
    return;
end

[f, dfdp, dfdv] = designObject.evaluate(parameters);

if isnan(f) || any(isnan(dfdp)) || any(isnan(dfdv))
    keyboard
end

%fs(end+1) = f;
%ps(end+1) = parameters;

if nargin > 2
    callback(designObject, parameters, f, dfdp, dfdv);
end

% Change the signs here to increase rather than decrease the objective function
dfdp = -dfdp;
f = -f;
