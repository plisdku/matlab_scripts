function jj = jacobian(func, parameters, varargin)
%
% jj = jacobian(func, parameters)
% jj = jacobian(func, parameters, deltas)
%

if nargin > 2
    deltas = varargin{1};
else
    deltas = 1e-5*parameters;
    deltas(deltas == 0) = 1e-9;
end

if ndims(deltas) ~= ndims(parameters)
    error('Deltas must have same size as parameters.');
end

if any(size(deltas) ~= size(parameters))
    error('Deltas must have same size as parameters.');
end

% Evaluate the function once.
f0 = func(parameters);
if size(f0,1) < size(f0,2)
    f0 = f0';
end

% Re-evaluate the function at perturbed positions to obtain a sensitivity.
jj = sparse(size(f0, 1), length(parameters)*size(f0,2));

for pp = 1:length(parameters)
    
    params1 = parameters;
    params1(pp) = params1(pp) + deltas(pp);
    f1 = func(params1);
    
    if size(f1,1) < size(f1,2)
        f1 = f1';
    end
    
    jj(:, (1+(pp-1)*size(f0,2)):(pp*size(f0,2))) = (f1-f0)/deltas(pp);
end



