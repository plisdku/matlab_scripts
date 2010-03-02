function valid = validatePMLParams(cellArray)
% Checks that cellArray contains only named parameters 'kappa', 'alpha' and
% 'sigma'.  Some valid cell arrays are:
%
% {'kappa', 5, 'sigma', 2, 'alpha', 0}
% {'kappa', 1}
% {'alpha', 2, 'kappa', 1}
% 
% etc.  Returns 1 for valid, 0 for invalid.

valid = 1;

X.kappa = '';
X.alpha = '';
X.sigma = '';
    
try
    X = parseargs(X, cellArray{:});
catch
    valid = 0;
end
