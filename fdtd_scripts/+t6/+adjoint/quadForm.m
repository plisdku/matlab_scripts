function X = quadForm(varargin)
% quadForm  Create quadratic objective function for optimization
%
% X         handle @(x) f(x), to be integrated against the measured field
%           (including dx measure, so it's resolution-independent!)
% Y, Z      like X
% XYZ       handle @(x,y,z) f(x,y,z) to be multiplied pointwise with field
% T         handle @(t) f(t), to be integrated against the measured field
%           (including dt measure, so it's resolution-independent!)
% Kernel    NfxNf matrix connecting field components, e.g. to calculate D*E
%           or ExH

X.X = [];
X.Y = [];
X.Z = [];
X.XYZ = [];
X.T = [];
X.Kernel = [];
X.InterpX = []; % these interpolation matrices aren't the user's problem!
X.InterpY = [];
X.InterpZ = [];
X.Factor = 1;

X = parseargs(X, varargin{:});
validateArguments(X);

%noms = fieldnames(X);
%for nn = 1:length(noms)
%    if isempty(X.(noms{nn}))
%        X = rmfield(X, noms{nn});
%    end
%end

end




function validateArguments(X)

tmpString = '';

% Validate X, Y, Z filters
for xyz = 'XYZ'
if ~isempty(X.(xyz))
    
    if iscell(X.(xyz))
        if ~exist('nFields', 'var')
            nFields = numel(X.(xyz));
            tmpString = xyz;
        elseif nFields ~= numel(X.(xyz))
            error('%s filter has different number of elements than %s filter',...
                xyz, tmpString);
        end
        
        %for nn = 1:numel(X.(xyz))
        %    if ~isa(X.(xyz){nn}, 'function_handle')
        %        error('%s filter must be a function handle or cell array of function handles', xyz);
        %    end
        %end
    else
        %if ~isa(X.(xyz), 'function_handle')
        %    error('%s filter must be a function handle or cell array of function handles', xyz);
        %end
    end
end
end

% Validate XYZ filter
if ~isempty(X.XYZ)
    if iscell(X.XYZ)
        if ~exist('nFields', 'var')
            nFields = numel(X.XYZ);
            tmpString = 'XYZ';
        elseif nFields ~= numel(X.XYZ)
            error('XYZ filter has different number of elements than %s filter',...
                tmpString);
        end
        
        for nn = 1:numel(X.XYZ)
            if ~isa(X.XYZ{nn}, 'function_handle')
                error('XYZ filter must be a function handle or cell array of function handles');
            end
        end
    else
        if ~isa(X.XYZ, 'function_handle')
            error('XYZ filter must be a function handle or cell array of function handles');
        end
    end
end

% Validate T filters

if ~isempty(X.T)
    if iscell(X.T)
        if ~exist('nFields', 'var')
            nFields = numel(X.T);
            tmpString = 'T';
        elseif nFields ~= numel(X.T)
            error('T filter has different number of elements than %s filter',...
                tmpString);
        end
        
        for nn = 1:numel(X.T)
            if ~isa(X.T{nn}, 'function_handle')
                error('T filter must be a function handle or cell array of function handles');
            end
        end
    else
        if ~isa(X.T, 'function_handle')
            error('T filter must be a function handle or array of function handles');
        end
    end
end

% Validate kernel

if ~isempty(X.Kernel)
    if exist('nFields', 'var')
        validateattributes(X.Kernel, {'numeric'}, {'size', [nFields nFields]});
    else
        validateattributes(X.Kernel, {'numeric'}, {'2d'}, 'quadraticObjective');
        nFields = size(X.Kernel, 1);
        if size(X.Kernel,1) ~= size(X.Kernel,2)
            error('Kernel must be a square matrix');
        end
    end
end

end