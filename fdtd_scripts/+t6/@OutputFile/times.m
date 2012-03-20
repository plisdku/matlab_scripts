function tt = times(obj, varargin)
% times()
% times('Field', 1)

X.Field = [];
X = parseargs(X, varargin{:});
validateArguments(obj, X);

if isempty(X.Field)
    X.Field = 1;
end

nn = obj.timesteps();
if ~iscell(nn)
    nn = {nn};
end
offset = obj.Fields{X.Field}.Offset(4);

tt = [];
for dd = 1:length(obj.Durations)
    tt = [tt, obj.Dt*(nn{dd} + offset)];
end


function validateArguments(obj, X)

if obj.numFields() > 1 && isempty(X.Field)
    error('More than one field is present in this file.  Please provide a Field.');
end

