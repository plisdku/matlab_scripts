function tt = times(obj, varargin)

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
    tt = [tt, obj.Dt*nn{dd} + offset];
end


function validateArguments(obj, X)

if obj.numFields() > 1 && isempty(X.Field)
    error('More than one field is present in this file.  Please provide a Field.');
end







%{
function tt = times(obj)

tt = cell(length(obj.Fields), 1);

nn = obj.timesteps();
if ~iscell(nn)
    nn = {nn};
end

for dd = 1:length(obj.Durations)
for ff = 1:length(obj.Fields)
    offset = obj.Fields{ff}.Offset(4);
    tt{ff} = [tt{ff}, reshape(obj.Dt*(nn{dd} + offset), 1, [])];
end
end

if length(tt) == 1
    tt = tt{1};
end
%}
