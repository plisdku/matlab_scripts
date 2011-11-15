function tt = times(obj)

tt = cell(length(obj.Durations), length(obj.Fields));

nn = obj.timesteps();
if ~iscell(nn)
    nn = {nn};
end

for dd = 1:length(obj.Durations)
for ff = 1:length(obj.Fields)
    offset = obj.Fields{ff}.Offset(4);
    tt{dd,ff} = obj.Dt*(nn{dd} + offset);
end
end

if length(tt) == 1
    tt = tt{1};
end