function tt = times(obj)

tt = cell(length(obj.Durations), length(obj.Fields));

nn = obj.timesteps();
if ~iscell(nn)
    nn = {nn};
end

for dd = 1:length(obj.Durations)
for ff = 1:length(obj.Fields)
    offset = t5.xml.fieldOffset(obj.Fields{ff}.Name);
    tt{dd,ff} = obj.Dt*(nn{dd} + offset(4));
end
end

if length(tt) == 1
    tt = tt{1};
end