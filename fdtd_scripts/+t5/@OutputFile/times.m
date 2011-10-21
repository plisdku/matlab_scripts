function tt = times(obj)

tt = cell(length(obj.Fields), 1);

nn = obj.timesteps();

for ff = 1:length(obj.Fields)
    offset = obj.Fields{ff}.Offset(4);
    tt{ff} = [tt{ff}, obj.Dt * (nn + offset)];
end

if length(tt) == 1
    tt = tt{1};
end
