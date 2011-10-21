function nn = timesteps(obj)

nn = [];

for dd = 1:length(obj.Durations)
    nn = [nn, obj.Durations{dd}.First:obj.Durations{dd}.Period:obj.Durations{dd}.Last];
end
