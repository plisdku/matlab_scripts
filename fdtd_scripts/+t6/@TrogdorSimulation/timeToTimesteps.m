function timesteps = timeToTimesteps(obj, duration, fieldTokens)
% timesteps = timeToTimesteps(duration, fieldTokens)
%
% duration should be a real-valued two-vector [t0 t1].

import t6.*

if ~iscell(fieldTokens)
    fieldTokens = {fieldTokens};
end

timesteps = [Inf, -Inf];

for ff = 1:numel(fieldTokens)
    offset = t6.fieldOffset(fieldTokens{ff}) .* [obj.Dxyz obj.Dt];
    
    nn = whichTimesteps(duration - offset(4), obj.Dt, obj.NumT);
    
    timesteps(1) = min(timesteps(1), nn(1));
    timesteps(2) = max(timesteps(2), nn(2));
end





function nn = whichTimesteps(duration, dt, numT)

nn = [floor(duration(1)/dt), ceil(duration(2)/dt)];

nn(1) = max(nn(1), 0);
nn(2) = min(nn(2), numT-1);

