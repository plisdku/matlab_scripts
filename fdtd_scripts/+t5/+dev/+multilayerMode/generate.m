function generate(direction, omega, dielectricEpsilon, modeIndex)
% generate     Create TM waveguide mode source for multilayer stack.
% Usage: generate(direction, frequency, dielectricEpsilon)
%        generate(direction, frequency, dielectricEpsilon, modeIndex)
%

tfsfreq_0;

posE = {afp.posEx, afp.posEy, afp.posEz};
posH = {afp.posHx, afp.posHy, afp.posHz};

[afpE afpH] = phasors(direction, omega, dielectricEpsilon, 1, posE, posH);

lambda = 2*pi*3e8./omega;
k0 = 2*pi/lambda;
kMin = 1.01*k0;
kMax = 3*k0;

n0 = 350;
srcRamp = @(n) 0.5*(1+erf(0.008*(n-n0)));

%% Transform to the time domain and write to file

fprintf('Writing fields to file.\n');
srcE = @(xyz,n) real(srcRamp(n) .* exp(-1i*n*afp.dt*omega) .* afpE{xyz});
srcH = @(xyz,n) real(srcRamp(n+0.5) .* exp(-1i*(n+0.5)*afp.dt*omega) .* afpH{xyz});

fh = fopen(afp.inputFile, 'w');

for tt = 1:afp.numT
    
    for xyz = 1:3
        fwrite(fh, srcH(xyz,tt), 'float32');
    end
    
    for xyz = 1:3
        fwrite(fh, srcE(xyz,tt), 'float32');
    end
    
end
fclose(fh);

