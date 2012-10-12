
%checkSmall = @(a) assert(abs(a) < 1e-5);
%checkRelativelyClose = @(a, b) checkSmall(norm(a-b)/norm(a));
%checkClose = @(a, b) checkSmall(norm(a-b));

lambda = 632.8;
omega = 2*pi/lambda;
k0 = 2*pi/lambda;

eps1 = 1;
eps2 = -10 + 1i;

boundaries = 0;
epsr = [eps1, eps2];
mur = ones(size(epsr));

k_spp = omega*sqrt(eps1*eps2/(eps1+eps2));
n_spp = k_spp/k0;

kMin = 1.01*k0;
kMax = 4*k0;

%%

[ks, t22, kvals] = tmm.modesTM(boundaries, epsr, mur, omega, kMin, kMax);

assert(length(ks) == 1);
assert(abs(ks-k_spp) < 1e-3 * abs(k_spp));

fprintf('Surface plasmon dispersion test: passed\n');

%% Normalization test

boundaries = [0 1000];
epsr = [1, 5+0.1i, 2];
mur = ones(size(epsr));

[ks, t22, kvals] = tmm.modes(boundaries, epsr, mur, omega, 1.01*k0, 3*k0, 'tm');
%%
for nn = 1:length(ks)

    k = ks(nn);

    % Plot a mode

    outPos = linspace(boundaries(1)-500, boundaries(end)+500);

    [hx ey ez bigT bigR mux epsy epsz txMat] = tmm.solveTM(boundaries, epsr, mur, omega, ...
        k, outPos, true);
    
    energy = trapz(outPos, hx.*ez);
    assert(real(energy) > 0.99);
    assert(real(energy) < 1.01);
end

fprintf('Mode normalization test: passed\n');

%% TE Normalization test

epsr = [1, 5+0.1i, 1];
[ks, t22, kvals] = tmm.modes(boundaries, epsr, mur, omega, 1.01*k0, 3*k0, 'te');
%%

k = ks(1);
[ex hy hz] = tmm.solveTE(boundaries, epsr, mur, omega, k, [], true);

%%

for nn = 1
    k = ks(nn);
    outPos = linspace(boundaries(1)-500, boundaries(end)+500);
    [ex hy hz] = tmm.solveTE(boundaries, epsr, mur, omega, k, outPos, true);
    
    figure(11); clf
    plot(outPos, real(hz), outPos, imag(hz));
    title('Hy')
    figure(12); clf
    plot(outPos, real(ex), outPos, imag(ex));
    title('Ex')
    pause
    
    energy = trapz(outPos, -ex.*hz);
    assert(real(energy) > 0.99);
    assert(real(energy) < 1.01);
end

%{
for nn = 1:length(ks)

    k = ks(nn);

    % Plot a mode

    outPos = linspace(boundaries(1)-500e-9, boundaries(end)+500e-9);

    [hx ey ez bigT bigR mux epsy epsz txMat] = tmm.solveTM(boundaries, epsr, mur, omega, ...
        k, outPos, true);

    vertPos = linspace(0, 1000e-9);
    vertPhase = exp(1i*k*vertPos);
    
    hx2 = transpose(hx)*vertPhase;
    
    figure(200); clf
    imagesc_centered(outPos*1e9, vertPos*1e9, real(hx2)');
    axis xy;
    colormap orangecrush
    colorbar
    title(sprintf('Mode %i, t_{22} = %2.4g', nn, txMat(2,2)))
    xlabel('X (nm)')
    ylabel('Y (nm)')
    
    figure(201); clf
    plot(kvals, abs(t22), '-', ...
        real(k), spline(kvals, abs(t22), real(k)), 'o')
    xlabel('Effective mode index')
    ylabel('|t_{22}|');
    title('t_{22} for proposed mode');

    pause

end

%% Plot t22 in the complex plane

kR = linspace(kMin, kMax);
kI = linspace(0*k0, 0.3*k0);
t22s = zeros(length(kR), length(kI));

for xx = 1:length(kR)
    for yy = 1:length(kI)
        
        [hx ey ez bigT bigR mux epsy epsz txMat] = tmm.solveTM(...
            boundaries, epsr, mur, omega, kR(xx)+1i*kI(yy), [], true);
        
        t22s(xx,yy) = txMat(2,2);
    end
    fprintf('Done with x = %i\n', xx);
end

%%

figure
clf
imagesc(kR/k0, kI/k0, log(abs(t22s))');
colorbar
colormap hot
axis xy
xlabel('R')
ylabel('I')

hold on

for nn = 1:length(ks)
    plot(real(ks(nn))/k0, imag(ks(nn))/k0, 'wo')
    plot(real(ks(nn))/k0, imag(ks(nn))/k0, 'kx')
end

%}