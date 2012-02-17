function isStable = testStability(zNumer, zDenom, dt, dxyz, numCells)
%testStability  Predict whether a material model is expected to be
% stable in FDTD

isStable = true;

invDxyz2 = 1 ./ (dxyz.^2);
invDxyz2(numCells <= 1) = 0;

kMin = 0;
kMax = 4*sum(invDxyz2);
delta = 1e-5;

%nn = 1; dd = 1;
kVals = linspace(kMin, kMax);

for k = kVals
    stabilityPoly = polyadd(conv([1 -2 1], zNumer), ...
        dt^2*conv([1 0], zDenom)*k);
    
    rr = roots(stabilityPoly);
    %{
    figure(99); clf
    plot(cos(linspace(0,2*pi)), sin(linspace(0, 2*pi)), 'Color', [0.7 0.7 1.0])
    hold on
    plot(real(rr(abs(rr) > 1+delta)), imag(rr(abs(rr) > 1+delta)), 'rx');
    plot(real(rr(abs(rr) <= 1+delta)), imag(rr(abs(rr) <= 1+delta)), 'go')
    %plot(rr, 'o')
    xlim([-2 2])
    ylim([-2 2])
    %}
    
    if any(abs(rr) > 1 + delta)
        isStable = false;
        %fprintf('Unstable\n');
    end
end