% test ExtrudePath


import t6.*
import t6.model.*

% Profile!
xfunc = @(p) [p(1) p(2) p(3) p(4)]';
yfunc = @(p) [p(5) p(6) p(7) p(8)]';

pathfunc = @(p, t) [0 0 0]' + t*[p(9) p(10) 0]' + t^2 * [p(11) p(12) 0]';
upfunc = @(p, t) [0 0 1]';


extrusion = ExtrudePath('X', xfunc, 'Y', yfunc, 'Path', pathfunc, ...
    'V', upfunc);

%%

elem = @(A) A{1};
p0 = [1 2 2 1 1 1 2 2 3 4 -1 1]';

m0 = elem(extrusion.meshes(p0));

%%


figure(20); clf
flatPatch('Faces', m0.faces, 'Vertices', m0.patchVertices, ...
    'FaceColor', 'r');
camlight right
axis image vis3d
close(20);

%%

Dv_obj = m0.jacobian;

delta = 1e-8;
for nn = 1:length(p0)
    
    p1 = p0;
    p1(nn) = p0(nn) + delta;
    
    m1 = elem(extrusion.meshes(p1));
    
    Dv = (m1.vertices - m0.vertices)/delta;
    
    assert(norm(Dv - Dv_obj(:,nn), inf) < 1e-3);
end

%%

fprintf('ExtrudePath Jacobian test PASSED\n');