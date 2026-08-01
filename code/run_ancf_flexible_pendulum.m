function out = run_ancf_flexible_pendulum(doPlot)
%RUN_ANCF_FLEXIBLE_PENDULUM Flexible pendulum benchmark using planar ANCF.
% The beam is initially horizontal. The root position is pinned, while the
% root slope remains free, so the body undergoes large rotation and elastic
% deformation under gravity.

if nargin < 1, doPlot = true; end

p.L = 1.0;          % m
p.ne = 4;
p.Le = p.L/p.ne;
p.E = 2.1e9;        % Pa
p.rho = 7800;       % kg/m^3
p.A = 2.0e-5;       % m^2
p.I = 1.67e-11;     % m^4
p.g = 9.81;         % m/s^2

nn = p.ne + 1;
ndof = 4*nn;
q0 = zeros(ndof,1);
for i = 1:nn
    x = (i-1)*p.Le;
    id = 4*(i-1) + (1:4);
    q0(id) = [x;0;1;0];
end
v0 = zeros(ndof,1);

[M,fg] = assembleConstantMatrices(p);
fixed = [1 2];
free = setdiff(1:ndof,fixed);
Mff = M(free,free);

z0 = [q0(free);v0(free)];
opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',2e-3);
[t,z] = ode15s(@rhs,[0 0.35],z0,opts);

qhist = zeros(numel(t),ndof);
qhist(:,fixed) = repmat(q0(fixed).',numel(t),1);
qhist(:,free) = z(:,1:numel(free));

out.t = t;
out.q = qhist;
out.tip = qhist(:,4*(nn-1)+(1:2));
out.parameters = p;
out.massMatrix = M;
out.gravityForce = fg;

if doPlot
    figure('Name','ANCF flexible pendulum');
    plot(out.tip(:,1),out.tip(:,2),'LineWidth',1.5); axis equal; grid on
    xlabel('Tip x [m]'); ylabel('Tip y [m]');
    title('ANCF flexible pendulum: tip trajectory');

    figure('Name','ANCF final configuration'); hold on; grid on; axis equal
    plotConfiguration(q0,p,'--');
    plotConfiguration(qhist(end,:).',p,'-');
    xlabel('x [m]'); ylabel('y [m]');
    legend('Initial','Final','Location','best');
    title(sprintf('Final configuration at t = %.3f s',t(end)));
end

    function dz = rhs(~,zv)
        q = q0;
        q(free) = zv(1:numel(free));
        v = zv(numel(free)+1:end);
        fi = assembleInternalForce(q,p);
        a = Mff \ (fg(free)-fi(free));
        dz = [v;a];
    end
end

function [M,fg] = assembleConstantMatrices(p)
nn = p.ne+1; ndof = 4*nn;
M = zeros(ndof); fg = zeros(ndof,1);
for e = 1:p.ne
    ed = elementDofs(e);
    qe0 = [0;0;1;0;p.Le;0;1;0];
    [Me,fge] = ancfBeam2D(qe0,p);
    M(ed,ed) = M(ed,ed) + Me;
    fg(ed) = fg(ed) + fge;
end
end

function fi = assembleInternalForce(q,p)
fi = zeros(size(q));
for e = 1:p.ne
    ed = elementDofs(e);
    [~,~,fie] = ancfBeam2D(q(ed),p);
    fi(ed) = fi(ed) + fie;
end
end

function ed = elementDofs(e)
% Node e and node e+1, four coordinates per node.
ed = [4*(e-1)+(1:4), 4*e+(1:4)];
end

function plotConfiguration(q,p,ls)
pts = [];
for e = 1:p.ne
    ed = elementDofs(e); qe = q(ed);
    xis = linspace(0,1,25);
    xy = zeros(2,numel(xis));
    for k = 1:numel(xis)
        S = ancfShape2D(xis(k),p.Le);
        xy(:,k) = S*qe;
    end
    if e > 1, xy = xy(:,2:end); end
    pts = [pts xy]; %#ok<AGROW>
end
plot(pts(1,:),pts(2,:),ls,'LineWidth',1.5);
end
