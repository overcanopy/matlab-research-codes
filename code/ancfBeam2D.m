function [Me,fg,fi,energy] = ancfBeam2D(qe,p)
%ANCFBEAM2D Planar two-node ANCF beam element.
% Uses Green axial strain and a simple curvature-vector bending energy.

L = p.Le; EA = p.E*p.A; EI = p.E*p.I; rhoA = p.rho*p.A;
% 5-point Gauss integration on xi in [0,1]
xg = [-0.906179845938664,-0.538469310105683,0,0.538469310105683,0.906179845938664];
wg = [0.236926885056189,0.478628670499366,0.568888888888889,0.478628670499366,0.236926885056189];

Me = zeros(8); fg = zeros(8,1); fi = zeros(8,1);
energy.kinetic = 0; energy.strain = 0; energy.potential = 0;
gvec = [0;-p.g];

for k = 1:numel(xg)
    xi = 0.5*(xg(k)+1); w = 0.5*wg(k)*L;
    [S,Sx,Sxx] = ancfShape2D(xi,L);
    rx = Sx*qe; rxx = Sxx*qe; r = S*qe;
    eps = 0.5*(rx.'*rx - 1);

    Me = Me + rhoA*(S.'*S)*w;
    fg = fg + rhoA*S.'*gvec*w;
    fi = fi + (EA*eps*Sx.'*rx + EI*Sxx.'*rxx)*w;
    energy.strain = energy.strain + (0.5*EA*eps^2 + 0.5*EI*(rxx.'*rxx))*w;
    energy.potential = energy.potential - rhoA*(gvec.'*r)*w;
end
end
