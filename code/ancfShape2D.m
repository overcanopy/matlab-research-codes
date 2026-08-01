function [S,Sx,Sxx] = ancfShape2D(xi,L)
%ANCFSHAPE2D Cubic planar ANCF beam interpolation matrices.
% xi is the element coordinate in [0,1]. The nodal coordinates are
% [x1;y1;x1_s;y1_s;x2;y2;x2_s;y2_s], where derivatives are with respect
% to the material coordinate x.

h1 = 1 - 3*xi^2 + 2*xi^3;
h2 = L*(xi - 2*xi^2 + xi^3);
h3 = 3*xi^2 - 2*xi^3;
h4 = L*(-xi^2 + xi^3);

h1x = (-6*xi + 6*xi^2)/L;
h2x = 1 - 4*xi + 3*xi^2;
h3x = ( 6*xi - 6*xi^2)/L;
h4x = -2*xi + 3*xi^2;

h1xx = (-6 + 12*xi)/L^2;
h2xx = (-4 + 6*xi)/L;
h3xx = ( 6 - 12*xi)/L^2;
h4xx = (-2 + 6*xi)/L;

S = [h1 0 h2 0 h3 0 h4 0;
     0 h1 0 h2 0 h3 0 h4];
Sx = [h1x 0 h2x 0 h3x 0 h4x 0;
      0 h1x 0 h2x 0 h3x 0 h4x];
Sxx = [h1xx 0 h2xx 0 h3xx 0 h4xx 0;
       0 h1xx 0 h2xx 0 h3xx 0 h4xx];
end
