%{
My name is Robert Richardson (robert.richardson@eng.ox.ac.uk) and I hold 
the MIT license for this code. An accompanying paper for this code entitled
"On-board monitoring of 2-D spatially-resolved temperatures
in cylindrical lithium-ion batteries: Part I. Low-order thermal modelling"
has been submitted for publication in the Journal of Power Sources. A 
preprint for this paper can be found at: http://arxiv.org

I would ask that you cite this paper if you want to use this code for 
your own research. For further details on the work of the Energy Power 
Group at Oxford, please see epg.eng.ox.ac.uk.
%}

function [x,Phi_x] = func_basis_fn (num,am,bm,ap,bp,N)
% This function first calculates the basis functions that satisfy specified
% boundary conditions, and then evaluates these functions at a set of
% Gauss-Lobatto points.
% Reference: page 148, Spectral methods algorithms, analysis and
% applications.

% am*u(-1) + bm*u'(-1) = 0,  ap*u(1) + bp*u'(1) = 0

theta = (0:num)/num*pi;
x = cos(theta);                      % Gauss-Lobatto points
Phi_x = zeros(N+1,num+1);

if (am^2+bm^2+ap^2+bp^2==0)
    for k = 0:1:num
        Phi_x(k+1,:) = cos(k*x);
    end;
else
    for k = 0:1:N
        DETk = 2*ap*am+((k+1)^2+(k+2)^2)*(am*bp-ap*bm)-2*bm*bp*(k+1)^2*(k+2)^2;
        ak = 4*(k+1)*(ap*bm+am*bp)/DETk;
        bk = (-2*am*ap+(k^2+(k+1)^2)*(ap*bm-am*bp)+2*bm*bp*k^2*(k+1)^2)/DETk;
        Phi_x(k+1,:) = cos(k*theta)+ak*cos((k+1)*theta)+bk*cos((k+2)*theta);
    end;
end;
