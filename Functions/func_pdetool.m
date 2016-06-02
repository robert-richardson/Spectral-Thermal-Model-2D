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

function [tlist, u, F, nels] = func_pdetool(p, tinc)
% This functions implements a finite element simulation of the model using
% Matlab PDE-Tool. Matlab pdetool must be installed for this to run.


%% Inputs
Q1 = p.q1/p.Vb; t1 = p.t1;
Q2 = p.q2/p.Vb; t2 = p.t2;
Q3 = p.q3/p.Vb;


%% Definition of PDE Coefficients
c = char(sprintf('%g*x', p.kr), sprintf('%g*x',p.kz));    
d = sprintf('%g*x', p.rho*p.cp);
f = sprintf('x*( %g*heaviside(%g-t) + %g*heaviside(t-%g)*heaviside((%g+%g)-t) + %g*heaviside(t-(%g+%g)) )',...
    Q1,t1,Q2,t1,t1,t2,Q3,t1,t2);
a = 0;
% c = [p.kr, p.kz]'; d = p.rho*p.cp; f = q; a = 0;


%% Geometry and Mesh
g = decsg([3 4 p.r1 p.r2 p.r2 p.r1 p.z1 p.z1 p.z2 p.z2]');
pg = pdeGeometryFromEdges(g);
hmax = .02*(p.r2-p.r1+p.z2-p.z1)/2;
[pmesh, e, t] = initmesh(g, 'Hmax', hmax);
nels = size(t,2);


%% Boundary Conditions
numberOfPDE = 1;
pb = pde(numberOfPDE);
leftQ  = @(thePde, loc, state) p.hl*loc.x;    leftG = @(thePde, loc, state) p.hl*loc.x*p.Tinfl;
rightQ  = @(thePde, loc, state) p.hr*loc.x;    rightG = @(thePde, loc, state) p.hr*loc.x*p.Tinfr;
topQ  = @(thePde, loc, state) p.ht*loc.x;    topG = @(thePde, loc, state) p.ht*loc.x*p.Tinft;
botQ  = @(thePde, loc, state) p.hb*loc.x;    botG = @(thePde, loc, state) p.hb*loc.x*p.Tinfb;


bLeft   = pdeBoundaryConditions(pg.Edges(4), 'g', leftG, 'q', leftQ);
bRight  = pdeBoundaryConditions(pg.Edges(2), 'g', rightG, 'q', rightQ);
bTop  = pdeBoundaryConditions(pg.Edges(3), 'g', topG, 'q', topQ);
bBot  = pdeBoundaryConditions(pg.Edges(1), 'g', botG, 'q', botQ);
pb.BoundaryConditions = [bRight bLeft bTop bBot];



%% Run
tlist = 0:tinc:p.tfin;
u0 = p.Tinit;                                       % initial T
u = parabolic(u0, tlist, pb, pmesh,e,t,c,a,f,d);    % Solve
F = pdeInterpolant(pmesh,t,u);                      % PDE Interpolant






