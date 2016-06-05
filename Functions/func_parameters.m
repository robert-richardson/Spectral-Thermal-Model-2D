function [p] = func_parameters(p)
% func_parameters Define model parameters
% This function defines the parameters used in the SG and FEM models and
% stores them in struct 'p'.

% Copyright (c) 2016 by Robert Richardson, Shi Zhao, David Howey
% and The Chancellor, Masters and Scholars of the University of Oxford.
% See the licence file LICENCE.txt for more information.

% Parameters
p.r1 = 0.008/2;                     % inner radius [m]
p.r2 = 0.064/2;                     % outer radius [m]
p.z1 = 0;                           % bottom coordinate [m]
p.z2 = 0.198;                       % top coordinate [m]
p.Vb = pi*(p.r2^2-p.r1^2)*p.z2;     % volume [m^2]
p.rho = 2118;                       % density [kg m^-3]
p.cp = 795;                         % specific heat capacity [J kg^-1 K^-1]
p.kr = 0.666;                       % radial thermal conductivity [W m^-1 K^-1]
p.kz = 66.6;                        % axial thermal conductivity [W m^-1 K^-1]
p.hl = 0;                           % left convection coefficient [W m^-2]
p.hr = 30;                          % right convection coefficient [W m^-2]
p.hb = 400;                         % bottom convection coefficient [W m^-2]
p.ht = 30;                          % top convection coefficient [W m^-2]

end
    
    
    
    
    
    
    
    
