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

function [p] = func_parameters(p)
% func_parameters Define model parameters
% This function defines the parameters used in the SG and FEM models and
% stores them in struct 'p'.

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
    
    
    
    
    
    
    
    
