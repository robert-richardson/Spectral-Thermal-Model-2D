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

function [p] = func_inputs(p)
%FUNC_INPUTS Defines model inputs
% This function defines the inputs to the SG and FEM models and
% stores them in struct 'p'.

% External temperatures at each edge
p.Tinfl = 18;                               % left (air) temperature [K]
p.Tinfr = 18;                               % right (air) temperature [K]
p.Tinfb = 3;                                % bottom (water) temperature [K]
p.Tinft = 18;                               % top (air) temperature [K]
p.Tinit = 18;                               % initial (cell) temperature [K]


% Heat generation pulse power profile
p.tfin = 2000;                              % end time [s]
p.t = 0:1:p.tfin;                           % total time span  [s]
p.q1 = 50;                                  % STAGE 1 heat gen. rate [W m^-3]
p.t1 = 150;                                 % STAGE 1 end time [s]
p.q2 = 1000;                                % STAGE 2 heat gen. rate [W m^-3]
p.t2 = 50;                                  % STAGE 2 end time [s]
p.q3 = 0;                                   % STAGE 3 heat gen. rate [W m^-3]
        
p.Q = zeros(length(p.t),1);                 % initialise
for i = 1:length(p.t);
  p.Q(i) = p.q1*heaviside(p.t1-p.t(i)) ...  % STAGE 1
    + p.q2*heaviside(p.t(i)-p.t1)*...       % STAGE 2
    heaviside((p.t1+p.t2)-p.t(i)) ...       % STAGE 2 (cont.)
    + p.q3*heaviside(p.t(i)-(p.t1+p.t2));   % STAGE 3
end

end

