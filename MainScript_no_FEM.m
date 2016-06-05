%-------------------------------------------------------------------------%
%   2-D Spectral-Galerkin Battery Thermal Model (without FEM comparison)  %
%-------------------------------------------------------------------------%

% This code implements a Chebyshev Spectral Galerkin 2D thermal model.
% It does NOT compare the result against FEM simulations. Use this script
% if you do not have Matlab PDEtool installed on your computer. Otherwise,
% use the standard Mainscript.m

% Guide to structures used:
% 'p':      stores parameters and inputs (common to both models)
% 'sg':     stores spectral-Galerkin model inputs/outputs
% 'fem':    stores Finite-Element model inputs/outputs
% 'mult':   higher level struct storing both model outputs at multiple time instances

% An accompanying paper for this code entitled
% "On-board monitoring of 2-D spatially-resolved temperatures
% in cylindrical lithium-ion batteries: Part I. Low-order thermal modelling"
% has been submitted for publication in the Journal of Power Sources.

% Copyright (c) 2016 by Robert Richardson, Shi Zhao, David Howey
% and The Chancellor, Masters and Scholars of the University of Oxford.
% See the licence file LICENCE.txt for more information.

%% Initialise
clear; 
clc;
close all;
addpath(genpath('./Functions'));


%% Parameters and Inputs
p = struct;
p = func_parameters(p);                     % append parameters to struct p
p = func_inputs(p);                         % append inputs to struct p


%% Spectral-Galerkin (SG) Model
sg.N = 3;                                   % number of basis functions in each direction
sg.Ns = sg.N^2;                             % total number of states
[A,B,C,E,rr,zz,Psi,Te,ye] = ...
    func_spectral_galerkin(p, sg.N);        % calculate state matrices

% Simulate
sg.Psi = Psi;
sg.t_step = 1;
A_d = expm((E\A)*sg.t_step);
B_d = (E\A)\((A_d-eye(size(A_d)))*(E\B));

sg.x = zeros(length(p.t),length(A));
sg.x(1,:) = sg.Psi'\(p.Tinit-Te);
sg.u = zeros(size(B,2),length(p.t));
sg.y = zeros(length(p.t),size(C,1));
for i = 1:length(p.t)
    sg.u(:,i) = [p.Q(i)/p.Vb; 1];
    sg.x(i+1,:) = A_d*sg.x(i,:)' + B_d*sg.u(:,i);
    sg.y(i,:) = C*sg.x(i,:)';
end

% Assign
sg.Tz1 = sg.y(:,1) + ye(1);
sg.Tz2 = sg.y(:,2) + ye(2);
sg.Tr1 = sg.y(:,3) + ye(3);
sg.Tr2 = sg.y(:,4) + ye(4);



%% Distribution at multiple times
time_instances = [60,240,800];

for ii = 1:length(time_instances)
    ts = time_instances(ii); 

    % SG
    sg.rr = rr;
    sg.zz = zz;
    uu = sg.Psi'*sg.x(ts+1,:)';
    nfac = sqrt(size(sg.Psi,2));
    uu = reshape(uu,length(rr)/nfac,length(zz)/nfac);
    Te = reshape(Te,length(rr)/nfac,length(zz)/nfac);
    sg.ts.Tr = uu(round(size(uu,2)/2),:);
    sg.ts.Tr = sg.ts.Tr + Te(round(size(Te,2)/2),:);
    sg.ts.T = reshape(uu,length(rr)/nfac,length(zz)/nfac);
    sg.ts.T = sg.ts.T + reshape(Te,length(rr)/nfac,length(zz)/nfac);


    % Assign values to a structure containg multiple (3) time instances
    mult.sg{ii} = sg;
    mult.ts{ii} = ts;

end




%% Plots
% Plots the results for the SG model and FE model, corresponding to 
% plots (a)-(d) in Figs. 3 & 4 of article:
% (a) Temperature profiles: T(top, centre) and T(bottom, centre) vs. time
% (b) Heat generation profile: Q vs. t
% (c) Temperature distributions at selected times: T(r) at y = H/2
% (d) 2-D Temperature contour plots at selected time

% First, set plotting preferences:
set(0,'defaultaxeslinewidth',1.0,'defaultlinelinewidth',1.4)

% --------------- Plot (a): Temperature Evolution --------------- %
figure('position',[100 100 630 440])
hold on;
ha = plot([0, p.tfin], [p.Tinfr p.Tinfr], 'k-', 'linewidth',1);   % air
hw = plot([0, p.tfin], [p.Tinfb p.Tinfb], 'b--','linewidth',1);   % water
hSG = plot(p.t,sg.Tz1,'r--');                                     % SG    bottom, centre
plot(p.t,sg.Tz2,     'r--');                                      %       top, centre

% Axis limits/labels/legend
xlim([0 p.tfin])
ylims = [-10 80];
ylim(ylims)
xlabel('Time (s)');
ylabel('Temperature (^\circC)');
set(gca,'fontsize',11);
str2 = sprintf('SG (N = %i)',sg.Ns);
legend([hSG ha hw],str2,'Air','Coolant','location','best');
legend('boxoff')

% Plot dashed lines and labels at selected times
for i = 1:length(time_instances)
    ltick = 8;
    plot([mult.ts{i}, mult.ts{i}],[ylims(1), ylims(1)+ltick],'k-','linewidth',1)
    text(mult.ts{i}+10,ylims(1),sprintf('\\it{%0.2i s}', mult.ts{i}),'fontsize',11,'verticalalignment','bottom')
end
str_ar1 = '\it{T_{t,c}}';
str_ar2 = '\it{T_{b,c}}';
annotation('textarrow',[0.5270 0.4492],[0.5523 0.5091], 'String', str_ar1)
annotation('textarrow',[0.2841 0.3508],[0.4136 0.4636], 'String', str_ar2)




% --------------- Plot (b): Heat Generation Profile --------------- %
figure('position',[100 100 630 150])
plot(p.t, p.Q);
xlabel('Time (s)');
ylabel('Q (W/m^3)');
set(gca,'fontsize',11)




% --------------- Plot (c): Temperature Distributions --------------- %
figure('position',[100 100 600 180])
for i = 1:length(time_instances)
    subplot(1,3,i);
    
    % Plot
    hold on;
    plot(mult.sg{i}.rr(1:nfac:end),mult.sg{i}.ts.Tr, 'r--');

    % Axis limits, labels
    xlim([p.r1 p.r2])
    xlabel('r (m)');
    if i == 1; ylabel('Temperature (\circC)'); end
    set(gca,'fontsize',11);

    % Text
    ylims = get(gca,'ylim'); 
    text(p.r1+0.05*(p.r2-p.r1),ylims(1)+0.02*(ylims(2)-ylims(1)),sprintf('\\it{t = %0.2i s}', mult.ts{i}),'fontsize',11, 'verticalalignment','bottom')
    set(gca,'fontsize',11)

end





% --------------- Plot (d): 2-D Contour Plots --------------- %
% Preliminaries
figure('position',[100 100 329 385])

% Contour Plot 2 (SG)
hold on
contourf(sg.rr(1:nfac:end),sg.zz(1:nfac),mult.sg{2}.ts.T,10);
contourf(-sg.rr(1:nfac:end),sg.zz(1:nfac),mult.sg{2}.ts.T,10);
rectangle('Position',[-p.r1 0 2*p.r1 p.z2],'facecolor',[.8 .8 .8],'edgecolor','none')
str = sprintf('SG (N = %i)', sg.Ns);
title(str,'fontweight','n')
pbaspect([2.5*(p.r2-p.r1),p.z2,1])
set(gca,'ytick',[])
xlabel('r (m)');
colorbar
set(gca,'fontsize',11)

% Time String
str = sprintf('\\it{t = %0.2i s}', mult.ts{2});
annotation('textbox', [0.85,0.0,0.1,0.1],'String', str, 'linewidth', 1, 'fontsize', 11);








