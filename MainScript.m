%-------------------------------------------------------------------------%
%               2-D Spectral-Galerkin Battery Thermal Model               %
%-------------------------------------------------------------------------%

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

% This code implements a Chebyshev Spectral Galerkin 2D thermal model and
% compares it against a finite element (FEM) model using Matlab PDEtool.
% Execution of the Mainscript.m file runs the simulation.

% Guide to structures used:
% 'p':      stores parameters and inputs (common to both models)
% 'sg':     stores spectral-Galerkin model inputs/outputs
% 'fem':    stores Finite-Element model inputs/outputs
% 'mult':   higher level struct storing both model outputs at multiple time instances


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



%% Finite Element (FEM) Model 
% ...uses Matlab PDEtool

% --- Run pde tool --- %
fem.ninc = 2000;
fem.t_step = p.tfin/fem.ninc;
[tlist, u, F, nels] = func_pdetool(p, fem.t_step);

% --- Parameters at all times--- %
nuP = 101;
nuPm = round(nuP/2);
fem.xP = linspace(p.r1,p.r2,nuP);
fem.uPr = evaluate(F, fem.xP,(p.z2/2)*ones(1,nuP));   % Evaluate pde interpolant
fem.yP = linspace(0,p.z2,nuP);
for ii = 1:length(fem.yP)
    fem.uP(ii,:,:) = evaluate(F, fem.xP, repmat(fem.yP(ii),1,length(fem.xP)));
end
fem.t = linspace(1,p.tfin,length(fem.uPr(1,2:end)));
fem.Tmax = max(fem.uPr(:,2:end));


% Assign model outputs
fem.Tz1(1,:) = fem.uP(1, nuPm,2:end);
fem.Tz2(1,:) = fem.uP(end, nuPm,2:end);
fem.Tr1(1,:) = fem.uP(nuPm, 1,2:end);
fem.Tr2(1,:) = fem.uP(nuPm, end,2:end);



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


    % FEM
    fem.ts.TPr = fem.uPr(:,round(fem.ninc*ts/p.tfin)+1);

    % Assign values to a structure containg multiple (3) time instances
    mult.fem{ii} = fem;
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
hFEM = plot(fem.t,fem.Tz1,   'k-');                               % FEM   bottom, centre
plot(fem.t,fem.Tz2,    'k-');                                     %       top, centre
hSG = plot(p.t,sg.Tz1,'r--');                                     % SG    bottom, centre
plot(p.t,sg.Tz2,     'r--');                                      %       top, centre

% Axis limits/labels/legend
xlim([0 p.tfin])
ylims = [-10 80];
ylim(ylims)
xlabel('Time (s)');
ylabel('Temperature (^\circC)');
set(gca,'fontsize',11);
str1 = sprintf('FEM (N = %i)',nels);
str2 = sprintf('SG (N = %i)',sg.Ns);
legend([hFEM hSG ha hw],str1,str2,'Air','Coolant','location','best');
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
    plot(fem.xP, mult.fem{i}.ts.TPr,'k');
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
figure('position',[100 100 630 385])
pdecontour = fem.uP(:,:,round(time_instances(2))/fem.t_step+1);
Tlim1 = min([mult.sg{2}.ts.T(:); pdecontour(:)]);                   % identify min. temperature
Tlim2 = max([mult.sg{2}.ts.T(:); pdecontour(:)]);                   % identify max. temperature

% Contour Plot 1 (FEM)
subplot(1,2,1);
hold on
contourf(fem.xP,fem.yP,pdecontour,10);                                                  % right contour
contourf(-fem.xP,fem.yP,pdecontour,10);                                                 % left contour
rectangle('Position',[-p.r1 0 2*p.r1 p.z2],'facecolor',[.8 .8 .8],'edgecolor','none')   % core mandrel
str = sprintf('FEM (N = %i)',nels);
title(str,'fontweight','n','fontsize',11);
pbaspect([2.5*(p.r2-p.r1),p.z2,1])
caxis([Tlim1 Tlim2])
xlabel('r (m)');
ylabel('z (m)');
set(gca,'fontsize',11)

% Contour Plot 2 (SG)
subplot(1,2,2);
hold on
contourf(sg.rr(1:nfac:end),sg.zz(1:nfac),mult.sg{2}.ts.T,10);
contourf(-sg.rr(1:nfac:end),sg.zz(1:nfac),mult.sg{2}.ts.T,10);
rectangle('Position',[-p.r1 0 2*p.r1 p.z2],'facecolor',[.8 .8 .8],'edgecolor','none')
str = sprintf('SG (N = %i)', sg.Ns);
title(str,'fontweight','n')
pbaspect([2.5*(p.r2-p.r1),p.z2,1])
caxis([Tlim1 Tlim2])
set(gca,'ytick',[])
xlabel('r (m)');
colorbar
set(gca,'fontsize',11)

% Time String
str = sprintf('\\it{t = %0.2i s}', mult.ts{2});
annotation('textbox', [0.85,0.0,0.1,0.1],'String', str, 'linewidth', 1, 'fontsize', 11);








