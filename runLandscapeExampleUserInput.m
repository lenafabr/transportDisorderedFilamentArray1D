%% Example runfile to test out functions to obtain analytic lanscape:
% Not done
%% 
% Generates cargo distribution profile, with energy landscape.
%% Define Kinetic and MT-related parameters. 

% Below are the only parameter inputs required for energy landscape calculations
% User defined kinetic parameters
param = struct();
param.ka = 100; % Attachment rate in 1/s
param.kd = 100; % Detachment rate in 1/s
param.D = .1;   % Diffusivity in \mu m^2/s
param.v = 1;   % Processive average velocity in \mu m/s
param.spacing = 0.01; % Spacing between axial points in the domain. 

% Bias ranges from 0 - 1. 
% Bias = 0 -> all cargos move processively to the minus ends of MT tips. 
% Bias = 1 -> all cargos move processively to the plus ends of MT tips. 
param.bias = 1;  

% User defined MT-related param ends snapping tolerance and domain length
MTparam = struct();
MTparam.domain = [0 6]; 
MTparam.snaptol = 0; % MT snapping tolerance, 0 means no snapping.  

%% User Defined MT configuration
MTparam.x_coords = {[1 2], [2 4], [4 5], [1.5,3.5]}; % Filament position
MTparam.marker_positions = {'left', 'right', 'left', 'right'}; % Which end of the filament  
%% Calculate for the effective number of cross-section. 
MTlen = cellfun(@(x) abs(x(2) - x(1)), MTparam.x_coords);
totalMTlen = sum(MTlen); % Total filament length
MTparam.ncross = totalMTlen/(MTparam.domain(2) -MTparam.domain(1)); 
%% Plot MT config
figure;
plotparam = struct;
plotMT(MTparam.x_coords,MTparam.marker_positions,MTparam,plotparam)
%% Obtain Landscape parameters
[singleLayerInfo, singleContInfo] = getEffectiveLandscapeContinuous(param,MTparam);
%% Plotting Code
plotParam = struct();
plotParam.fontsize = 15;
%%
figure
% Plotting filament configuration. 
subplot(5,1,1)
plotMT(MTparam.x_coords,MTparam.marker_positions,MTparam,param)

title(sprintf('Bias = %g, D = %g, v = %g, ka = %g, kd = %g', ...
    param.bias, param.D, param.v, param.ka, param.kd), ...
    'Interpreter', 'latex', 'FontSize', plotParam.fontsize);

set(gcf,'color','w');
% set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',fontsize)
ax = gca;
hold(ax, 'on')
set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',plotParam.fontsize)
xlim([MTparam.domain(1) MTparam.domain(2)])

% Ploting cargo densities
subplot(5,1,2)
hold on
plot(singleContInfo.bounds, singleContInfo.distrib, '-', 'MarkerSize', 20,'LineWidth',2);
xlim([MTparam.domain(1) MTparam.domain(2)])

ylabel('$\rho(x)$', 'Interpreter','latex','FontSize',plotParam.fontsize)
ax = gca;
hold(ax, 'on')
set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',plotParam.fontsize)

% Plotting energy landscape
subplot(5,1,3)
% plot(NodeLocations, Energies, 'LineWidth', 2, 'LineStyle','-')
plot(singleContInfo.bounds, singleContInfo.energy, '-','LineWidth', 2)

xlim([MTparam.domain(1) MTparam.domain(2)])

ylabel('$\epsilon (x)$ ', 'Interpreter','latex','FontSize',plotParam.fontsize)

ax = gca;
hold(ax, 'on')
set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',plotParam.fontsize)

%Plotting Effective Diffusivity
subplot(5,1,4)

stairs(singleLayerInfo.bounds, [singleLayerInfo.Deff;singleLayerInfo.Deff(end)],'LineWidth', 2)
yline(param.D, 'LineWidth', 2, 'LineStyle','--', 'Color','k');

xlim([MTparam.domain(1) MTparam.domain(2)])
ylabel('$D_{\mathrm{eff}}$ ', 'Interpreter','latex','FontSize',plotParam.fontsize)
ax = gca;
hold(ax, 'on')
set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',plotParam.fontsize)

% Plotting Effective Velocity
subplot(5,1,5)
stairs(singleLayerInfo.bounds, [singleLayerInfo.Veff;singleLayerInfo.Veff(end)] ,'LineWidth', 2);
xlim([MTparam.domain(1) MTparam.domain(2)])
ylabel('$v_{\mathrm{eff}}$', 'Interpreter','latex','FontSize',plotParam.fontsize)
ax = gca;
hold(ax, 'on')
set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',plotParam.fontsize)
yline(0, 'LineWidth', 2, 'LineStyle','--', 'Color','k');


xlabel('Axial position', 'Interpreter','latex','FontSize',plotParam.fontsize)
