%% Example runfile to test out functions to obtain analytic lanscape:
% Generates a random MT config, given MT parameters, and plot their energy
% landscape, distribution, effective diffusivity and velocity profile. 
%% Define Kinetic and MT-related parameters. 

% Below are the only parameter inputs required for energy landscape calculations
% User defined kinetic parameters
param = struct();
param.ka = 100; % Attachment rate in 1/s
param.kd = 100; % Detachment rate in 1/s
param.D = .1;   % Diffusivity in \mu m^2/s
param.v = 1;   % Processive average velocity in \mu m/s
param.spacing = 0.01; % Spacing between axial points in the domain. 

% Bias ranges from [-1,1]. 
% Bias = -1 -> all cargos move processively to the minus ends of MT tips. 
% Bias = 0 -> all cargos move processively to the minus ends of MT tips. 
% Bias = 1 -> all cargos move processively to the plus ends of MT tips. 
param.bias = 1;  

% User defined MT-related param ends snapping tolerance and domain length
MTparam = struct();
MTparam.domain = [0 5]; 
MTparam.snaptol = 0; % MT snapping tolerance, 0 means no snapping.  

%% Randomly scattered MT array
% Additional paramters when constructing randomly scattered MT configurations

MTparam.avglen = 1; % Average MT-length
MTparam.ncross = 2; % Average number of MT per cross-section
MTparam.rng = 80; % Comment out this line for random configuration.

% In addition to the domain, define where you would want the random MT configuration to begin and end within the domain.
MTparam.randomStart = 1; 
MTparam.randomEnd = 5; 
% Do you want the MT to be constrained to no net polarity (plus ends on right = plus ends on left of MT)?
MTparam.constrain = 'false';
[x_coords, marker_positions] = getRandomMTConfig(MTparam);
%% Define MTparam struct to include the xcoords and marker positions for the microtubular configuration
MTparam.x_coords = x_coords;
MTparam.marker_positions = marker_positions;

%% Plot MT config
figure;
% verbose = true;
plotparam = struct;
plotMT(MTparam.x_coords,MTparam.marker_positions,MTparam,plotparam)

%% Obtain Landscape parameters
[singleLayerInfo, singleContInfo] = getEffectiveLandscapeContinuous(param,MTparam);
%% Plotting Code
plotParam = struct();
plotParam.fontsize = 15;
plotParam.linspace = 100;
%%
figure
subplot(5,1,1)

plotMT(x_coords,marker_positions,MTparam,param)

title(sprintf('Bias = %g, D = %g, v = %g, ka = %g, kd = %g', ...
    param.bias, param.D, param.v, param.ka, param.kd), ...
    'Interpreter', 'latex', 'FontSize', plotParam.fontsize);

set(gcf,'color','w');
% set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',fontsize)
ax = gca;
hold(ax, 'on')
set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',plotParam.fontsize)
xlim([MTparam.domain(1) MTparam.domain(2)])

subplot(5,1,2)
hold on
plot(singleContInfo.bounds, singleContInfo.distrib, '-', 'MarkerSize', 20,'LineWidth',2);
% legend_entries = sprintf('Bias = %.2f', param.bias);
xlim([MTparam.domain(1) MTparam.domain(2)])

% legend(legend_entries, 'Location', 'Best','FontSize',10); 
ylabel('Cargo Density', 'Interpreter','latex','FontSize',plotParam.fontsize)
ax = gca;
hold(ax, 'on')
set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',plotParam.fontsize)

subplot(5,1,3)
% plot(NodeLocations, Energies, 'LineWidth', 2, 'LineStyle','-')
plot(singleContInfo.bounds, singleContInfo.energy, '-','LineWidth', 2)

xlim([MTparam.domain(1) MTparam.domain(2)])

ylabel('Energy ', 'Interpreter','latex','FontSize',plotParam.fontsize)

ax = gca;
hold(ax, 'on')
set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',plotParam.fontsize)


subplot(5,1,4)

stairs(singleLayerInfo.bounds, [singleLayerInfo.Deff;singleLayerInfo.Deff(end)],'LineWidth', 2)
yline(param.D, 'LineWidth', 2, 'LineStyle','--', 'Color','k');

set(gcf,'Color','w')
xlim([MTparam.domain(1) MTparam.domain(2)])
ylabel('Effective Diffusivity ', 'Interpreter','latex','FontSize',plotParam.fontsize)
ax = gca;
hold(ax, 'on')
set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',plotParam.fontsize)

subplot(5,1,5)

stairs(singleLayerInfo.bounds, [singleLayerInfo.Veff;singleLayerInfo.Veff(end)] ,'LineWidth', 2);
set(gcf,'Color','w')
xlim([MTparam.domain(1) MTparam.domain(2)])
ylabel('Effective Velocity ', 'Interpreter','latex','FontSize',plotParam.fontsize)
ax = gca;
hold(ax, 'on')
set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',plotParam.fontsize)
yline(0, 'LineWidth', 2, 'LineStyle','--', 'Color','k');


xlabel('Axial position', 'Interpreter','latex','FontSize',plotParam.fontsize)
