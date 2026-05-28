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

%% 
% User Defined MT configuration (Can comment out part A or B to test functionality)
%% A. Example case with x_coords (To test Case A mute block B, to test Case B mute Block A)
% Ensure that min(x_coords_input) >  MTparam.domain(1) & 
% and max(x_coords_input) <  MTparam.domain(2)

x_coords_input = {[1 2], [2 4], [4 5], [1.5,3.5]}; % Filament position
direction_input = {'left', 'right', 'left', 'right'}; % Which end of the filament  

% x_coords_input = {[1 5]};
% direction_input = {'right'};
%
% x_coords_input = {[1 2], [2 4],[2 4],[2 4],[2 4],[2 4], [4 5]};
% direction_input = {'left', 'right', 'left','right', 'left','right', 'left'};

%% Define MTparam struct to include the xcoords and marker positions for the microtubular configuration
MTparam.x_coords = x_coords_input;
MTparam.marker_positions = direction_input;
%% Set markers to calculate 
% param.calcIPR = 'True';
% param.calcIPRcut = 'False';
%% Plot MT config
% figure;
% % verbose = true;
% plotparam = struct;
% plotMT(MTparam.x_coords,MTparam.marker_positions,MTparam,plotparam)

%% Obtain Landscape parameters
[singleLayerInfo, singleContInfo] = getEffectiveLandscapeContinuous(param,MTparam);
%% Plotting Code

plotParam = struct();
plotParam.fontsize = 15;
plotParam.linspace = 100;

%%
MidNodeLocations = [];
for i = 1:length(NodeLocations) - 1
    MidNodeLocations(i) = (NodeLocations(i+1) + NodeLocations(i))/2;
end

XPoints = [];
VPoints = [];
DPoints = [];

for i = 1:length(NodeLocations)-1
    xSegment = linspace(NodeLocations(i), NodeLocations(i+1), plotParam.linspace);
    VSegment = repmat(Veff(i), size(xSegment));
    DSegment = repmat(Deff(i), size(xSegment));
    if i > 1
        xSegment = xSegment(2:end);
        VSegment = VSegment(2:end);
        DSegment = DSegment(2:end);
    end
    
    XPoints = [XPoints, xSegment];
    VPoints = [VPoints, VSegment];
    DPoints = [DPoints, DSegment];
end

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
plot(NodeLocations, CargoDistribution, '-', 'MarkerSize', 20,'LineWidth',2);
legend_entries = sprintf('Bias = %.2f', param.bias);
xlim([MTparam.domain(1) MTparam.domain(2)])

Reflect_points = [NodeLocations(1); NodeLocations(end)];
xline(Reflect_points, '--r','HandleVisibility','off');

legend(legend_entries, 'Location', 'Best','FontSize',10); 
ylabel('Cargo Density', 'Interpreter','latex','FontSize',plotParam.fontsize)
ax = gca;
hold(ax, 'on')
set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',plotParam.fontsize)

subplot(5,1,3)

plot(XPoints, VPoints,'LineWidth', 2, 'LineStyle','-')
% plot(XPoints, VPoints,'.-', 'MarkerSize',10)
xlim([MTparam.domain(1) MTparam.domain(2)])
ylabel('Effective V', 'Interpreter','latex','FontSize',plotParam.fontsize)
yline(0, '--r','HandleVisibility','off');

ax = gca;
hold(ax, 'on')
set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',plotParam.fontsize)


subplot(5,1,4)

plot(XPoints, DPoints, 'LineWidth', 2, 'LineStyle','-')
% plot(XPoints, DPoints,'.-', 'MarkerSize',10)

xlim([MTparam.domain(1) MTparam.domain(2)])
ylabel('Effective D', 'Interpreter','latex','FontSize',plotParam.fontsize)

% int_v = [];
% 
% for i = 1:length(spaced_bounds)-1
%     if i == 1
%         int_v(1) = 0 ;
%     else
%         int_v(i) = int_v(i-1) + Veff_spaced_bounds(i)*(spaced_bounds(i+1)-spaced_bounds(i));
%     end
% end
% plot(MidNodeLocations, int_v, 'LineWidth', 2, 'LineStyle','-')
% xlim([0 max(allbounds)])
% ylabel('Integral v', 'Interpreter','latex','FontSize',fontsize)


ax = gca;
hold(ax, 'on')
set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',plotParam.fontsize)

subplot(5,1,5)

% plot(NodeLocations, Energies, 'LineWidth', 2, 'LineStyle','-')
plot(NodeLocations, Energies, '.-')

xlim([MTparam.domain(1) MTparam.domain(2)])

ylabel('Energy ', 'Interpreter','latex','FontSize',plotParam.fontsize)

ax = gca;
hold(ax, 'on')
set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',plotParam.fontsize)

xlabel('x positions', 'Interpreter','latex','FontSize',plotParam.fontsize)

%%
%% Calculating Eigenvalues for relaxation time
[sortedEigVal,sortedEigVec] = getRelaxTime(NodeLocations,kplus, kminus);
Relax_time = -1./sortedEigVal;

initialDist = ones(length(NodeLocations),1)/length(NodeLocations);
coeff = sortedEigVec\initialDist;

tmax = 10000;     
timeVec = linspace(0, tmax, 2000);  

figure;
% hold on
for tIndex = 1:length(timeVec)
    tCurrent = timeVec(tIndex);

    p_t = zeros(length(NodeLocations),1);
    for mode = 1:length(NodeLocations)
        p_t = p_t + coeff(mode) * sortedEigVec(:, mode) * exp(sortedEigVal(mode) * tCurrent);
    end
    p_t = p_t/trapz(NodeLocations,p_t);

    plot(NodeLocations, p_t);
    xlabel('xpos');
    ylabel('Cargo Probability Density');
    title(sprintf('Probability Distribution at t = %0.1f s', tCurrent));
    % ylim([0, 1]);
    
    % Pause briefly to create the animation effect
    drawnow;
    pause(0.1);  % adjust pause duration if needed
end