%% Hybrid-Gillespie Simulations.
% clear all
%%
% Construct Random MT configuration.
MTparam = struct();
MTparam.avglen = 1; % Average MT-length
MTparam.ncross = 7; % Average number of MT per cross-section
MTparam.domain = [1 7]; % Domain start and end points

% MT construction RNG. 
MTparam.rng = 379; % Example configuration
% MTparam.rng = 'shuffle'; % Random Configuration

% In addition to the domain, define where you would want the random MT 
% configuration to begin and end within the domain.
MTparam.randomStart = 1;
MTparam.randomEnd = 7;

% Construct random MT configuration
[x_coords, marker_positions] = getRandomMTConfig(MTparam);

MTparam.x_coords = x_coords;
MTparam.marker_positions = marker_positions;
MTparam.snaptol = 0.1; % MT snapping tolerance, 0 means no snapping.  

param = struct();
% attachment and detachment rate (assumed independent of direction)
param.kd = 100;
param.ka = 100;
param.D = 0.1; % diffusivity
param.v = 1; % walking speed

% -----------
% parameters specifically for the stochastic simulation
% -----------
param.npart = 1000; % number of paraticles
param.maxevent = 1e10; % max number of events for each particle

% scaling prefactor forhow much smaller than the event times to make
% the discrete diffusive timestep
param.dtscale = 1E-2;
% param.dtscale = 2;

% save times for particle positions
param.savetimes = linspace(0,1e5,1e4);
% param.savetimes = linspace(0,1e4,1e2);

% how much info to print out
param.verbose = 0; % 0 provides summary, else, provides detailed state and event times
param.printevery = 1e4; % event interval for printing simulation details

tic
pgroup = localizationSims(param,MTparam);
toc
%% Plot simulation movie
figure
fontsize = 20;
plotparam = struct;
[~,allbounds,MTpos,plusdir,mindir,nMTdomain,MTdomain,MTstartbound,MTendbound] = ...
    getConfigWallsModified_full(x_coords,marker_positions,MTparam.snaptol,MTparam.domain);
plotMT(MTparam.x_coords,MTparam.marker_positions,MTparam,plotparam)

hold all
colormap(jet)
pshow = 1:2:pgroup.npart; % which particles to show
[pos0,ind] = sort(pgroup.savepos(pshow,1)); % sort particles by initial position
phandle = scatter(pos0,pgroup.savestate(pshow(ind),1,2),50,1:length(pshow),'filled');

% Show Diffusive on the top
phandle.YData(ind) = length(plusdir) + 2;
xlim([MTparam.domain(1) MTparam.domain(2)])
set(gcf,'color','w');
set(gca,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',fontsize)
xlabel('Domain Length ($\mu$m)','Interpreter', 'latex','FontSize',fontsize);
set(gca, 'YTick', []);

yl = ylim;
ylim([0 yl(2)+3])
title(sprintf('time = %g',param.savetimes(1)))
pause(1);

% Plot the rest of the movie according to the savetimes. 
for tc = 1:length(param.savetimes)
    phandle.XData = pgroup.savepos(pshow(ind),tc);
    diff_ind = find(pgroup.savestate(pshow(ind),tc,2) == 0);
    phandle.YData = pgroup.savestate(pshow(ind),tc,2);
    phandle.YData(diff_ind) = length(plusdir) + 2;
    pause(0.1);
    title(sprintf('time = %g',param.savetimes(tc)))
    drawnow
end


hold off

%% Histogramming particle distribution at the end of the simulation
figure;
XData = pgroup.savepos(:,end);
histogram(XData,'BinWidth', 0.1, 'Normalization','pdf',...
         'FaceColor',[0.6, 0.6, 0.6],'DisplayName','new sims');
set(gcf,'color','w');
set(gca,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',fontsize)
