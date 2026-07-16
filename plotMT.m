function plotMT(x_coords,marker_positions,MTparam_in,plotParam_in, verbose)
% Function plotMT plots MT config

%%Inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% param_in : structure containing kinetic parameters
%   .bias : scalar, bias value
%   .D    : scalar, diffusion coefficient
%   .v    : scalar, velocity
%   .ka   : scalar, attachment rate
%   .kd   : scalar, detachment rate rate
%   .spacing : scalar, break down node locations further into integer
%   multiple of size spacing. Keep large number (not inf) to not implement
%   spacing. 
% MTparam_in : structure containing MT-related parameters
%   .snaptol : scalar, x-position MT ends snapping tolerance. Default = 0.1
%   .domain: 2x1 matrix of scalars. Entire system size. 
%   .x_coords: Cell array of 2x1 matrices, for MT tips location. 
%   .marker_positions: Cell array of strings, 'right' and 'left' indicates MT plus ends
%   on the right and left hand side respectively.
%%Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Veff_spaced_bounds:   1xn-1 matrix of scalars, effective velocities at
%   the regions between nodes.
% Deff_spaced_bounds : 1xn-1 matrix of scalars, effective diffusivity at
%   the regions between nodes.
% nodeLocations: nx1 matrix of scalars, x-locations of nodes. 
% energy_plot : 1xn matrix of scalars, energies at each node. 
% normalized_density : 1xn matrix of scalars, cargo probability density at each node. 
% kplus: 1xn matrix of scalars, transition rates to hop from one node to
%   a node to the right of it.
% kminus: 1xn matrix of scalars, transition rates to hop from one node to
%   a node to the left of it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
if nargin < 5
    verbose = false;
end
%%
defaultPlotParam = struct();
defaultPlotParam.fontsize = 20;

param = struct();
% param.domain = []
if (exist('plotParam_in','var'))
    param = copyStruct(plotParam_in,param,'addnew',true);
else
    param = defaultPlotParam;
end

fields = fieldnames(defaultPlotParam);
for i = 1:length(fields)
    field = fields{i};
    if ~isfield(param, field)
        if verbose == true
            warning(['No user-defined "', field, '" provided. Using default value: ', num2str(defaultPlotParam.(field)), '.']);
        end
        param.(field) = defaultPlotParam.(field);
    end
end
%%
defaultMTParam = struct();
defaultMTParam.snaptol = 0.1;
defaultMTParam.domain = NaN;

MTparam = struct();
MTparam = copyStruct(MTparam_in,MTparam,'addnew',true);

fields = fieldnames(defaultMTParam);
for i = 1:length(fields)
    field = fields{i};
    if ~isfield(MTparam, field)
        if verbose == true
        warning(['No user-defined "', field, '" provided. Using default value: ', num2str(defaultMTParam.(field)), '.']);
        end
        MTparam.(field) = defaultMTParam.(field);
    end
end
%%
[~,allbounds,MTpos,plusdir,mindir,~,~,~,~] = getConfigWallsModified_full(x_coords,marker_positions,MTparam.snaptol);
if (isnan(MTparam.domain))
    MTparam.domain = [min(allbounds) max(allbounds)];
end

% clf
plotMTposModified(MTpos,plusdir,mindir, allbounds)
ylim([0 length(plusdir)+1] )
xlim(MTparam.domain)
set(gcf,'color','w');
set(gca, 'YTick', []);
set(gca,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',param.fontsize)
% title(sprintf('MT Config, RNG = %g, Avglen = %g, ncross = %g', rng().Seed,avglen,ncross),'Interpreter', 'latex','FontSize',fontsize); %unmute for sweeps

end