function [singleLayerInfo, singleContInfo] =getEffectiveLandscapeContinuous(param_in,MTparam_in)
% Function getEffectiveLandscapeContinuous takes in kinetic and geometric
% MT parameters to provide information of the cargo energy landscape,
% distribution and inverse IPR, which are:
% evaluated at region boundaries (singleLayerInfo)
% evaluated at x-position values with spacing param.spacing (singleContInfo)
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
%   .snaptol : scalar, x-position MT ends snapping tolerance. Default = 0 (
%   no snapping)
%   .domain: 2x1 matrix of scalars. Entire system size. 
%   .x_coords: Cell array of 2x1 matrices, for MT tips location. 
%   .marker_positions: Cell array of strings, 'right' and 'left' indicates MT plus ends
%   on the right and left hand side respectively.
%%Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% singleLayerInfo : Evaluated the landscape at region boundaries
%   .bounds   : col vec (m,1), x-positions of region boundaries (MT endpoints)
%   .energy   : col vec (m,1), energy evaluated at bounds.
%   .Veff     : col vec (m-1, 1), effective velocities at regions (between boundaries).  
%   .Deff     : col vec (m-1, 1), effective diffusivities at regions.  
%   .n        : col vec (m-1, 1), MT cross-sectional number at regions.  
% singleContInfo  : Evaluated the landscape at x-position values spaced by param.spacing
%   .bounds   : col vec (p,1), x-positions of region boundaries (MT endpoints)
%   .energy   : col vec (p,1), energy evaluated at xvals.
%   .distrib  : col vec (p,1), energy evaluated at bounds.
%   .n        : col vec (p,1), energy evaluated at bounds.
%   .regbounds: col vec (m,1), x-positions of region boundaries (MT endpoints) 
%   .Veff     : col vec (m-1, 1), effective velocities at regions (between boundaries).
%   .Deff     : col vec (m-1, 1), effective diffusivities at regions.
%   .invIPR   : scalar, inverse IPR, provides information on the support of
%   the probability distribution.
%   .invIPRcut: scalar, inverse IPR, evaluated for the probability distribution 
%   at xvals 1 MT away from the domains. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defining our input kinetic and MT-parameters 
defaultParam = struct();
defaultParam.ka = 100;
defaultParam.kd = 100;
defaultParam.D = 0.1;
defaultParam.v = 1;

% default to big numbers for spacing -> No spacing. 
defaultParam.spacing = 0.1;
defaultParam.bias = 1;
defaultParam.calcIPR = 'True';
defaultParam.calcIPRcut = 'True';
% Copy structure of the inputs
param = struct();
param = copyStruct(param_in,param,'addnew',true);

% Check if the input fieldnames missing that does not exist in the defaults
fields = fieldnames(defaultParam);
for i = 1:length(fields)
    field = fields{i};
    if ~isfield(param, field)
%         warning(['No user-defined "', field, '" provided. Using default value: ', num2str(defaultParam.(field)), '.']);
        param.(field) = defaultParam.(field);
    end
end
%% Redefine Bias parameter, to go from -1 -> 1 to 0 -> 1;
param.bias = (param.bias + 1)/2;
%% 
defaultMTParam = struct();
defaultMTParam.snaptol = 0;

% Copy structure of the inputs
MTparam = struct();
MTparam = copyStruct(MTparam_in,MTparam,'addnew',true);

% Check if the input fieldnames missing that does not exist in the defaults
fields = fieldnames(defaultMTParam);
for i = 1:length(fields)
    field = fields{i};
    if ~isfield(MTparam, field)
        warning(['No user-defined "', field, '" provided. Using default value: ', num2str(defaultMTParam.(field)), '.']);
        MTparam.(field) = defaultMTParam.(field);
    end
end
%% 
x_coords = MTparam.x_coords;
marker_positions = MTparam.marker_positions;

% Code will not function when spacing is set to inf or very large values. 
if (param.spacing == inf) || (param.spacing > 1e5)
    error('Large spacing, please double check.')
end
%% Check for if x_coords or marker positions are defined. 
if ~exist('x_coords', 'var') || isempty(x_coords)
    error('No x_coords being passed into the function, please doluble check')
end
if  ~exist('marker_positions', 'var') || isempty(marker_positions)
    error('No Marker positions being passed into the function, please doluble check')
end
%% Check for User input error in x_coords
if length(marker_positions) > length (x_coords)
    warning('Marker positions has more inputs than x_coords plese double check your inputs')
elseif length(marker_positions) < length (x_coords)
    warning('Marker positions has less inputs than x_coords plese double check your inputs')
end
if any(cellfun(@(x) x(1) > x(2), x_coords))
    error('Please double check your notation for the x_coords. First number within each cell should be smaller than the second one')
end
%% Obtain Relevant Filament information
[~,allbounds,MTpos,plusdir,mindir,nMTdomain,nMTplus,nMTminus] = getConfigWallsModified(x_coords,marker_positions,MTparam.snaptol);
%% Check for User input error

if isempty(MTparam.domain)
    error('No MT domain defined.')
end
if MTparam.domain(1) >= MTparam.domain(2)
    error('Syntax for MTparam.domain is the first input is smaller than the second one.')
end
if min(allbounds) <= MTparam.domain(1)
    error('For energy landscape calculations please have the minimum MT end be greater than the minimum of the MT Domain')
elseif max(allbounds) >= MTparam.domain(2)
    error('For energy landscape calculations please have the maximum MT end be less than the maximum of the MT Domain')
end

%% Obtain effective velocities and effective diffusivities
%[vels, Deff] = GetEffectiveVelocityAndDiffusivityBiased(param.D,param.v,param.ka,param.kd, param.bias, nMTdomain,MTdomain, plusdir, mindir);
[vels, Deff] = GetEffectiveVelocityAndDiffusivityBiased(param.D,param.v,param.ka,param.kd, param.bias, nMTdomain,nMTplus,nMTminus);
%% Account for the 0 MT regions at domain ends. 
allbounds = [MTparam_in.domain(1);allbounds; MTparam_in.domain(2)];
nMTdomain = [0; nMTdomain; 0];
vels = [0, vels, 0];
Deff = [param.D, Deff, param.D];
%% Construct the axial position values, given the spacing input. 
npt = ceil(MTparam.domain(end)-MTparam.domain(1))/param.spacing;
xvals = linspace(MTparam.domain(1),MTparam.domain(end),npt);
%% Construct Effective kinetic parameters Deff and Veff for uniformly spaced nodes 
% and evaluate the energy landscape + cargo distribution. 
[energy,normProbDistrib,IPR,IPRcut,energybounds,numMT] = getLandscapeContinuous(param,MTparam.ncross, Deff, vels,allbounds,nMTdomain,xvals,MTpos);

%% change everything to column vectors (X-())
vels = vels(:);
Deff = Deff(:);
energy = energy(:);
normProbDistrib = normProbDistrib(:);

%% Outputs are column vectors. 

% Info contains regions information 
singleLayerInfo = struct();
singleLayerInfo.Veff = vels;
singleLayerInfo.Deff = Deff;
singleLayerInfo.n = nMTdomain;
singleLayerInfo.bounds = allbounds;
%singleLayerInfo.distrib = normProbDistrib;
singleLayerInfo.energy = energybounds(:);
%singleLayerInfo.invIPR = invIPR;

% Cont contains contiuous values
singleContInfo = struct();
singleContInfo.Veff = vels;
singleContInfo.Deff = Deff;
singleContInfo.n = numMT(:);
singleContInfo.bounds = xvals(:);
singleContInfo.regbounds = allbounds;
singleContInfo.distrib = normProbDistrib;
singleContInfo.energy = energy(:);
singleContInfo.invIPR = 1/IPR;
singleContInfo.invIPRcut = 1/IPRcut;
end