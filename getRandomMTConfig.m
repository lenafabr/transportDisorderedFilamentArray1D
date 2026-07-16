function [x_coords, marker_positions] = getRandomMTConfig(MTparam_in)

% Function getEffectiveLandscape serves to calculate cargo probability density 
% and transition rates from a given MT configuration and kinetic parameters
%%Inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MTparam_in : structure containing MT-related parameters
%   .randomStart : scalar, x-position in MT start location
%   .randomEnd: scalar, x-position in MT end location
%   .snaptol: scalar, x-position MT ends snapping tolerance. Default = 0.1
%   .avglen: scalar, average length of each MT. 
%   .ncross : scalar, number of average MT cross-section
%   .rng: int, random number generator scalar
%%Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x_coords: Cell array of 2x1 matrices, for MT tips location. 
% marker_positions: Cell array of strings, 'right' and 'left' indicates MT plus ends
%   on the right and left hand side respectively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultParam = struct();
defaultParam.randomStart = 1;
defaultParam.randomEnd = 4;
defaultParam.snaptol = 0;

defaultParam.avglen = 1;% avg MT length 
defaultParam.ncross = 6; % typical number of MTs in cross-section
defaultParam.rng = 'shuffle';
defaultParam.constrain = true; % force exactly half of MTs to point right

param = struct();
% param.domain = []
param = copyStruct(MTparam_in,param,'addnew',true);


fields = fieldnames(defaultParam);
for i = 1:length(fields)
    field = fields{i};
    if ~isfield(param, field)
        % warning(['No user-defined "', field, '" provided. Using default value: ', num2str(defaultParam.(field)), '.']);
        param.(field) = defaultParam.(field);
    end
end
%%
rng(param.rng) %setup rng
randomScatteredMTLength = param.randomEnd - param.randomStart;

N = round((randomScatteredMTLength) * param.ncross/param.avglen); % total number of MTs
L = N*param.avglen; % total length of MTs

%% Distribute MT
[x_coords,marker_positions] = distributeMTUniform1D(N,L,randomScatteredMTLength,0.5,param.constrain);
for i = 1:length(x_coords)
    x_coords{i} = x_coords{i} + param.randomStart;
end
% x_coords = x_coords + param.randomStart;
% %% Distribute MT
% [x_coords_unmodified,marker_positions_unmodified] = distributeMTUniform1D(N,L,randomScatteredMTLength,0.5);
%% Add buffer such that the MT locations that intiially start at the origin will start at randomStart. 
% [x_coords, marker_positions] = addBufferRegion(x_coords_unmodified, marker_positions_unmodified, param);
end