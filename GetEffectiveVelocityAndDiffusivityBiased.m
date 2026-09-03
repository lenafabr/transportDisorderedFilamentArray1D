function [vels, Deff] = GetEffectiveVelocityAndDiffusivityBiased(D,v,ka,kd, fraction_kin, nMTdomain,nMTplus,nMTminus)
% Function GetEffectiveVelocityAndDiffusivityBiased calculates effective
%   velocity and diffusivities, incoorporated to consider bias. 
%%Inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D : scalar, Free diffusivity
% v : scalar, average processive velocity once attached
% ka : attachment rate
% kd : detachment rate. 
% fraction_kin: scalar, bias values that ranges from 0 - 1. 1-> all cargos
%   move to plus ends
% nMTdomain: n x 1 scalar, number of MT per domain
% MTdomain: 1 x n cell array, each cell contains the MT corresponding to
%   that region/domain.
% plusdir: n x 1 logic, 1 if second column of MTpos is a plus end
% mindir: n x 1 logic, 1 if first column of MTpos is a plus end
%%Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vels : 1 x n scalar, effective velocities in each region
% Deff : 1 x n scalar, effective diffusivities in each region.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make it all col vec;
nMTdomain = nMTdomain(:);
nMTplus = nMTplus(:);
nMTminus = nMTminus(:);

fraction_dyn = 1 - fraction_kin;

% nMTminus = nMTdomain-nMTplus; % # pointing left in each domain

fraction_MT_plus = nMTplus./nMTdomain;
fraction_MT_minus = nMTminus./nMTdomain;
alpha1 = fraction_kin * fraction_MT_plus + fraction_dyn * fraction_MT_minus;
alpha2 = fraction_kin * fraction_MT_minus + fraction_dyn * fraction_MT_plus;
dalpha = alpha1-alpha2;

%vels are the effective velocities.
vels = v*dalpha .* (nMTdomain*ka./(nMTdomain*ka+kd));

% Here Dvals are modified such that its D_eff from notes
Deff = (D * kd) ./ (nMTdomain*ka + kd) + (v^2 * nMTdomain*ka .* (-(nMTdomain*ka).^2 .* ((dalpha).^2 - 1) ...
    - 2 * nMTdomain*ka * kd .* ((dalpha).^2 - 1) + kd^2)) ./ (kd * (nMTdomain*ka + kd).^3);

ind0 = (nMTdomain==0);
Deff(ind0) = D;
vels(ind0) = 0;
%% Make it to row vectors
vels = vels';
Deff = Deff';

end