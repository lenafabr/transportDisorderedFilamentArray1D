function [energy,probdistrib,IPR,IPRcut,energybounds,numMT] = ...
    getLandscapeContinuous(param, ncross,Deff, vels,allbounds,nMTdomain,xvals,MTpos)
% Takes kinetic and geometric parameters, and solves for the energy
% landscape, probability distribution and calculates for the IPR.

% Used at runLandscapeExampleRandConfig.m

%%%%%%%%%%%%%%
% Everything assumed to be row vectors here!
% Input:
% param: kinetic parameters
% ncross: average number of MT per cross-section. 
% Deff, vels: effective diffusivity and velocity in each region
% allbounds: boundary points (including domain edges), in order
% nMTdomain: number of MTs in each domain
% xvals: axial positions at which to evaluate landscape
% MTpos: x-positions of MT start and end points, sorted by position.
% if xvals are out of bounds, the energy will default to the value at the
% domain ends

% Output:
% energy = energy values at all the xvals
% probdistrib = normalized prob distribution at each x value
% IPR = inverse participation ratio
% IPRcut = inverse participation ratio with regions corresponding to first
% and last MTs left off.
% energybounds = energy just before each boundary (for testing purposes)
% numMT = number of MTs corresponding to each x value
%%%%%%%%%%%%%%

if (size(allbounds,1)>size(allbounds,2)) % make into row vector
    allbounds = allbounds';
end
if (size(nMTdomain,1)>size(nMTdomain,2)) % make into row vector
    nMTdomain = nMTdomain';
end
if (size(xvals,1)>size(xvals,2)) % make into row vector
    xvals = xvals';
end

% get region lengths:
reglen = diff(allbounds);

% get the Pe numbers
Pe = vels.*reglen./Deff;
% cumulative Pe
cumPe = [0,cumsum(Pe)];

% get discontinuity coefficients across region boundaris
z1 = 1 + (param.ka/param.kd)*nMTdomain;
z = (param.ka*nMTdomain + param.kd)/(param.ka*ncross + param.kd);
entropicTerm = log(z);

% For each position, find which region it belongs in
nbound = length(allbounds);
regind = interp1(allbounds,1:nbound,xvals,'previous');
regind(xvals<=allbounds(1)) = 1; % deal with value directly on domain edge;
regind(xvals>=allbounds(end)) = nbound-1; % deal with value directly on domain edge;

numMT = nMTdomain(regind);

% x values relative to start of region
xrel = xvals - allbounds(regind);

% energy at all xvals
energy = -cumPe(regind) - xrel.*vels(regind)./Deff(regind) - entropicTerm(regind);
% energies at region boundaries specifically (right before each boundary)
energybounds = -cumPe - [0,entropicTerm];

%% Solving for the normalization constant
% % Original:
% % integral of exp(energy) over each region
% cumPe_reg = cumPe(1:end-1);
% distint = Deff./vels.*(exp(Pe) - 1).*exp(cumPe_reg).*z;
% ind0 = (vels==0);
% distint(ind0) = reglen(ind0).*exp(cumPe_reg(ind0)).*z(ind0);
% 
% % normalization constant
% normConst = sum(distint);

% Large SumPe Case:
% integral of exp(energy) over each region
cumPe_reg = cumPe(1:end-1);
logdistint = log(abs(Deff./vels)) + log(abs(exp(Pe) - 1)) + cumPe_reg + log(z);
ind0 = (vels==0);
logdistint(ind0) = log(reglen(ind0)) + cumPe_reg(ind0) + log(z(ind0));
logdistintMax = max(logdistint);
logNormConst = log(sum(exp(logdistint - logdistintMax))) + logdistintMax;
% logIPR = logdistint2 - 2*log(normConst);

% normalization constant
normConst = exp(logNormConst);

logprobdistrib = -energy - logNormConst;
probdistrib = exp(logprobdistrib);
%% probability distribution
% % Original normalized distribution (issues when cumsumPe gets big)
% probdistrib = exp(-energy)/normConst;

%% Original IPR calculations, but will explode for large cumulative Pe. 
% integrate the distribution squared
% Original IPR calculations (issues when cumsumPe gets big)
% distint2 = Deff./(2*vels).*(exp(2*Pe) - 1).*exp(2*cumPe_reg).*z.^2;
% distint2(ind0) = reglen(ind0).*exp(2*cumPe(ind0)).*z(ind0).^2;

% IPR = sum(distint2)/normConst^2;
%% Large cumsum Pe case:
% Never encountered a case of very large Pe, so not used for now. 
% specifically for this term: log(abs(exp(2*Pe) - 1)) used below
% largePeVal = 2*Pe > 100;
% normalPe = ~largePeVal;
% log2Pe(largePeVal) = 2*Pe(largePeVal);
% log2Pe(normalPe) = log(abs(exp(2*Pe(normalPe)) - 1));

% Take the log of distinct2.
logdistint2 = log(abs(Deff./(2*vels))) + log(abs(exp(2*Pe) - 1)) + 2*cumPe_reg + 2*log(z);
logdistint2(ind0) = log(reglen(ind0)) + 2*cumPe_reg(ind0) + 2*log(z(ind0));

logdistint2Max = max(logdistint2);
% remove e^max, take sum, then add back the factorized e^max, then log back 
logsumdistint2 = log(sum(exp(logdistint2 - logdistint2Max))) + logdistint2Max;

logIPR = logsumdistint2 - 2*logNormConst;
% logIPR = log(sum(distint2))- 2*log(normConst);
IPR = exp(logIPR);


% cut off the regions corresponding to a single MT length on either size
cutind = interp1(allbounds,1:nbound,[MTpos(1,2),MTpos(end,1)],'nearest');
useind = cutind(1)+1:cutind(2)-2; 
if isempty(useind)
    warning('Unable to calculate InvIPRcut, No valid MT indices after removing 1 MT from the ends.')
    IPRcut = NaN;
else
    logdistintcut = logdistint(useind);
    maxCut = max(logdistintcut);
    logNormConstCut = maxCut + log(sum(exp(logdistintcut - maxCut)));
    normConstcut = exp(logNormConstCut);

    logdistint2cut = log(abs(Deff./(2*vels))) + log(abs(exp(2*Pe) - 1)) + 2*cumPe_reg + 2*log(z);
    logdistint2cut(ind0) = log(reglen(ind0)) + 2*cumPe(ind0) + 2*log(z(ind0));
    logdistint2cut = logdistint2cut(useind);
    logdistint2Max_cut = max(logdistint2cut);
    logsumdistint2cut = log(sum(exp(logdistint2cut - logdistint2Max_cut))) + logdistint2Max_cut;
    logIPRcut = logsumdistint2cut - 2*logNormConstCut;
    IPRcut = exp(logIPRcut);
end


end