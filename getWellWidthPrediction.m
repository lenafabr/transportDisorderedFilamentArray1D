function wellwidth = getWellWidthPrediction(params,MTparam,cutoff,whichlim)
% estimate the well width from kinetic and configurational parameters
% cutoff defines how high to go to get out of the well
% whichlim = 1: assume width is less than a correlation length
% whichlim = 2: assume width is greater than a correlation length
% whichlim = 0 (default): pick the appropriate limit
% Includes factor of 2 for symmetric well(going out to either side)

% by default, pick the appropriate limit
if (~exist('whichlim'))
    whichlim = 0;
end

beta = params.v/params.kd;
gamma = sqrt(params.D/(MTparam.ncross*params.ka));
Lmt = MTparam.avglen;
ncross = MTparam.ncross;

% stdev of individual Pe
sigPe = Lmt*beta/(sqrt(2)*ncross^1.5*(beta^2 + gamma^2));
% correlated Pe step size
sPe = (2*ncross*sigPe)*sqrt(1/3 + 1/2/ncross + 1/4/ncross^2);

if (whichlim==1)
    wellwidth = (3)^(1/3)*(cutoff/(sigPe*2*ncross))^(2/3)*Lmt;
elseif (whichlim==2)
   
    % factor of 2 because we want the entire symmetric well
    % stretching out to both sides of the global minimum
    wellwidth = (1/3*(cutoff/sPe)^2  + 2/3*(cutoff/sPe))*Lmt;

elseif (whichlim==0)
    error('not done yet')
end

end