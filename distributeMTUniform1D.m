function [x_coords,marker_positions] = distributeMTUniform1D(nMT,Ltot,Ldom,f,constrain)
% generate a random distribution of MTs
% uniformly distributed
% with a specified fraction pointing in each direction
% the MTs are assumed to all be equal in length
% input:
% nMT = number of MTs
% Ltot = total length of all MTs
% Ldom = domain length
% f = fraction pointed in the positive direction
% output:
% x_coords = cell array of 1x2 coordinates for each MT
% marker_positions = cell array of 'left' or 'right' for each MT

% constrain = if true, the force exactly fraction f to point right
% otherwise, each MT independently chooses left or right with appropriate
% probability
if (~exist('constrain'))
    constrain = true;
end

% length of each MT;
MTlen = Ltot/nMT;

% distribute left end-points of the MTs
leftpts = rand(nMT,1)*(Ldom-MTlen);
% place right end-points
rightpts = leftpts + MTlen;
allcoords = [leftpts rightpts];

% decide which way the MTs point
if (constrain)
    % how many point right.
    nright = round(f*nMT);
    % which ones point right
    rightind = datasample(1:nMT,nright,'Replace',false);

    MTdirs = zeros(1,nMT);
    MTdirs(rightind) = 1;
else
    % sample each MT independently
    uvals = rand(1,nMT);
    isright = uvals<f;
    MTdirs = zeros(1,nMT);
    MTdirs(isright) = 1;
end

% convert to cell arrays
for mc = 1:nMT
    x_coords{mc} = allcoords(mc,:);
    if (MTdirs(mc)==0)
        marker_positions{mc} = 'left';
    else
        marker_positions{mc} = 'right';
    end
end