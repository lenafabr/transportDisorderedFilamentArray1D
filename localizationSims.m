function pgroup = localizationSims(param_in,MTparam_in,pgroup)
% Function localizationSims takes in kinetic and geometric
% MT parameters to simulate cargo moving around a MT array
%%Inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% param_in : structure containing kinetic and simulation parameters
%   .D          : scalar, diffusion coefficient
%   .v          : scalar, velocity
%   .ka         : scalar, attachment rate
%   .kd         : scalar, detachment rate rate
%   .npart      : scalar, number of particles
%   .maxevent   : scalar, max number of events allowed for a particle
%   .dtscale    : scalar, factor to scale dt. 
%   .savetimes  : row vec, savetimes for particle positions 
%   .verbose    : scalar, 0 = print minimal info, else print detailed info
%   .printevery : scalar, event interval for printing simulation details. 
% MTparam_in : structure containing MT-related parameters
%   .snaptol : scalar, x-position MT ends snapping tolerance. Default = 0 (
%   no snapping)
%   .domain: 2x1 matrix of scalars. Entire system size. 
%   .x_coords: Cell array of 2x1 matrices, for MT tips location. 
%   .marker_positions: Cell array of strings, 'right' and 'left' indicates MT plus ends
%   on the right and left hand side respectively.
% pgroup = optional, start with an existing particle group (Initial condition)
% if not provided, start uniformly distributed in diffusive state
%%Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prgroup : Simulation data containing particle position and times
%   .npart     : scalar, number of particles
%   .curtime   : col vec (npart,1), current time of particle
%   .pos       : col vec (npart,1), current particle position
%   .state     : col vec (npart,1), current state of particle (states defined
%   below)
%   .MTind     : col vec (npart,1), current MT particle is attached to.
%   .savepos   : col vec (npart, p), position of all particle at all .savetimes.  
%   .savestate : col vec (npart,p,2), state and MT particle is attached to
%   for all particles, at all .savetimes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% snap domain endpoints, convert to more convenient datatypes
% MTdomain = list of MT indices in each domain
% plusdir = list of logicals for whether each MT is pointed in the +
% direction
% MTpos and plusdir are both sorted based on the left-hand coord of each MT
nMT = length(MTparam_in.x_coords);
[~,allbounds,MTpos,plusdir,~,nMTdomain,MTdomain,MTstartbound,MTendbound] = ...
    getConfigWallsModified_full(MTparam_in.x_coords,MTparam_in.marker_positions,MTparam_in.snaptol,MTparam_in.domain);
%% set up convenient arrays to store information about the regions

% MTbound = MT that particles can attach to, at that given region
MTbound = cell(1, length(allbounds));
MTbound{1} = MTdomain{1}; 
nMTbound(1) = length(MTbound{1});

for rc = 1:length(MTdomain)
    MTs = MTdomain{rc};

    % list of + and - MTs in this region
    MTplusreg{rc} = MTs(plusdir(MTs));
    MTminusreg{rc} = MTs(~plusdir(MTs));
    nMTplus(rc) = length(MTplusreg{rc});
    nMTminus(rc) = length(MTminusreg{rc});

    % MTs touching or crossing boundary above this region
    if (rc == length(MTdomain))
        MTbound{rc+1} = MTdomain{rc};
    else
        MTbound{rc+1} = union(MTdomain{rc},MTdomain{rc+1});
    end
    nMTbound(rc+1) = length(MTbound{rc+1});

    % total on rate in each region
    onrates(rc) = param_in.ka*nMTdomain(rc);

    % length of each region
    reglen(rc) = allbounds(rc+1)-allbounds(rc);
end

nreg = length(allbounds)-1; % total number of regions
reglen = reglen(:); % region lengths for all region
% -----------
%% define and initialize a structure describing a group of particles
if (~exist('pgroup','var'))
    pgroup = struct();
    % number of particles
    pgroup.npart = param_in.npart;

    % initialize the particles

    % particle positions
    % initialize uniformly distributed
    pgroup.pos = rand(pgroup.npart,1)*(MTparam_in.domain(2) - MTparam_in.domain(1)) + MTparam_in.domain(1);
    
    % initialize particle states:
    % (0) diffusive
    % (1) moving to the right
    % (-1) moving to the left
    % (2) stuck at the positive tip of a MT
    pgroup.state = zeros(pgroup.npart,1);
    % which MT are they on? 0 means not on any
    pgroup.MTind = zeros(pgroup.npart,1);

    % saved positions at specified timepoints
    pgroup.savepos = NaN*zeros(pgroup.npart,length(param_in.savetimes));
    % save the state and which MT is is on
    pgroup.savestate = NaN*zeros(pgroup.npart,length(param_in.savetimes),2);

    % current time for each particle
    pgroup.curtime = zeros(pgroup.npart,1);
end

%% cycle through events for all particles

% for each particle, what is the index of the next savetime
nextsave = ones(pgroup.npart,1);
% add an extra element to savetimes to simplify bookkeeping
savetimes = param_in.savetimes;
%
% boundary index will only be defined for particles sitting on a MT tip
curbound = nan(param_in.npart,1);

% finished simulation to the last savepoint
finished = false(param_in.npart,1);
tic
for eventc = 1:param_in.maxevent
    % stop if all particles are done
    if (all(finished)); break; end

    if (param_in.verbose == 0 & mod(eventc,param_in.printevery)==0)
        disp(sprintf('event #, # particles left, total time: %d, %d, %g', eventc, nnz(~finished), min(pgroup.curtime)))
        disp(sprintf('cumulative time spent = %g', toc))
    end

    % check that no particles are out of range of the MTs
    % (this is mostly for debugging purposes)
    notdiff = pgroup.state>0;
    outofrange = (pgroup.pos(notdiff) > MTpos(pgroup.MTind(notdiff),2)+eps | ...
        pgroup.pos(notdiff) < MTpos(pgroup.MTind(notdiff),1)-eps);
    if (any(outofrange))
       tmp = 1:pgroup.npart;
       tmp = tmp(notdiff);
       badpc = tmp(outofrange);
       error(sprintf('some particle went out of range! %d',badpc(1)))
    end
    if (any(isnan(pgroup.pos)))
        error('got NaN in pos')
    end

    % keep track of the type of event
    % 0 = unbind
    % 1 = bind
    % 2 = hit MT end
    % 3 = hit domain boundary
    % 5 = hit savetime
    eventtype = nan(1,pgroup.npart);

    dtevent = zeros(pgroup.npart,1);
    % time until next savepoint for each particle
    dtsave(~finished) = savetimes(nextsave(~finished)) - pgroup.curtime(~finished)';
    dtsave = dtsave(:);
    % what region are the particles currently in?
    % deal with precision issues
    curreg = interp1(allbounds,1:(nreg+1),pgroup.pos,'previous');

    % particles are exactly on right boundary of the region and should be
    % shifted past it
    onboundR = (abs(pgroup.pos - allbounds(curreg+1))<10*eps) & (pgroup.state==1 | pgroup.state==2) & curreg < nreg;
    curreg(onboundR) = curreg(onboundR)+1;
    % particles are exactly on left boundary of the region and should be
    % shifted past it
    onboundL = (abs(pgroup.pos - allbounds(curreg))<10*eps) & (pgroup.state==-1 | pgroup.state==2) & curreg >1;
    curreg(onboundL) = curreg(onboundL)-1;

    if (param_in.verbose>1 | param_in.verbose < 0)
        if (param_in.verbose>0)
            disp([eventc pgroup.pos' pgroup.state' pgroup.MTind' curreg'])
        else
            pc = -param_in.verbose;
            disp([eventc pgroup.pos(pc)' pgroup.state(pc)' pgroup.MTind(pc)' curreg(pc)'])
        end
    end

    % keep track of which particles have finished the step
    donestep = false(pgroup.npart,1);
    donestep(finished) = true; % some particles may already be finished

    % need to update savepos for these particles?
    dosave = false(pgroup.npart,1);

    % ------------  ----------------------
    %% diffusive particles
    % ------------------------------------
    % get particles currently in diffusive state
    p0 = find(pgroup.state == 0 & ~donestep);

    % get appropriately small diffusive timestep for each particle
    L = reglen(curreg(p0)); % lengths of current domains
    % diffusive time is smaller than the time to cover the domain length
    dtdiff = L.^2/(2*param_in.D)*param_in.dtscale;
    % dtdiff = dtdiff';
    % on-rates for these particles
    kon = onrates(curreg(p0));
    kon = kon(:);

    % decide whether particles come on before timestep is over
    ontimes = exprnd(1./kon);
    % mark particles that came on before anything else happens
    onfirst = (ontimes<dtdiff & ontimes<dtsave(p0));
    savefirst = (dtsave(p0)<ontimes & dtsave(p0)< dtdiff);
    dtevent(p0(onfirst)) = ontimes(onfirst);
    dtevent(p0(savefirst)) = dtsave(p0(savefirst));
    dtevent(p0(~onfirst & ~savefirst)) = dtdiff(~onfirst & ~savefirst);

    % propagate in space
    dx = randn(length(p0),1).*sqrt(2*param_in.D*dtevent(p0));
    % deal with reflecting boundaries of full domain
    newpos = pgroup.pos(p0) + dx;
    newpos(newpos < MTparam_in.domain(1)) = 2*MTparam_in.domain(1) - newpos(newpos < MTparam_in.domain(1));
    newpos(newpos > MTparam_in.domain(2)) =   2*MTparam_in.domain(2) - newpos(newpos > MTparam_in.domain(2));
    pgroup.pos(p0) = newpos;


    % update MT
    p0on = p0(onfirst);
    for pc = p0on'
        % which MT did it jump onto?
        MTc = randi(nMTdomain(curreg(pc)));
        pgroup.MTind(pc) = MTdomain{curreg(pc)}(MTc);
    end
    % update state
    pgroup.state(p0on) = plusdir(pgroup.MTind(p0on))*2-1;

    % if escaped region and came on during this timestep, keep within the
    % MT range 
    % passed end of MT to the right or left
    MTind = pgroup.MTind(p0on);
    passR = (pgroup.pos(p0on) > MTpos(MTind,2));
    passL = (pgroup.pos(p0on) < MTpos(MTind,1));

    pgroup.pos(p0on(passR)) = MTpos(MTind(passR),2);
    pgroup.pos(p0on(passL)) = MTpos(MTind(passL),1);

    hitendR = (passR & pgroup.state(p0on)==1);
    hitendL = (passL & pgroup.state(p0on)==-1);
    hitend = hitendR | hitendL;
    pgroup.state(p0on(hitend)) = 2;
    eventtype(p0on(hitend)) = 2; % hit MT boundary

    % update what boundary you are on
    curbound(p0on(hitendR)) = MTendbound(MTind(hitendR));
    curbound(p0on(hitendL)) = MTstartbound(MTind(hitendL));

%    pgroup.pos(p0on) = min(MTpos(pgroup.MTind(p0on),2),pgroup.pos(p0on));
%    pgroup.pos(p0on) = max(MTpos(pgroup.MTind(p0on),1),pgroup.pos(p0on));

    % need to update savepos for these particles
    dosave(p0(savefirst)) = true;

    % keep track that these have already finished the step
    donestep(p0) = true;
    eventtype(p0on) = 1; % bind
    eventtype(p0(~savefirst & ~onfirst)) = 3; % hit domain boundary
    eventtype(p0(savefirst)) = 5; % hit savetime

    % -----------------
    %% particles walking to the right or left
    % -----------------
    p1p = find(pgroup.state == 1 & ~donestep);
    p1m = find(pgroup.state == -1 & ~donestep);
    p1  = [p1p; p1m];

    if (~isempty(p1))
        dxbound_p = allbounds(curreg(p1p)+1) - pgroup.pos(p1p);
        dxbound_m = pgroup.pos(p1m) - allbounds(curreg(p1m));
        dxbound = [dxbound_p' dxbound_m'];
        dtwalk = dxbound/param_in.v;
        dtwalk = dtwalk(:);

        % sample times for the unbinding event
        kintimes = exprnd(1./param_in.kd, size(p1));

        % kinetic event (unbind) happened before savetime and before
        % reaching end of domain
        kinfirst = (kintimes < dtwalk) & (kintimes < dtsave(p1));
        dtevent(p1(kinfirst)) = kintimes(kinfirst);

        savefirst = (dtsave(p1) < kintimes & dtsave(p1) < dtwalk);
        dtevent(p1(savefirst)) = dtsave(p1(savefirst));
        eventtype(p1(savefirst)) = 5; % hit savetime

        p1d = p1(~kinfirst & ~savefirst);
        dtevent(p1d) = dtwalk(~kinfirst & ~savefirst);
        eventtype(p1d) = 3; % hit domain boundary

        % propagate in space
        pgroup.pos(p1) = pgroup.pos(p1) + param_in.v*pgroup.state(p1).*dtevent(p1);

        % unbind from the MTs
        pgroup.state(p1(kinfirst)) = 0;
        pgroup.MTind(p1(kinfirst)) = 0;
        eventtype(p1(kinfirst )) = 0; % unbind

        % reached end of the MT that it is on
        MTind = pgroup.MTind(p1d);
        % hit a left or right end
        hitendR = (plusdir(MTind) & MTendbound(MTind) == curreg(p1d) + 1);
        hitendL = (~plusdir(MTind) & MTstartbound(MTind) == curreg(p1d));
        hitend = hitendR | hitendL;
        pgroup.state(p1d(hitend)) = 2;
        eventtype(p1d(hitend)) = 2; % hit MT boundary

        % update what boundary you are on
        curbound(p1d(hitendR)) = MTendbound(MTind(hitendR));
        curbound(p1d(hitendL)) = MTstartbound(MTind(hitendL));

        dosave(p1(savefirst)) = true;
        donestep(p1) = true;
    end

    % -----------------
    %% particles stuck at MT end
    % -----------------
    p2 = find(pgroup.state == 2 & ~donestep)';
    p2 = p2(:);
    if (~isempty(p2))

        % sample times for the unbinding events
        kintimes = exprnd(1./param_in.kd, size(p2));
        
        % kinetic event (unbind) happened before savetime and before
        % reaching end of domain
        kinfirst = (kintimes < dtsave(p2));
        dtevent(p2(kinfirst)) = kintimes(kinfirst);

        savefirst = (dtsave(p2) < kintimes);
        dtevent(p2(savefirst)) = dtsave(p2(savefirst));
        eventtype(p2(savefirst)) = 5; % hit savetime

        % unbind from the MTs
        pgroup.state(p2(kinfirst)) = 0;
        pgroup.MTind(p2(kinfirst)) = 0;
        eventtype(p2(kinfirst)) = 0; % unbind

        dosave(p2(savefirst)) = true;
        donestep(p2) = true;
    end

    %% update save positions as needed
    psave = find(dosave);
    saveind = sub2ind(size(pgroup.savepos),psave,nextsave(psave));
    pgroup.savepos(saveind) = pgroup.pos(psave);
    % update the state and which MT the particle is on
    saveind1 = sub2ind(size(pgroup.savestate),psave,nextsave(psave),ones(size(psave)));
    pgroup.savestate(saveind1) = pgroup.state(psave);
    % save which MT they are on
    saveind2 = sub2ind(size(pgroup.savestate),psave,nextsave(psave),2*ones(size(psave)));
    pgroup.savestate(saveind2) = pgroup.MTind(psave);

    nextsave(psave) = nextsave(psave)+1;

    % update current time for each particle
    pgroup.curtime(~finished) = pgroup.curtime(~finished) + dtevent(~finished);

    % check which particles have finished
    finished(nextsave > length(savetimes)) = true;

    if (param_in.verbose>0 | param_in.verbose < 0)
        % display for all particles
        if (param_in.verbose > 0)
            disp(['N, events, times: ' sprintf('%d, ', eventc) sprintf('%d ',eventtype) sprintf('%g ',pgroup.curtime)])
        else
            % display for a particular particle only   
            pc = -param_in.verbose;
            disp(['N, events, times: ' sprintf('%d, ', eventc) sprintf('%d ',eventtype(pc))...
                sprintf('%g %g ',pgroup.curtime(pc),pgroup.pos(pc))])
        end
    end
end

end
