function [wall,allbounds,MTpos,plusdir,mindir,nMTdomain,MTdomain,MTstartbound,MTendbound,MTsortind] = ...
    getConfigWallsModified_fullnosort(x_coords,marker_positions,snaptol,domainends)
% Function getConfigWallsModified serves to sort out the x_coord and
%   markerposition inputs
%%Inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x_coords: Cell array of 2x1 matrices, for MT tips location. 
% for each MT, the x_coords are assumed listed from low value (start) to
% high value (end)
% marker_positions: Cell array of strings, 'right' and 'left' indicates MT plus ends
%   on the right and left hand side respectively.
% snaptol : scalar, x-position MT ends snapping tolerance. Default = 0.1
% domainends: 2x1 optional input, if provided, create empty regions on either
% side so that the domain stretches between these endpoints
%%Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wall : n x2 logical array, 'wall' to the left, means no non-diffusive particles can move left from
%   that bound
% allbounds : n x 1 scalar, locations of each MT ends
% MTpos : n x 2 scalar, x-locations of each MT start and end pos at each row
% plusdir: n x 1 logic, 1 if second column of MTpos is a plus end
% mindir: n x 1 logic, 1 if first column of MTpos is a plus end
% nMTdomain: n x 1 scalar, number of MT per domain
% MTdomain: 1 x n cell array, each cell contains the MT corresponding to
%   that region/domain. 
% MTstartbound: for each MT, which is its leftmost boundary
% MTendbound: for each MT, which is its rightmost boundary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do we bother to calculate MTdomain?
getMTdomain = nargout>6;

% nMT = length(x_coords);
% 
% % Logical indexing to filter out 'none' marker positions
% valid_indices = ~strcmp(marker_positions, 'none');
% % Filter x_coords based on valid indices
% filtered_x_coords = x_coords(valid_indices);
% % Compute nMT based on valid entries
% nMT_filtered = length(filtered_x_coords);

% covert to more convenient datatypes
% nx2 matrix of start and and points for each MT
% sort by starting (leftmost)  points
% [MTpos,ind] = sortrows(vertcat(x_coords{:}));
% MTsortind = ind;
% [MTpos_filtered, ind_filtered] = sortrows(vertcat(filtered_x_coords{:}));
% [~, indices] = ismember(MTpos_filtered, MTpos, 'rows');
% 
% % list of logical values. 1 when MT points in the positive direction
% plusdir = cellfun(@(x) strcmp(x,'right'),marker_positions(1:nMT));
% plusdir = plusdir(ind)';
% 
% mindir = cellfun(@(x) strcmp(x,'left'),marker_positions(1:nMT));
% mindir = mindir(ind)';

nMT = length(x_coords);

% Valid MTs are those not marked as 'none'
valid = ~strcmp(marker_positions(:), 'none');

% Keep MTpos in original x_coords order
MTpos = nan(nMT,2);
MTpos(valid,:) = vertcat(x_coords{valid});

% MT row index is now original x_coords index
MTsortind = (1:nMT).';

% Direction arrays also stay in original order
plusdir = strcmp(marker_positions(:), 'right');
mindir  = strcmp(marker_positions(:), 'left');

% Inactive MTs should not have direction
plusdir(~valid) = false;
mindir(~valid) = false;
%%
% plusdir_filtered = plusdir(indices);
% mindir_filtered = mindir(indices);
activeMTpos = MTpos(valid,:);

if isempty(activeMTpos)
    error('No active MTs left.')
end
%% ----------------
% define boundary positions (intermediate between highest and lowest for
% nearby values
if abs(snaptol) > 0
    allbounds_low = uniquetol(activeMTpos(:),snaptol/max(abs(activeMTpos(:))),'lowest');
    allbounds_high = uniquetol(activeMTpos(:),snaptol/max(abs(activeMTpos(:))),'highest');
    allbounds = (allbounds_low+allbounds_high)/2;
else
    allbounds = unique(activeMTpos(:));
end
% 
% snapped1 = interp1(allbounds,allbounds, MTpos(:,1), 'nearest', 'extrap');
% snapped2 = interp1(allbounds,allbounds, MTpos(:,2), 'nearest', 'extrap');
% 
% % coordinates should always go left to right
% MTpos(:,1) = min(snapped1,snapped2);
% MTpos(:,2) = max(snapped1,snapped2);
%% Snapping logic above does not account for snaptol, instead snaps to nearest boundary.
%% adjust endpoints to the snapped bounds (old slow version)
% if snaptol == 0
%     flipIdx = MTpos(:,1) > MTpos(:,2);
%     MTpos(flipIdx,:) = fliplr(MTpos(flipIdx,:));
%     % for mc = 1:nMT
%     %     if (MTpos(mc,1) > MTpos(mc,2))
%     %         MTpos(mc,:) = fliplr(MTpos(mc,:));
%     %     end
%     % end
% else
%     for mc = 1:nMT
%         dists1 = abs(MTpos(mc,1)-allbounds);
%         % find a boundary within tolerance and snap to it
%         ind = find(dists1<snaptol,1);
%         if ~isempty(ind1)
%             MTpos(mc,1) = allbounds(ind1);
%         end
% 
%         dists2 = abs(MTpos(mc,2)-allbounds);
%         % find a boundary within tolerance and snap to it
%         ind = find(dists2<snaptol,1);
%         if ~isempty(ind2)
%             MTpos(mc,2) = allbounds(ind2);
%         end
% 
%         % Reorder coordinates to always go from left to right
%         if (MTpos(mc,1) > MTpos(mc,2))
%             MTpos(mc,:) = fliplr(MTpos(mc,:));
%         end
%     end
% end
%% Snap MT endpoints to nearby boundaries
if snaptol == 0
    flipIdx = valid & MTpos(:,1) > MTpos(:,2);
    MTpos(flipIdx,:) = fliplr(MTpos(flipIdx,:));
else
    for mc = find(valid).'
        dists1 = abs(MTpos(mc,1) - allbounds);
        ind1 = find(dists1 < snaptol,1);
        if ~isempty(ind1)
            MTpos(mc,1) = allbounds(ind1);
        end

        dists2 = abs(MTpos(mc,2) - allbounds);
        ind2 = find(dists2 < snaptol,1);
        if ~isempty(ind2)
            MTpos(mc,2) = allbounds(ind2);
        end

        if MTpos(mc,1) > MTpos(mc,2)
            MTpos(mc,:) = fliplr(MTpos(mc,:));
        end
    end
end
%x_coords_fix = num2cell(MTpos,2)';
%plotMTConfig(x_coords_fix,marker_positions(1:4),params)
%%
% valid__MT_indices = zeros(size(allbounds) - 1); % One less than allbounds since it's for intervals
% interval_index = {};
% % Process each interval
% for i = 1:length(allbounds) - 1
%     % Get the start and end of the current interval
%     % start_bound = x_coords{i}(1);    
% 
%     % Find the interval index in allbounds
%     interval_index{i} = find(MTpos(:,2) >= allbounds(i+1) & MTpos(:,1) <= allbounds(i));
% 
%     plusdir_at_interval = plusdir(interval_index{i});
%     mindir_at_interval = mindir(interval_index{i});
%     % Update valid_indices based on marker_positions
%     if (any(plusdir_at_interval) == 1) || (any(mindir_at_interval) == 1)
%         valid__MT_indices(i) = 1;
%     else
%         valid__MT_indices(i) = 0;
%     end
% end

%% Go through each bound and decide if it a 'wall' in each direction
% 'wall' to the left, means no non-diffusive particles can move left from
% that bound
nbound = length(allbounds);
wall = true(nbound,2);
MTdomain = cell(1,nbound-1,1);
nMTdomain = zeros(nbound-1,1);
MTendbound = zeros(size(MTpos,1),1);
MTstartbound = zeros(size(MTpos,1),1);

for bc = 1:nbound

    % microtubules with left end on this boundary
    Lend = valid& abs(MTpos(:,1)-allbounds(bc))<=snaptol + eps ;
    % microtubules with right end on this boundary
    Rend = valid& abs(MTpos(:,2)-allbounds(bc))<=snaptol + eps;
    % microtubules that cross over this boundary
    Cross = valid& MTpos(:,1)<allbounds(bc)-(snaptol + eps) & MTpos(:,2)>allbounds(bc)+ (snaptol + eps);

    % not a wall left if any minus-directed MTs have right end here or
    % cross over
    wall(bc,1) = ~any(mindir & (Cross | Rend));
    
    % not a wall right if any plus-directed MTs have left end here or
    % cross over
%     wall(bc,2) = ~any(plusdir_filtered' & (Cross | Lend)); 
    wall(bc,2) = ~any(plusdir & (Cross | Lend));

    if (bc < nbound)
        % number of microtubules in the domain beyond this bound
       % if valid__MT_indices(bc) > 0
            nMTdomain(bc) = nnz((Cross |Lend));
            
            if (getMTdomain)
                % list of MT indices in this domain
                MTdomain{bc} = find(Cross | Lend);
            end
       % end
    end

    % for each MT what boundary does it start and end on (does not consider
    % directionality, just what is the left-most and right-most endpoint)
    MTstartbound(Lend) = bc;
    MTendbound(Rend) = bc;    
end

if (exist('domainends','var'))
    % create empty regions on either side
    wall = [1 1; wall; 1 1];
    allbounds = [domainends(1); allbounds; domainends(2)];
    nMTdomain = [0; nMTdomain; 0];
    if (getMTdomain)
        MTdomain = {zeros(0,1), MTdomain{:}, zeros(0,1)};
    end

    MTstartbound = MTstartbound+1;
    MTendbound = MTendbound+1;
end