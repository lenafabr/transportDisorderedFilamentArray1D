function plotMTposModified(MTpos, plusdir, mindir, allbounds)
N = size(MTpos, 1);

ycoords = 1:N;

% Conditionally plot for plusdir or mindir
for mc = 1:N
    if plusdir(mc)
        plot(MTpos(mc, :), ycoords(mc) * ones(1, 2),'Color',[0.5 0 0]);
        hold on;
    elseif mindir(mc)
        plot(MTpos(mc, :), ycoords(mc) * ones(1, 2),'Color',[0 0 0.5]);
        hold on;
    end
end

% plusdir
% scatter(MTpos(plusdir, 1), ycoords(plusdir), 30, 'k');
scatter(MTpos(plusdir, 2), ycoords(plusdir), 30, 'k', 'filled');

% mindir
scatter(MTpos(mindir, 1), ycoords(mindir), 30, 'k', 'filled');
% scatter(MTpos(mindir, 2), ycoords(mindir), 30, 'k');

% If not mindir or plusdir its an artificial boudndary
neitherDir = ~(plusdir | mindir); 
% scatter(MTpos(neitherDir, 1), ycoords(neitherDir), 30, 'k','x');
% scatter(MTpos(neitherDir, 2), ycoords(neitherDir), 30, 'k', 'x'); 

% Plot allbounds as vertical dashed lines
for bc = 1:length(allbounds)
    % plot(allbounds(bc) * ones(1, 2), [0, N], 'k:','LineWidth',0.3);
end

hold off;
end
