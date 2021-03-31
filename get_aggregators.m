% Function input ind_sets is a cell array, each element a vector showing 
% the set of indices that need to be aggregated over

% These need to be identified in separate groups to allow a different row
% of agg to correspond to each group

function [sel, agg] = get_aggregators(ind_sets, nstates)

nrows = length(ind_sets);

% Selector
sel = zeros(nstates);
for ii = 1:nrows
    sel(ind_sets{ii},:) = 1;
end
sel(logical(eye(size(sel)))) = 0;
sel = sparse(sel);

% Aggregator
agg = zeros(nrows, nstates);
for ii = 1:nrows
   agg(ii,ind_sets{ii}) = 1;
end
agg = sparse(agg);