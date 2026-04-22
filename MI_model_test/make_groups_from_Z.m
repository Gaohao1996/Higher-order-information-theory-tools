function [group_id, edges] = make_groups_from_Z(Z, n_groups, method)

if nargin < 3
    method = 'quantile';
end

Z = Z(:);

switch lower(method)
    case 'quantile'
        edges = quantile(Z, linspace(0,1,n_groups+1));
        edges(1) = edges(1) - 1e-12;
        edges(end) = edges(end) + 1e-12;
        edges = unique(edges, 'stable');

    case 'equal'
        edges = linspace(min(Z), max(Z), n_groups+1);
        edges(1) = edges(1) - 1e-12;
        edges(end) = edges(end) + 1e-12;

    otherwise
        error('Unknown grouping method.');
end

if numel(edges) < 2
    group_id = ones(size(Z));
    return;
end

group_id = zeros(size(Z));

for g = 1:(numel(edges)-1)
    idx = Z > edges(g) & Z <= edges(g+1);
    group_id(idx) = g;
end

end