function group_id = make_groups(z, nbins, method)

z = z(:);

switch lower(method)
    case 'quantile'
        q = quantile(z, linspace(0,1,nbins+1));
        q(1) = q(1) - eps;
        q(end) = q(end) + eps;

        % 防止重复边界导致discretize报错
        q = unique(q, 'stable');
        if numel(q) < 2
            group_id = ones(size(z));
        else
            group_id = discretize(z, q);
        end

    case 'equal'
        edges = linspace(min(z), max(z), nbins+1);
        edges(1) = edges(1) - eps;
        edges(end) = edges(end) + eps;
        group_id = discretize(z, edges);

    case 'kmeans'
        group_id = kmeans(z, nbins, 'Replicates', 10);

    otherwise
        error('Unknown grouping method.');
end

end