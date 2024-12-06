function [idx, phi] = put_cache(phi, z, k, idx, d)
% - request not satisfied from cache, do complete eval, set phi.scaled
% - further optimisations are possible in this case but not without
%   adversely affecting code readability...

siz  = size(z);
maxk = max(k);

phi_list = eval_pade(z, k, d);

s = struct('siz', siz, 'nrm_z', norm(z, inf), ...
    'maxk', maxk, 'z', z, 'phi', cat(2, phi_list{1}{:}));

phi.scaled = struct('siz', 0, 'nrm_z', 0, 'maxk', 0, 'z', [], 'phi', []);
for j = 2:numel(phi_list),
    z = z / 2;
    phi.scaled(j-1) = struct('siz',   siz,          ...
        'nrm_z', norm(z, inf), ...
        'maxk',  maxk,         ...
        'z', z, 'phi', cat(2, phi_list{j}{:}));
end

if isempty(phi) || ~isfield(phi, 'LRU'),
    phi.LRU.list = [];
    phi.LRU.idx  = [];
end

idx = abs(idx);
n = numel(phi.LRU.list);
if (0 < idx) && (idx <= n),
    % - new phi-list computed for additional \phi_\ell functions
    % - replace existing phi.LRU z-entry
    
    phi.LRU.list(phi.LRU.idx(idx)) = s;
    phi.LRU.idx = phi.LRU.idx([idx, 1:idx-1, idx+1:n]);
else
    % - either not found or insufficient \phi_\ell functions
    % - prepend the newly computed phi functions to existing list
    
    phi.LRU.list = [s, phi.LRU.list];
    phi.LRU.idx  = [1, 1 + phi.LRU.idx];
end

idx = 1;        % entry found at phi.LRU.idx(1)

end