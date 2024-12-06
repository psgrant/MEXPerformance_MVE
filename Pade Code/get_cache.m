function [phi, varargout] = get_cache(phi, idx, k, bound)
% - idx > 0
% - numel(phi.LRU.list) > 0
% - idx < numel(phi.LRU.list) || ~isempty(phi.scaled(end).z)

[nl, ns] = deal(numel(phi.LRU.list), numel(phi.scaled));
varargout = cell([1, numel(k)]);
if idx <= nl + ns;
    % request directly satisfied from phi.LRU or phi.scaled
    
    if idx <= nl,
        p = phi.LRU.list(phi.LRU.idx(idx));
        phi.LRU.idx = phi.LRU.idx([idx, 1:idx-1, idx+1:nl]);
    else
        idx = idx - nl;
        p = phi.scaled(idx);
        phi.LRU.list = [phi.scaled(idx), phi.LRU.list];
        phi.LRU.idx  = [1, 1 + phi.LRU.idx];
    end
    
    siz  = p.siz;
    maxk = p.maxk;
    
    phi_vals = mat2cell(p.phi, siz(1), repmat(siz(2), [1, maxk]));
else
    % z = 2^p * phi.scaled(end).z, but outside of phi.scaled
    % must do a bit of squaring
    
    s  = idx - (nl + ns) + 1;    % number of squarings needed
    p  = phi.scaled(1);
    
    siz  = p.siz;
    maxk = p.maxk;
    Id   = speye(siz);
    pv   = mat2cell(p.phi, siz(1), repmat(siz(2), [1, maxk]));
    
    [phi_vals, d{1:s}] = square_pade(p.z, Id, s, pv);
    z = 2^s .* p.z;
    
    phi.LRU.list = [struct('siz', siz,            ...
        'nrm_z', norm(z, inf), ...
        'maxk', max(k),        ...
        'z', z,                ...
        'phi', cat(2, phi_vals{:})), ...
        phi.LRU.list];
    phi.LRU.idx  = [1, 1 + phi.LRU.idx];
end

% maintain upper bound on number of cached entries
if numel(phi.LRU.list) > bound,
    phi.LRU.list = phi.LRU.list(phi.LRU.idx(1:bound));
    phi.LRU.idx  = 1:bound;      % data copying normalises LRU index
end

% return values
[varargout{1:numel(k)}] = deal(phi_vals{k});

end