function idx = find_cache(phi, z, k, atol, rtol)
idx = 0;

if ~isempty(phi)
    nrm_z = norm(z, inf);
    siz   = size(z);
    siz_c = numel(phi.LRU.list);
    
    % likely non-optimal, but use linear search for now...
    j = 0;
    while (idx == 0) && (j < siz_c),
        j = j + 1;
        p = phi.LRU.list(phi.LRU.idx(j));
        
        if match_cache(p, z, nrm_z, siz, atol, rtol),
            if p.maxk >= k,
                idx = j;            % all requirements satisfied
            else
                idx = -j;           % correct z, too few phi's
            end
        end
    end
    
    if idx == 0,
        j = match_scaled(phi.scaled, z, k, atol, rtol);
        
        % relies on sign(j)==0 for j==0
        idx = sign(j) * (siz_c + abs(j));
    end
end

end