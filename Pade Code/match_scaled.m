function idx = match_scaled(phi, z, k, atol, rtol)
idx = 0;

if ~isempty(phi),
    [s, z] = scaled_arg(z);
    sgn = 1;
    
    if match_cache(phi(end), z, norm(z, inf), numel(z), atol, rtol),
        if phi(end).maxk < k, sgn = -1; end
        
        if s <= numel(phi),
            idx = numel(phi) - s;
        else
            idx = s;
        end
        
        idx = sgn * idx;
    end
end

end