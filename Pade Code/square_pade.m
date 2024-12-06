function varargout = square_pade(z, Id, s, phi)
% Undo scaling (squaring).
%
% Formulae derived by W. Wright.

[varargout{1:s+1}] = deal(phi);         % prealloc
Exp = z*phi{1} + Id;                    % exponential

for m = 1:s,
    varargout{end-m}{1} = (Exp + Id) * phi{1} / 2;
    
    i = [0, 1];
    p = ~i;
    for k = 2:numel(phi),
        i = i + p;
        v = phi{i(1)} * phi{i(2)};
        
        c = 2;
        a = mod(k, 2);
        ell = floor(k / 2);
        for j = k : -1 : ell + 1 + a,
            v = v + c.*phi{j};
            c = c / (k + 1 - j);
        end
        
        % odd-numbered phi's need special coeff in \phi_{\ell+1} term
        if a > 0, v = v + phi{ell+1}./prod(1:ell); end
        
        varargout{end-m}{k} = v / 2^k;
        p = ~p;
    end
    
    [phi{:}] = deal(varargout{end-m}{:});
    Exp = Exp * Exp;
end
end
