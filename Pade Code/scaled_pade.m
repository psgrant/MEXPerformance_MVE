function varargout = scaled_pade(z, Id, d, k)
evn = 2 * (0 : floor(d / 2));
d_evn = numel(evn) - 1;
[s_evn, r_evn] = partition_polynomial(d_evn);

odd = 1 + 2*(0 : floor((d - 1) / 2));
d_odd = numel(odd) - 1;
[s_odd, r_odd] = partition_polynomial(d_odd);

Z = cell([s_evn + 1, 1]);       % s_evn >= s_odd
P = cell([2, 1]);

% assume d >= 2
Z{1} = Id;
Z{2} = z * z;

% we only need even powers of z
for j = 3:s_evn + 1, Z{j} = Z{j-1} * Z{2}; end

[num, den] = pade_cof(d, k);
varargout  = cell([k, 1]);

% PadÈ eval algorithm due to Higham in
%   The Scaling and Squaring Method for the Matrix Exponential Revisited
for ell = 1:k,
    N =       mat_pol(num{ell}(evn+1), Z, d_evn, s_evn, r_evn);
    N = N + z*mat_pol(num{ell}(odd+1), Z, d_odd, s_odd, r_odd);
    
    D =       mat_pol(den{ell}(evn+1), Z, d_evn, s_evn, r_evn);
    D = D + z*mat_pol(den{ell}(odd+1), Z, d_odd, s_odd, r_odd);
    
    varargout{ell} = D \ N;
end
end