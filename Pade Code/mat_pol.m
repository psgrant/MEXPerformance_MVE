% Evaluate the (matrix) polynomial
%
%   p(z) = \sum_{j=0}^d b_j z^j
%
% in a somewhat optimised fashion.  (s,r) partitions (0:d) into r
% balanced sets of s terms with a possible d-rs extra terms at the end.
%
% The number of multiplications is minimised if s \approx sqrt(d), but
% the choice of s is up to the caller.  We explicitly assume d >= s*r,
% and note that s == 1 corresponds to the traditional Horner rule.
%
% Reference:
%  - ``Evaluating Matrix Polynomials'', section 11.2.4 (pp. 568-569) of
%    ``Matrix Computations'' by G.H. Golub and C.F. van Loan (3rd ed).

function p = mat_pol(b, Z, d, s, r)
j = d;
k = d - s*r;

p = b(j + 1) * Z{k + 1};

% modified Horner rule, backwards accumulation (high -> low degree)
while j > 0,
    while k > 0,
        j = j - 1;
        k = k - 1;
        p = p + b(j + 1)*Z{k + 1};
    end
    
    % prepare accumulation run only if there are any runs left
    if j > 0,
        k = s;
        p = p * Z{k + 1};
    end
end
end