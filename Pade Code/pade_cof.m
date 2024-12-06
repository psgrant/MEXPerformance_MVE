function [N, D] = pade_cof(d, k)
% Re-normalised (d,d)-PadÈ coefficients for the \phi_\ell functions.
%
%%%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% This code is based on the explicit (d,d)-PadÈ formulae derived by W.
% Wright, but uses recurrence relations to reduce the computational
% complexity.  In particular, the coefficients n_i of the numerator
% polynomial N_d^\ell are defined by
%
%    n_i = \sum_{j=0}^i a_{ij} = sum(tril(A), 2)_i
%
% for i=0:d, in which
%
%    a_{ij} = (2d + ell - j)! (-1)^j / (j! (d-j)! (ell + i - j)!)
%           = -(d+1 - j) (ell+1+i - j) / (j (2d+ell+1 - j)) a_{i,j-1}
%
% for j=1:i.  Similar recurrence relations may be derived for the other
% coefficients, and for the denominator polynomial D_d^\ell.
%
% We note that roundoff errors unfortunately affects the accuracy of the
% coefficients.  However, as the errors are generally in the order of
% 1-5 ULP, we do not implement more accurate evaluation routines at this
% time.

n1 = prod(d + 1 : 2*d + 1);     % (2d + 1)! / d!
d1 = n1;

i = 1:d;
[J, I] = meshgrid(i);   % MESHGRID gives wrong order for this purpose

N = cell([k, 1]);
D = cell([k, 1]);
A = zeros(d + 1);

ell = 1;
while ell <= k,
    A(:, 1) = n1 .* cumprod([1, 1 ./ (ell + i)]) .';
    A(2:end, 2:end) = - (d + 1 - J) .* (ell + 1 + I - J) ./ ...
        ((2*d + ell + 1 - J) .* J);
    
    N{ell} = sum(tril(cumprod(A, 2)), 2);
    D{ell} = d1 .* cumprod([1, -(d + 1 - i) ./ (i .* (2*d + ell + 1 - i))]) .';
    
    ell = ell + 1;
    
    n1 = n1 * (2*d + ell) / ell;
    d1 = d1 * (2*d + ell);
end
end


