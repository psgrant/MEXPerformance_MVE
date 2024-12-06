function B = flteq(x, y, varargin)
% FLTEQ - Determine floating point equality for arrays of DOUBLE or
%         COMPLEX.
%
% SYNOPSIS:
%   B = flteq(x, y);
%   B = flteq(x, y, atol);
%   B = flteq(x, y, atol, rtol);
%
% DESCRIPTION:
%   Using direct equality operators (== or ~=) may not be appropriate
%   for matrices of class DOUBLE or COMPLEX.  FLTEQ implements a weaker
%   sense of equality with exact equality definitions overridable by the
%   user.
%
%   Two objects, X and Y, are deemed equal if and only if
%    - ALL(SIZE(X) == SIZE(Y)), and
%    - (ALL(ABS(X - Y) < atol) or
%       ALL(ABS(X - Y) < rtol.*ABS(Y)))
%   with `atol' and `rtol' being absolute and relative tolerances
%   respectively.
%
%   For complex arrays, separate checks are made for the real and
%   imaginary parts.
%
% PARAMETERS:
%   x, y - Objects to check for equality.
%   atol - Absolute tolerance.  OPTIONAL.  DEFAULT VALUE = 1.0e-6.
%   rtol - Relative tolerance.  OPTIONAL.  DEFAULT VALUE = 1.0e-7.
%
% RETURNS:
%   B    - Boolean status indicating whether x is equal to y or not.
%          Possible values are  TRUE  and  FALSE.
%
% SEE ALSO:
%   RELOP, REAL, IMAG, TRUE, FALSE.

% DUPLICATION NOTE:
%
%   This function is an exact duplicate of FLOATEQUALS of the Expint
%   package.  The duplication is made in order to create a
%   self-contained PHIPADE function for use in other projects.  Any
%   change made to FLTEQ should be replicated in FLOATEQUALS if the
%   latter is available.

error(nargchk(2, 4, nargin));

if issparse(x) || issparse(y),
    % only work on non-zero elements...
    [ix, jx, x] = find(x);
    [iy, jy, y] = find(y);
    
    B = (numel(ix) == numel(iy)) && ...
        all(ix == iy) && all(jx == jy);
else
    B = true;                    % assume equality by default
end

sx = size(x);
sy = size(y);

if B && all(sx == sy),
    [atol, rtol] = deal(1.0e-6, 1.0e-7);
    
    if nargin > 2, atol = abs(varargin{1}); end
    if nargin > 3, rtol = abs(varargin{2}); end
    
    % Straighten out  x  and  y  for multi-D cases
    if all(sx > 1), x = x(:); y = y(:); end
    
    [xc, yc] = deal(real(x), real(y));
    xc = abs(xc - yc);
    a = all(xc < atol);
    r = all(xc < rtol.*abs(yc));
    
    if ~isreal(x) || ~isreal(y),
        [xc, yc] = deal(imag(x), imag(y));
        xc = abs(xc - yc);
        a = a & all(xc < atol);
        r = r & all(xc < rtol.*abs(yc));
    end
    
    B = a | r;
else
    B = false;
end
end