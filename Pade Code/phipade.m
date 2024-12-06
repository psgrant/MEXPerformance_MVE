function varargout = phipade(z, k, varargin)
% PHIPADE - Evaluate phi functions using diagonal PadÈ approximations.
%
% SYNOPSIS:
%    phi_k                     = phipade(z, k);
%    phi_k                     = phipade(z, k, d);
%   [phi_1, phi_2, ..., phi_k] = phipade(...);
%
% DESCRIPTION:
%   This function evaluates phi functions needed in exponential
%   integrators using diagonal PadÈ approximants.
%   We define the phi functions according to the integral representation
%
%      \phi_k(z) = \frac{1}{(k - 1)!} \int_0^1 e^{z (1-x)} x^{k-1} dx
%
%   for k=1, 2, ...
%
% PARAMETERS:
%   z - Evaluation point.  Assumed to be one of
%         - 1D vector, treated as the main diagonal of a diagonal matrix
%         - sparse diagonal matrix
%         - full or sparse matrix
%   k - Which phi function(s) to evaluate.
%       Index (integer) of the (highest) phi function needed.
%   d - Degree of diagonal PadÈ approximant.  OPTIONAL.
%       Default value: d = 7.
%
% RETURNS:
%    phi_k                     =      \phi_k(z)
%   [phi_1, phi_2, ..., phi_k] = DEAL(\phi_1(z), \phi_2(z), ..., \phi_k(z))
%
% NOTES:
%   When computing more than one phi function, it is the caller's
%   responsibility to provide enough output arguments to hold all of the
%   \phi_k function values.
%
%   For efficiency reasons, PHIPADE caches recently computed function
%   values.  The caching behaviour is contingent on the WANTCACHE
%   function and may be toggled on or off as needed.
%
% SEE ALSO:
%   WANTCACHE.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.17 $  $Date: 2005/10/12 16:28:00 $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local user configuration.
CACHE_BND = 4;          % upper bound on number of cache entries
PADE_DEGR = 7;          % degree of (d,d)-PadÈ approximant
atol = 1.0e-12;         % FLTEQ absolute error tolerance
rtol = 5.0e-13;         % FLTEQ relative error tolerance

% Uncomment to remove WANTCACHE function dependency (ie. make PHIPADE
% run in a self contained environment).
wantcache = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No user changeable parameters below this line.

% Note: Short circuiting operators are critical in this statement.  In
%       particular we need the NUMEL check before the k==NARGOUT check
%       because the operators only accept scalar (logical) operands.
%
%       We furthermore note the inconsistency of NARGOUT.  NARGOUT is a
%       peculiar function and statements of the form
%
%               nargout (op) <some number>
%
%       generates errors due to too many input args to NARGOUT.
%       However, such statements *ARE* allowed within IF statements...

arg_ok = (0 == nargout) || (numel(k) == nargout) || (k == nargout);
if ~arg_ok,
    error('phipade:nargout', ...
        'Inconsistent number of output arguments');
end

% support the
%        [phi{1:k}] = phipade(z, k);
% syntax
if (nargout > 1) && (numel(k) == 1), k = 1:k; end

if (nargin > 2) && isnumeric(varargin{1}),
    d = varargin{1};
else
    d = PADE_DEGR;
end

% treat vectors as (sparse) diagonal matrices
if sum(size(z) > 1) == 1,
    n = length(z);
    z = spdiags(z(:), 0, n, n);
end

if wantcache,
    % main cache data structure
    % see *_cache() functions for operational detail...
    persistent phi_cache;
    
    idx = find_cache(phi_cache, z, max(k), atol, rtol);
    if idx < 1,
        [idx, phi_cache] = put_cache(phi_cache, z, max(k), idx, d);
    end
    
    nk = numel(k);
    [phi_cache, varargout{1:nk}] = get_cache(phi_cache, idx, k, CACHE_BND);
else
    pv = eval_pade(z, max(k), d);
    [varargout{1:numel(k)}] = deal(pv{1}{k});
end

end