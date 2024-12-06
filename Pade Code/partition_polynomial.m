function [s, r] = partition_polynomial(d)
% Assume there is always at least one term
s = max(floor(sqrt(floor(d / 2))), 1);
r = floor(d / s);
end