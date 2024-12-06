function [s, z] = scaled_arg(z)
% Scaling to obtain well-behaved PadÈ approximations.
% Use DOUBLE to handle VPA calculations too
s = max(0, nextpow2(norm(double(z), inf)/4));
z = z ./ 2^s;

end