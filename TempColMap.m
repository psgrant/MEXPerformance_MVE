function [map] = TempColMap(N)
% WATER Returns the water colormap as a set of normalized RGB colors.
%   Inputs:
%       - N: Number of colors in the colormap (optional, default = 128).
%   Output:
%       - map: Water colormap, a matrix of size N-by-3, where each row
%       represents an RGB color value.

% Check if the number of input arguments is less than 1
if nargin < 1
    N = 128;  % Default value for N
end

% Define the interpolation points for the colormap
% InterpPoints = [247,251,255
%                 222,235,247
%                 198,219,239
%                 158,202,225
%                 107,174,214
%                 66,146,198
%                 33,113,181
%                 8,81,156
%                 8,78,137]/255;

InterpPoints = [222,235,237
    254,224,210
    252,187,161
    252,146,114
    251,106,74
    239,59,44
    203,24,29
    153,0,13]/255;




% Create a vector of x values corresponding to the interpolation points
x = 1:length(InterpPoints);

% Create a vector of equally spaced x values for interpolation
xq = linspace(1, length(InterpPoints), N)';

% Perform linear interpolation to generate the colormap
map = interp1(x, InterpPoints(:, 1:3), xq);