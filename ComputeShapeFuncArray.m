
function [ShapeFuncArray] = ComputeShapeFuncArray()
% Compute the values for the Ni terms for the bilinear triangular prism
% shape function. The function generates the ShapeFuncArray matrix, where
% each row represents the Ni values for a specific point in the prism. The
% ShapeFuncArray is then used to interpolate values at different locations
% within the prism.

% Define the local coordinate array for the bilinear triangular prism
LocalCoordsArray = [5/12 1/6 -0.5;
    5/12 5/12 -0.5;
    1/6 5/12 -0.5;
    5/12 1/6 0.5;
    5/12 5/12 0.5;
    1/6 5/12 0.5;
    5/24 5/24 0;
    7/12 5/24 0;
    5/24 7/12 0];

% Extract the x, y, and z coordinates from the local coordinate array
xA = LocalCoordsArray(:, 1);
yA = LocalCoordsArray(:, 2);
zA = LocalCoordsArray(:, 3);

% Initialize the ShapeFuncArray matrix with zeros
ShapeFuncArray = zeros(length(LocalCoordsArray), 6);

% Plotting setup (optional)
% clf
% hold on

% Iterate over each point in the local coordinate array
for i = 1:length(LocalCoordsArray)
    
    % Extract the x, y, and z coordinates for the current point
    x = xA(i);
    y = yA(i);
    z = zA(i);
    
    % Compute the Ni values for the current point and store them in the
    % ShapeFuncArray
    for N = 1:6

        switch N
            case 1
                Q = (1 - x - y) .* (1 - z);
            case 2
                Q = x .* (1 - z);
            case 3
                Q = y .* (1 - z);
            case 4
                Q = (1 - x - y) .* (1 + z);
            case 5
                Q = x .* (1 + z);
            case 6
                Q = y .* (1 + z);
            otherwise
                error('Something is messed up')
        end
        
        ShapeFuncArray(i, N) = Q;
    end

end
ShapeFuncArray = ShapeFuncArray';