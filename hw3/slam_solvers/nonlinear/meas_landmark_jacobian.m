% MEAS_LANDMARK_JACOBIAN
% 16-831 Fall 2016 - *Stub* Provided
% Compute the Jacobian of the measurement function
%
% Arguments: 
%     rx    - robot's x position
%     ry    - robot's y position
%     lx    - landmark's x position
%     ly    - landmark's y position
%
% Returns:
%     H     - Jacobian of the measurement fuction
%
function H = meas_landmark_jacobian(rx, ry, lx, ly)

x = (lx-rx);
y = (ly-ry);
H = [        y/(x^2 + y^2),       -x/(x^2 + y^2),      -y/(x^2 + y^2),       x/(x^2 + y^2);
      -x/(x^2 + y^2)^(1/2), -y/(x^2 + y^2)^(1/2), x/(x^2 + y^2)^(1/2), y/(x^2 + y^2)^(1/2)];

