% ERROR_NONLINEAR
% 16-831 Fall 2016 - *Stub* Provided
% Computes the total error of all measurements (odometry and landmark)
% given the current state estimate
%
% Arguments: 
%     x       - Current estimate of the state vector
%     odom    - Matrix that contains the odometry measurements
%               between consecutive poses. Each row corresponds to
%               a measurement. 
%                 odom(:,1) - x-value of odometry measurement
%                 odom(:,2) - y-value of odometry measurement
%     obs     - Matrix that contains the landmark measurements and
%               relevant information. Each row corresponds to a
%               measurement.
%                 obs(:,1) - idx of pose at which measurement was 
%                   made
%                 obs(:,2) - idx of landmark being observed
%                 obs(:,3) - theta of landmark measurement
%                 obs(:,4) - d of landmark measurement
%     sigma_o - Covariance matrix corresponding to the odometry
%               measurements
%     sigma_l - Covariance matrix corresponding to the landmark
%               measurements
% Returns:
%     err     - total error of all measurements
%
function err = error_nonlinear(x, odom, obs, sigma_odom, sigma_landmark)
% Extract useful constants 
n_poses = size(odom, 1) + 1;               % +1 for prior on the first pose
n_landmarks = max(obs(:,2));

n_odom = size(odom, 1);
n_obs  = size(obs, 1);

% Dimensions of state variables and measurements (all 2 in this case)
p_dim = 2;                                  % pose dimension
l_dim = 2;                                  % landmark dimension
o_dim = size(odom, 2);                      % odometry dimension
m_dim = size(obs(1, 3:end), 2);    % landmark measurement dimension

% A matrix is MxN, b is Nx1
N = p_dim*n_poses + l_dim*n_landmarks;
M = o_dim*(n_odom+1) + m_dim*n_obs;         % +1 for prior on the first pose

% Initialize error
err = 0;
hx = zeros(2*(n_poses+n_obs),1);
%size(hx)
sigma_o = sqrt(sigma_odom(1));
sigma_l = sqrt(sigma_landmark(1));

%calculate h*x-z
%hx is the current state vector run through the actual odom finding equations
for i = 1:n_poses-1
    rx1 = x(2*i-1);
    ry1 = x(2*i);
    rx2 = x(2*i+1);
    ry2 = x(2*i+2);
    hx(2+2*i-1:2+2*i) = meas_odom(rx1, ry1, rx2, ry2)/sigma_o;
end

l_offset = p_dim*(n_poses);

for i = 1:n_obs
    %lookup r and L x and y coords based on observation data
    x_idx = obs(i,1)*2-1;
    l_idx = obs(i,2)*2+l_offset-1;
    
    rx = x(x_idx);
    ry = x(x_idx+1);
    lx = x(l_idx);
    ly = x(l_idx+1);
    
    %h = meas_landmark(rx1, ry1, rx2, ry2);
    m = meas_landmark(rx, ry, lx, ly);
    m(1) = wrapToPi(m(1));
    hx(l_offset+2*i-1:l_offset+2*i) = m/sigma_l;
end

%z = b from create_Ab
z = cat(1, [0;0],reshape(odom', [],1))/sigma_o;
z = cat(1, z, reshape(obs(:,3:4)',[],1))/sigma_l;

for i = 1:2:o_dim*n_obs
   z(l_offset+i) = z(l_offset+i);
end
%z
%zsize=size(z)
%hxsize=size(hx)


err = sum((hx-z).^2);
