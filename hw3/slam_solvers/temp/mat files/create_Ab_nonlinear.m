% CREATE_AB_NONLINEAR
% 16-831 Fall 2016 - *Stub* Provided
% Computes the A and b matrices for the 2D nonlinear SLAM problem
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
%     A       - MxN matrix
%     b       - Mx1 vector
%
function [As, b] = create_Ab_nonlinear(x, odom, obs, sigma_o, sigma_l)
n_poses = size(odom, 1) + 1;                % +1 for prior on the first pose
n_landmarks = max(obs(:,2));

n_odom = size(odom, 1);
n_obs  = size(obs, 1);

% Dimensions of state variables and measurements (all 2 in this case)
p_dim = 2;                                  % pose dimension
l_dim = 2;                                  % landmark dimension
o_dim = size(odom, 2);                      % odometry dimension
m_dim = size(obs(1, 3:end), 2);             % landmark measurement dimension

% A matrix is MxN, b is Mx1
N = p_dim*n_poses + l_dim*n_landmarks;
M = o_dim*(n_odom+1) + m_dim*n_obs;         % +1 for prior on the first pose

% Initialize matrices
A = zeros(M, N);
b = zeros(M, 1);
sigma_o = sqrt(sigma_o);
sigma_l = sqrt(sigma_l);

% Add odometry and landmark measurements to A, b - including prior on first
% pose
hx = zeros(2*(n_poses+n_obs),1);

%upper left quadrant - pose/odom
pose_block = eye(p_dim)/sigma_o;
A(1:2,1:2) = eye(2);
for i = 3:2:p_dim*n_poses
    A(i:i+1,i:i+1) = pose_block;
    A(i:i+1,i-2:i-1) = -1*pose_block;
    rx1 = x(i-2);
    ry1 = x(i-1);
    rx2 = x(i);
    ry2 = x(i+1);
    hx(i:i+1) = meas_odom(rx1, ry1, rx2, ry2)/sigma_o(1,1);
end

%lower left quadrant - pose/landmark measurement
r_offset = o_dim*(n_odom+1);
c_offset = p_dim*(n_poses);

for i = 1:n_obs
    x_idx = obs(i,1)*2-1;
    l_idx = obs(i,2)*2+c_offset-1;
    
    rx = x(x_idx);
    ry = x(x_idx+1);
    lx = x(l_idx);
    ly = x(l_idx+1);

    J = meas_landmark_jacobian(rx,ry,lx,ly)/sigma_l(1);
    J1 = J(:,1:2); 
    J2 = J(:,3:4); 
    
    m = meas_landmark(rx, ry, lx, ly);
    m(1) = wrapToPi(m(1));
    hx(r_offset+2*i-1:r_offset+2*i) = m/sigma_l(1,1);

    l = obs(i,:);
    A(r_offset+2*i-1:r_offset+2*i, l(1)*2-1:l(1)*2) = J1;
    A(r_offset+2*i-1:r_offset+2*i, c_offset+l(2)*2-1:c_offset+l(2)*2) = J2;
end

z = cat(1, [0;0],reshape(odom', [],1))/sigma_o(1,1);
z = cat(1, z, reshape(obs(:,3:4)',[],1)/sigma_l(1,1));

b = z - hx;

% Make A a sparse matrix 
As = sparse(A);
