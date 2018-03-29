% CREATE_AB_LINEAR
% 16-831 Fall 2016 - *Stub* Provided
% Computes the A and b matrices for the 2D linear SLAM problem
%
% Arguments: 
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
%                 obs(:,3) - x-value of landmark measurement
%                 obs(:,4) - y-value of landmark measurement
%     sigma_o - Covariance matrix corresponding to the odometry
%               measurements
%     sigma_l - Covariance matrix corresponding to the landmark
%               measurements
% Returns:
%     A       - MxN matrix
%     b       - Mx1 vector
%
function [As, b] = create_Ab_linear(odom, obs, sigma_o, sigma_l)


% Useful Constants
n_poses = size(odom, 1) + 1; % +1 for prior on the first pose
n_landmarks = max(obs(:,2));

n_odom = size(odom, 1);
n_obs  = size(obs, 1);

% Dimensions of state variables and measurements (all 2 in this case)
p_dim = 2;
l_dim = 2;
o_dim = size(odom, 2);
m_dim = size(obs(1, 3:end), 2);

% A matrix is MxN, b is Mx1
N = p_dim*n_poses + l_dim*n_landmarks;
M = o_dim*(n_odom+1) + m_dim*n_obs;     % +1 for prior on the first pose

% Initialize matrices
A = zeros(M, N);
b = zeros(M, 1);
sigma_o = sqrt(sigma_o);
sigma_l = sqrt(sigma_l);

% Add odometry and landmark measurements to A, b - including prior on first
% pose

%upper left quadrant - pose/odom
pose_block = eye(p_dim)/sigma_o;
A(1:2,1:2) = eye(2);
for i = 3:2:p_dim*n_poses
    A(i:i+1,i:i+1) = pose_block;
    A(i:i+1,i-2:i-1) = -1*pose_block;
end

%lower left quadrant - pose/landmark measurement
r_offset = o_dim*(n_odom+1);
c_offset = p_dim*(n_poses);
land_block = eye(l_dim)/sigma_l
for i = 1:n_obs
    l = obs(i,:);
    A(r_offset+2*i-1:r_offset+2*i, l(1)*2-1:l(1)*2) = -land_block;
    A(r_offset+2*i-1:r_offset+2*i, c_offset+l(2)*2-1:c_offset+l(2)*2) = land_block;
end

b = cat(1, [0;0],reshape(odom', [],1))/sigma_o(1,1);
b = cat(1, b, reshape(obs(:,3:4)',[],1)/sigma_l(1,1));


% Make A a sparse matrix 
As = sparse(A);
