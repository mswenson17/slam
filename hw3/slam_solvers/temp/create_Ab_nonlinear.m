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
function [As, b] = create_Ab_nonlinear(x, odom, obs, sigma_o, sigma_l)
%% Extract useful constants which you may wish to use
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

%% Initialize matrices
A = zeros(M, N);
b = zeros(M, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Your code goes here %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ho = [-1  0  1  0;
      0  -1  0  1];
  
% Hm_l = [1  0;
%         0  1];
%     
% Hm_r = [-1  0;
%          0 -1];

%Calculating square-root inverse of sigma_landmark
covs = diag(sigma_l);
sqrt_covs = sqrt(covs);
reciprocal_sqrt_covs = 1./sqrt_covs;
sigma_landmark_factor = diag(reciprocal_sqrt_covs);

%Calculating square-root inverse of sigma_odom
covs = diag(sigma_o);
sqrt_covs = sqrt(covs);
reciprocal_sqrt_covs = 1./sqrt_covs;
sigma_odom_factor = diag(reciprocal_sqrt_covs);

%Filling vector 'b'
j=1;
for i=3:2:(o_dim*(n_odom+1))
    
    state_var = x(i-2:i-2+4);
    
    h_xi0 = meas_odom(state_var(1), state_var(2), state_var(3), state_var(4));
    
    z_i = [odom(j,1) ; odom(j,2)];
    b_i = sigma_odom_factor * (z_i - h_xi0);

    b(i:i+1, 1) = b_i;
    
    j = j+1;    
end

j=1;
for i=(o_dim*(n_odom+1))+1:2:size(b,1)
    
    id_p = obs(j,1);
    id_l = obs(j,2);
    
    l_ind = (p_dim*n_poses) + ((id_l-1)*2) + 1;
    p_ind = ((id_p-1)*2) + 1;
    
    state_var_l = x(l_ind:l_ind+1);
    state_var_p = x(p_ind:p_ind+1);
    
    h_xi0 = meas_landmark(state_var_p(1), state_var_p(2), state_var_l(1), state_var_l(2));
    
    z_i = [obs(j,3) ; obs(j,4)];
    diff_z_h = [ wrapToPi(z_i(1) - h_xi0(1)) ; z_i(2) - h_xi0(2) ];
    b_i = sigma_landmark_factor * (diff_z_h);
    
%     b(i, 1) = obs(j,3);
%     b(i+1, 1) = obs(j,4);
    b(i:i+1, 1) = b_i;

    j = j+1;
end

%Filling vector 'A'
A(1:2, 1:2) = sigma_odom_factor*[1 0 ; 0 1];
for i=1:2:(p_dim*n_poses)-2
    A(i+2:i+1+2, i:i+3) = sigma_odom_factor*Ho;
end

j = (o_dim*(n_odom+1))+1;
for i=1:size(obs, 1)
    id_p = obs(i,1);
    id_l = obs(i,2);
    
    l_ind = (p_dim*n_poses) + ((id_l-1)*2) + 1;
    p_ind = ((id_p-1)*2) + 1;
    
    state_var_l = x(l_ind:l_ind+1);
    state_var_p = x(p_ind:p_ind+1);
    
    H = meas_landmark_jacobian(state_var_p(1), state_var_p(2), state_var_l(1), state_var_l(2));
    Hm_l = H(:, 3:4);
    Hm_r = H(:, 1:2);
    
    A(j:j+1, l_ind:l_ind+1) = sigma_landmark_factor*Hm_l;
    A(j:j+1, p_ind:p_ind+1) = sigma_landmark_factor*Hm_r;
    
    j=j+2;
end

%% Make A a sparse matrix 
As = sparse(A);
