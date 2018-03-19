syms x y theta alpha d p P
syms sigma_x2 sigma_y2 sigma_alpha2
p = [x;y;theta]

x_pre = p + [d*cos(theta); d*sin(theta); alpha];
latex(x_pre)
J_m = jacobian(x_pre, [x,y,theta])

G_t = J_m 

pose_cov = diag([0.02^2, 0.02^2, 0.1^2]);
control_cov = diag([sigma_x2, sigma_y2, sigma_alpha2])
P_pre =  G_t*pose_cov*G_t' + control_cov
    %J_m = [1, 0 , -d*sin(theta);
           %0, 1 , d*cos(theta);
           %0, 0   1];


syms r n_r beta n_beta l_x l_y
l = p(1:2) + [(r+n_r)*cos(beta+n_beta+theta);
              (r+n_r)*sin(beta+n_beta+theta)]



syms delta_x delta_y q
delta_x1 = l_x-x
delta_y1 = l_y-y

r_est = sqrt((delta_x1)^2 + (delta_y1)^2);
beta_est = atan2(delta_y1, delta_x1)-theta;

pose1 = [r_est; beta_est]

pose = subs(pose1, l_x-x, delta_x);
pose = subs(pose, l_y-y, delta_y)
pose(2) = -theta + atan2(delta_y,delta_x)
pose = subs(pose, delta_x^2 + delta_y^2, q)

H_p = (1/q)*[-sqrt(q)*delta_x -sqrt(q)*delta_y  0 sqrt(q)*delta_x sqrt(q)*delta_y;
              delta_y -delta_x  -q -delta_y delta_x]


H_l = jacobian(l, [r,beta])
H_l = subs(H_l, n_r, 1) 
H_l = subs(H_l, n_beta, 0) 
