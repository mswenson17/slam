syms l_y l_x r_x r_y theta d  x y
assume(l_y, 'real')
assume(l_x, 'real')
assume(r_y, 'real')
assume(r_x, 'real')

h_m = [atan(l_y - r_y,l_x-r_x);
       sqrt((l_y - r_y)^2+(l_x-r_x)^2)]

H_m = simplify(jacobian(h_m, [r_x, r_y, l_x, l_y]))
H_m_simp = subs(subs(H_m,l_y - r_y, y),l_x - r_x, x)
    
