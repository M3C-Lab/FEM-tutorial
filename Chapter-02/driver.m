clear all; clc;

kappa = 1.0;   % isotropic and homogeneous heat conductivity

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D( n_int_xi, n_int_eta );

% mesh generation
n_en   = 4;      % quad element
n_el_x = 80;
n_el_y = 80;
n_el   = n_el_x * n_el_y;

n_np_x = n_el_x + 1;
n_np_y = n_el_y + 1;
n_np   = n_np_x * n_np_y;

hh_x = 1 / n_el_x;
hh_y = 1 / n_el_y;

% generate the coordinates of nodal points
x_coor = zeros(n_np, 1);
y_coor = zeros(n_np, 1);

for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny - 1) * n_np_x + nx;
    x_coor(index) = hh_x * (nx - 1);
    y_coor(index) = hh_y * (ny - 1);
  end
end

% generate the IEN array
IEN = zeros(n_el, n_en);
for ey = 1 : n_el_y
  for ex = 1 : n_el_x
    ee = (ey - 1) * n_el_x + ex;
    IEN(ee, 1) = (ey - 1) * n_np_x + ex;
    IEN(ee, 2) = (ey - 1) * n_np_x + ex + 1;
    IEN(ee, 3) = ey * n_np_x + ex + 1;
    IEN(ee, 4) = ey * n_np_x + ex;
  end
end

% ID array
ID = zeros(n_np, 1);

counter = 1;
for ny = 2 : n_np_y - 1
  for nx = 2 : n_np_x - 1
    index = (ny-1) * n_np_x + nx;
    ID(index) = counter;
    counter = counter + 1;
  end
end

n_eq = counter - 1;

LM = ID(IEN);

% Start the assembly
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

for ee = 1 : n_el
  k_ele = zeros(n_en, n_en);
  f_ele = zeros(n_en, 1);
  
  x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );
  
  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
    
    for aa = 1 : n_en
      Na = Quad(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      
      Na_x = ( Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
      Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
      
      f_ele(aa) = f_ele(aa) + weight(ll) * Na * f(x_l, y_l) * detJ;
      
      for bb = 1 : n_en
        [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
        Nb_x = ( Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
        
        k_ele(aa,bb) = k_ele(aa,bb) + weight(ll) * kappa * (Na_x * Nb_x + Na_y * Nb_y) * detJ;
      end % end of bb loop
    end % end of aa loop
  end % end of quadrature loop
  
  % Assembly of the global stiffness matrix and load vector
  for aa = 1 : n_en
    PP = LM(ee, aa);
    if PP > 0
      F(PP) = F(PP) + f_ele(aa);
      for bb = 1 : n_en
        QQ = LM(ee, bb);
        if QQ > 0
          K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
        end
      end
    end
  end
  
end % end of element loop

% Solve the linear system for the displacement vector
temp = K \ F;

% enforce the Dirichlet boundary nodal value
disp = zeros(n_np, 1);

for AA = 1 : n_np
  PP = ID(AA);
  if PP > 0
    disp(AA) = temp(PP);
  else
    % insert the Dirichlet value
  end
end

% error calculation
error_top = 0.0; error_bot = 0.0;

for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );
  u_ele = disp( IEN(ee, 1:n_en) );
  
  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0; u_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
      u_l = u_l + u_ele(aa) * Quad(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

    u_ext = exact_sol(x_l, y_l);
    
    error_top = error_top + weight(ll) * detJ * ( u_ext - u_l )^2;
    error_bot = error_bot + weight(ll) * detJ * u_ext^2;
  end
end

error = sqrt(error_top) / sqrt(error_bot);

% EOF