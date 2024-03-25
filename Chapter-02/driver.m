clear all; clc;

kappa = 1.0;   % isotropic and homogeneous heat conductivity

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D( n_int_xi, n_int_eta );

% mesh generation
n_en   = 4;      % quad element
n_el_x = 8;
n_el_y = 6;
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

n_eq = counter;

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
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
    end
  end
end














% EOF