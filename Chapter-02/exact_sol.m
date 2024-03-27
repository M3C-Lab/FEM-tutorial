function [val] = exact_sol(xx,yy)

val = xx*xx*(1.0-xx)*yy*(1.0-yy)*(1.0-yy);