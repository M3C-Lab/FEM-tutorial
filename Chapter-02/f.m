function [val] = f(xx, yy)

val = 2*xx^2*(2*yy - 2)*(xx - 1) + 2*xx^2*yy*(xx - 1) + 4*xx*yy*(yy - 1)^2 + 2*yy*(xx - 1)*(yy - 1)^2;

% dist = sqrt( (xx-0.5).^2 + (yy-0.5).^2 );
% if dist < 0.05
%   val = 1.0;
% else
%   val = 0.0;
% end