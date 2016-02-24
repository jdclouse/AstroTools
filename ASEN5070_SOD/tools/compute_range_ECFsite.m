function range = compute_range_ECFsite( inrtl_pos, ecf_site, theta )
%compute_range_ECFsite Summary of this function goes here
%   Detailed explanation goes here
% fcnPrintQueue(mfilename('fullpath'))

x = inrtl_pos(1);
y = inrtl_pos(2);
z = inrtl_pos(3);
xs = ecf_site(1);
ys = ecf_site(2);
zs = ecf_site(3);

range = sqrt(...
    (x-(xs*cos(theta)-ys*sin(theta)))*(x-(xs*cos(theta)-ys*sin(theta))) + ...
    (y-(xs*sin(theta)+ys*cos(theta)))*(y-(xs*sin(theta)+ys*cos(theta))) + ...
    (z-zs)*(z-zs));

end

