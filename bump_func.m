function p_h = bump_func(z,h)

p_h = zeros(size(z));

p_h(z<h) = 1;

i = z>=h & z<=1;
p_h(i) = 0.5*(1+cos(pi*(z(i)-h)/(1-h)));