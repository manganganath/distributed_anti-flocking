function p_h = action_func(z,d,c1)

z=z/d;
p_h = zeros(size(z));

i = z<=1;
p_h(i) = -(c1*0.5*pi/d)*sin(0.5*pi*(z(i)+1));