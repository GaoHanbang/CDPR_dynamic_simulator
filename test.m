syms phi theta psi
% w(expressed in body frame) = a* d_orientation
a = [1 0 -sin(theta);
     0 cos(phi) sin(phi)* cos(theta);
     0 -sin(phi) cos(phi)* cos(theta)];
simplify(inv(a));
% w(expressed in inertia frame) = a* d_orientation
b = [cos(theta) * cos(psi) -sin(psi) 0;
    cos(theta)* sin(psi) cos(psi) 0
    -sin(theta) 0 1];
inv_b = simplify(inv(b));
d = 0.3;


% [             cos(psi)/cos(theta),              sin(psi)/cos(theta), 0]
% [                       -sin(psi),                         cos(psi), 0]
% [cos(psi)*tan(theta), sin(psi)*tan(theta), 1]


plot(times,theta);
hold on 
plot(all_times,vd_generated);
title('simulated orientation w.r.t. time');
grid on ; xlabel('time (s)') ; ylabel('orientation (rad)') ;
legend('\phi','\theta','\psi');
