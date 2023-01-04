% This script generates the minimum jerk trajectory through the initial, final and given viapoint
% reference: https://www.jneurosci.org/content/jneuro/5/7/1688.full.pdf
% page 15
function res = cal_via_jerk_traj

syms t1 t
t1 = sym('t1', 'real'); 
tf = 6;
% initial point [0;0];
x0 = 0.0;
y0 = 0.0;

x1 = 0.3 ; % viapoint
y1 = 0.1 ;

xf = 0.6 ; % final point
yf = 0.0 ;

%% To solve the viapoint time t1 such that x(t1) = x1; y(t2)= 
A2 = expand(xf*(120*t1^5 - 300*t1^4 + 200*t1^3) - 20*x1) ;
A22 = expand(yf*(120*t1^5 - 300*t1^4 + 200*t1^3) - 20*y1);
A1 = expand(xf*(300*t1^5 - 1200*t1^4 + 1600*t1^3) + t1^2*(-720*xf + 120*x1) - x1*(300*t1 - 200)) ;
A11 = expand(yf*(300*t1^5 - 1200*t1^4 + 1600*t1^3) + t1^2*(-720*yf + 120*y1) - y1*(300*t1 - 200)) ;
A3 = 60*t1^7 - 210*t1^6 + 240*t1^5 - 90*t1^4 ;
A4 = 60*t1^3 - 30*t1^2 - 30*t1^4 ;
Result = expand(A3*A2*A2 + A3*A22*A22 + (t1^3)*A1*A2*A4 + (t1^3)*A11*A22*A4) ;
%numerical solution
solvet = vpasolve(Result, t1);
[size_solvet,~] = size(solvet);
% solution should be 0<t1<1
for i = 1:size_solvet
    if (solvet(i)>0 && solvet(i)<1)
        sol = double(solvet(i));
        break;
    end
end
% fprintf('solution of t1 is %f\n', sol);
%% 
sA2 =  xf*(120*sol^5 - 300*sol^4 + 200*sol^3) - 20*x1 ;
sA22 = yf*(120*sol^5 - 300*sol^4 + 200*sol^3) - 20*y1;
sA1 =  xf*(300*sol^5 - 1200*sol^4 + 1600*sol^3) + sol^2*(-720*xf + 120*x1) - x1*(300*sol - 200) ;
sA11 = yf*(300*sol^5 - 1200*sol^4 + 1600*sol^3) + sol^2*(-720*yf + 120*y1) - y1*(300*sol - 200) ;
sA3 = 60*sol^7 - 210*sol^6 + 240*sol^5 - 90*sol^4 ;
sA4 = 60*sol^3 - 30*sol^2 - 30*sol^4 ;

c1x =  sA1/ (tf^5*sol^2* (1-sol)^5);
pi1x = sA2/ (tf^5*sol^5* (1-sol)^5);
c1y =  sA11/ (tf^5*sol^2* (1-sol)^5);
pi1y = sA22/ (tf^5*sol^5* (1-sol)^5);

%% express the trajectory
syms x_(t) y_(t) xs(t) ys(t) x(t) y(t) vx(t) vy(t)
x_(t) = tf^5/720*(pi1x*(sol^4*(15*(t/tf)^4 - 30*(t/tf)^3)+ sol^3*(80*(t/tf)^3 - 30* (t/tf)^4) - 60* (t/tf)^3* sol^2 + 30* (t/tf)^4* sol - 6* (t/tf)^5) + c1x*(15*(t/tf)^4- 10* (t/tf)^3 - 6 *(t/tf)^5)) +x0  ;
y_(t) = tf^5/720*(pi1y*(sol^4*(15*(t/tf)^4 - 30*(t/tf)^3)+ sol^3*(80*(t/tf)^3 - 30* (t/tf)^4) - 60* (t/tf)^3* sol^2 + 30* (t/tf)^4* sol - 6* (t/tf)^5) + c1y*(15*(t/tf)^4- 10* (t/tf)^3 - 6 *(t/tf)^5)) +y0  ;
xs(t) = pi1x* tf^5* (t/tf - sol)^5/120;
ys(t) = pi1y* tf^5* (t/tf - sol)^5/120;
x(t) = piecewise(t/tf< sol,x_(t), t/tf> sol,x_(t)+xs(t));
y(t) = piecewise(t/tf< sol,y_(t), t/tf> sol,y_(t)+ys(t));
vx(t) = diff(x(t));
vy(t) = diff(y(t));

res(t) = [x(t);y(t); vx(t);vy(t)];
%% to plot the results
% t = linspace(0,tf,100);
% plot(x(t),y(t))
% axis equal
% % plot(t,x(t))
% % hold on
% % plot(t,vx(t))
end