clc
clear
x0 = [0.2;0];
x01 = [0.2;0.2];
x02 = [0.8;0.2];
xf = [0.8;0];
t_max = 2;
interval_2 = cal_via_jerk_traj;
s = linspace(0,10,1000);
figure(1);
plot(s,mytraj(s,interval_2));
legend('x','y','v_x','v_y');title('Trajectory with respect to time(s)');
grid on ; xlabel('Time (s)') ;
figure(2);
res = mytraj(s,interval_2);
plot(res(:,1),res(:,2))
title('Trajectory in 2D phase');
grid on ; xlabel('position (m)') ; ylabel('position (m)') ;
function pos_vel = mytraj(t,interval_2)

    pos_vel=zeros(max(size(t)),4);   

    for i=1:length(t)       
        pos_vel(i,:)= calc_traj(t(i),interval_2)';  
    end

end

function pos_vel = calc_traj(t,interval_2)
    x0 = [0.2;0];
    x01 = [0.2;0.2];
    x02 = [0.8;0.2];
    xf = [0.8;0];
    t_max = 2;
    if(t>=0 && t< 2)
        pos_vel = cal_jerk_traj(t,x0,x01,t_max) ;
    elseif(t>=2 && t<8)
        pos_vel = interval_2(t-2)+ [0.2;0.2;0;0];
    else
        pos_vel = cal_jerk_traj(t-8,x02,xf,t_max) ;
    end
    fprintf('simulating process: t = %f s \n',t)
end


function res = cal_jerk_traj(t,x0,xf,t_max) 
    xd = x0 + (6*(t/t_max).^5 -15* (t/t_max).^4 + 10* (t/t_max).^3)*(xf-x0);
    vd = 1/t_max* (30*(t/t_max).^4 -60* (t/t_max).^3 + 30* (t/t_max).^2)*(xf-x0);
    res = [xd; vd];
    
end