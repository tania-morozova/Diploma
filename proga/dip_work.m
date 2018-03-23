clear
clc
%%
%nobody extincts?
r = [1,1,1,1]';
b = [1,1,1,0]';
c = [0,1,1,1]';
u_min1 = [1,0];
u_max1 = [2,0];

u_min = [1,0];
u_max = [2,2];

x_0 = [10,1,3,1]';

t_0 = 0;
t_1 = 10;

%%

f = @(t,x,u) [x(1).*(r(1) + u(1) - b(1).*x(2)); ...
            x(2).*(-r(2) - b(2).*x(3) + c(2).*x(1));...
            x(3).*(-r(3) + u(2) - b(3).*x(4) + c(3).*x(2));...
            x(4).*(-r(4) + c(4).*x(3))];
    
P = @(u) [(r(2)*c(4) + b(2)*r(4))/(c(4)*c(2)); ...
     (r(1) + u(1))/b(1); ...
     r(4)/c(4); ...
     (c(3)*(r(1) + u(1))+(u(2) - r(3))*b(1))/(b(1)*b(3))];

P_min = P(u_min);
ok = P_min(4) > 0 %requirement for parameters

K = @(x,u) KK(x,u,P,b,c);
K1 = @(x) K(x,u_min);
K2 = @(x) K(x,[u_max(1),u_min(2)]);
K3 = @(x) K(x,[u_min(1),u_max(2)]);
K4 = @(x) K(x,u_max);

delta_t = 10^(-5);

%%

P_curr = P(u_min);
f_s = @(t,x)f_synth(t,x,u_min,u_max,f,P,r,b,c);

options = odeset('Events',@(t,x)events_func(t,x,P_curr));

sol = ode45(@(t,x) f_s(t,x), t_0:delta_t:t_1, x_0, options);
solut = sol.y;      %trajectory
time = sol.x;       
      %for display of switch moments
mom_switch1 = [];
mom_switch2 = [];
time_switch1 = [];
time_switch2 = [];
x_switch1 = [];
x_switch2 = [];
if sol.ie == 1
    mom_switch1 = numel(time); 
    time_switch1 = time(end);
    x_switch1 = solut(:,end);
else
    mom_switch2 = numel(time);
    time_switch2 = time(end);
    x_switch2 = solut(:,end);
end

num = 0;
while sol.ie 
    if t_1-time(end)<delta_t
        break;
    end
    sol = ode45(@(t,x) f_s(t,x), time(end):delta_t:t_1, solut(:,end), options);
    
    solut = [solut, sol.y];
    time = [time, sol.x];
    if sol.ie == 1
        if isempty(mom_switch1)
            mom_temp = 1;
        else
            mom_temp = mom_switch1(end);
        end
        
        mom_switch1 = [mom_switch1, mom_temp + numel(sol.x)];
        time_switch1 = [time_switch1, time(end)];
        x_switch1 = [x_switch1, solut(:,end)];
    else
        if isempty(mom_switch2)
            mom_temp = 1;
        else
            mom_temp = mom_switch2(end);
        end
        mom_switch2 = [mom_switch2, mom_temp + numel(sol.x)];
        time_switch2 = [time_switch2, time(end)];
        x_switch2 = [x_switch2, solut(:,end)];
    end
    num = num+1;
end

mom_switch1 = mom_switch1(1:end-1);
time_switch1 = time_switch1(1:end-1);
x_switch1 = x_switch1(:,1:end-1);
mom_switch2 = mom_switch2(1:end-1);
time_switch2 = time_switch2(1:end-1);
x_switch2 = x_switch2(:,1:end-1);

%%
len = numel(solut(1,:));
K_ev1 = zeros(1,len);
K_ev2 = zeros(1,len);
K_ev3 = zeros(1,len);
K_ev4 = zeros(1,len);

K1_lim = K1(P(u_min));
K2_lim = K2(P([u_max(1),u_min(2)]));
K3_lim = K3(P([u_min(1),u_max(2)]));
K4_lim = K4(P(u_max));

for i = 1:len
    K_ev1(i) = K1(solut(:,i));
    K_ev2(i) = K2(solut(:,i));
    K_ev3(i) = K3(solut(:,i));
    K_ev4(i) = K4(solut(:,i));
end

figure
hold on
plot(time, K_ev1, 'r', time_switch1, K_ev1(mom_switch1), 'c*', ...
    time_switch2, K_ev1(mom_switch2), 'm*');
plot(time, K_ev2, 'b', time_switch1, K_ev2(mom_switch1), 'c*', ...
    time_switch2, K_ev2(mom_switch2), 'm*');
plot(time, K_ev3, 'g', time_switch1, K_ev3(mom_switch1), 'c*', ...
    time_switch2, K_ev3(mom_switch2), 'm*');
plot(time, K_ev4, 'k', time_switch1, K_ev4(mom_switch1), 'c*', ...
    time_switch2, K_ev4(mom_switch2), 'm*');
se = [time(1) time(end)];
plot(se, [K1_lim K1_lim], 'r--', se, [K2_lim K2_lim], 'b--', ...
    se, [K3_lim K3_lim], 'g--', se, [K4_lim K4_lim], 'k--');
xlabel('t');
ylabel('K');
legend('K_1','sw1','sw2','K_2','sw1','sw2','K_3','sw1','sw2','K_4','sw1','sw2');
grid minor
hold off

%%
P_curr1 = P(u_min);
P_curr2 = P(u_max);

figure
ax1 = subplot(2,2,1);
plot(time, solut(1,:), time_switch1, x_switch1(1,:),'co',...
    time_switch2, x_switch2(1,:),'mo', ...
    [min(time) max(time)], [P_curr1(1), P_curr1(1)], 'r');
xlabel('t');
ylabel('x_1');
legend('x_1(t)', 'switches1','switches2','P_1');
grid minor    

ax2 = subplot(2,2,2);
plot(time, solut(2,:), time_switch1, x_switch1(2,:),'co',...
    time_switch2, x_switch2(2,:),'mo')%, ...
    %[min(time) max(time)], ones(1,2)*P_curr1(2), 'r', ...
    %[min(time) max(time)], ones(1,2)*P_curr2(2), 'r');
xlabel('t');
ylabel('x_2');
legend('x_2(t)', 'switches1','switches2'),%'P_2');
grid minor  

ax3 = subplot(2,2,3);
plot(time, solut(3,:), time_switch1, x_switch1(3,:),'co', ...
    time_switch2, x_switch2(3,:),'mo',...
    [min(time) max(time)], ones(1,2)*P_curr1(3),'r');
xlabel('t');
ylabel('x_3');
legend('x_3(t)', 'switches1','switches2','P_3');
grid minor  

ax4 = subplot(2,2,4);
plot(time, solut(4,:), time_switch1, x_switch1(4,:),'co',...
    time_switch2, x_switch2(4,:),'mo')%, ...
    %[min(time) max(time)], ones(1,2)*P_curr1(4), 'r', ...
    %[min(time) max(time)], ones(1,2)*P_curr2(4), 'r');
xlabel('t');
ylabel('x_4');
legend('x_4(t)', 'switches1','switches2')%,'P_4');
grid minor  

linkaxes([ax4,ax3,ax2,ax1],'x'); 
