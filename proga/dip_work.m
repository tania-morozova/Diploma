clear
clc
%%
%nobody extincts?
r = [1,1,1,1]';
b = [1,1,1,0]';
c = [0,1,1,1]';
u_min1 = [1,0];
u_max1 = [2,0];

u1_min = 1;
u1_max = 2;
u2_min = 0;
u2_max = 2;

x_0 = [10,1,1,1]';

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

P_curr = P([u1_min, u2_min]);
ok = P_curr(4) > 0 %requirement for parameters

K = @(x,u) KK(x,u,P,b,c);
K1 = @(x) K(x,[u1_min,u2_min]);
K2 = @(x) K(x,[u1_max,u2_min]);

delta_t = 10^(-5);


%%

P_curr = P([u1_min,u2_min]);
if x_0(1) > P_curr(1)
    u_curr = [u1_min,u2_min];
elseif x_0(1) < P_curr(1)
    u_curr = [u1_max,u2_min];
else disp('x_0 = P_1');
end

f(0,x_0,u_curr)
f_s = @(t,x)f_synth(t,x,u_min,u_max,f,P_curr);
%%

options1 = odeset('Events',@(t,x)events_func(t,x,P,[u2_min,u2_max]));
options2 = odeset('Events',@(t,x)events_func2(t,x,P,[u2_min,u2_max]));

sol = ode45(f_s, t_0:delta_t:t_1, x_0, options1);
solut = sol.y;      %trajectory
time = sol.x;       %time
mom_switch = numel(time);       %for display of switch moments
time_switch = time(end);
x_switch = solut(:,end);

num = 0;
while sol.ie == 1
    u_curr = [u1_min + u1_max - u_curr(1), u2_min];    
        
    sol = ode45(f_s, time(end):delta_t:t_1, solut(:,end), options1);
    
    solut = [solut, sol.y];
    time = [time, sol.x];
    mom_switch = [mom_switch, mom_switch(end) + numel(sol.x)];
    time_switch = [time_switch, time(end)];
    x_switch = [x_switch, solut(:,end)];
    num = num+1;
end

num
if sol.ie == 2
    u_curr = [u_curr(1),u2_max];
    
    sol = ode45(f_s, time(end):delta_t:t_1, solut(:,end), options2);
    solut = [solut, sol.y];
    time = [time, sol.x];
    mom_switch = [mom_switch, mom_switch(end) + numel(sol.x)];
    time_switch = [time_switch, time(end)];
    x_switch = [x_switch, solut(:,end)];
    num = num+1;
end

mom_switch = mom_switch(1:end-1);
time_switch = time_switch(1:end-1);
x_switch = x_switch(:,1:end-1);

%%
len = numel(solut(1,:));
K_ev1 = zeros(1,len);
K_ev2 = zeros(1,len);

for i = 1:len
    K_ev1(i) = K1(solut(:,i));
    K_ev2(i) = K2(solut(:,i));
end

figure
plot(time, K_ev1, 'r', time, K_ev2, 'g', time_switch, K_ev1(mom_switch), 'b*',...
    time_switch, K_ev2(mom_switch), 'b*')
xlabel('t');
ylabel('K');
legend('K_1','K_2', 'switches');
grid minor

%%
P_curr1 = P(u_min);
P_curr2 = P(u_max);

figure
ax1 = subplot(2,2,1);
plot(time, solut(1,:), time_switch, x_switch(1,:),'go',...
    [min(time) max(time)], [P_curr1(1), P_curr1(1)], 'r');
xlabel('t');
ylabel('x_1');
legend('x_1(t)', 'switches','P_1');
grid minor    

ax2 = subplot(2,2,2);
plot(time, solut(2,:), time_switch, x_switch(2,:),'go', ...
    [min(time) max(time)], ones(1,2)*P_curr1(2), 'r', ...
    [min(time) max(time)], ones(1,2)*P_curr2(2), 'r');
xlabel('t');
ylabel('x_2');
legend('x_2(t)', 'switches','P_2');
grid minor  

ax3 = subplot(2,2,3);
plot(time, solut(3,:), time_switch, x_switch(3,:),'go', ...
    [min(time) max(time)], ones(1,2)*P_curr1(3),'r');
xlabel('t');
ylabel('x_3');
legend('x_3(t)', 'switches','P_3');
grid minor  

ax4 = subplot(2,2,4);
plot(time, solut(4,:), time_switch, x_switch(4,:),'go', ...
    [min(time) max(time)], ones(1,2)*P_curr1(4), 'r', ...
    [min(time) max(time)], ones(1,2)*P_curr2(4), 'r');
xlabel('t');
ylabel('x_4');
legend('x_4(t)', 'switches','P_4');
grid minor  

linkaxes([ax4,ax3,ax2,ax1],'x'); 
