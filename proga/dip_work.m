clear
clc
%%
%nobody extincts?
r = [1,1,1,1]';
b = [1,1,1,0]';
c = [0,1,1,1]';
u_min = 1;
u_max = 2;

x_0 = [2,2.39,1,1.5]';

t_0 = 0;
t_1 = 15;

%%

ok = u_min - (r(3)*b(1) - c(3)*r(1))/c(3) > 0 %requirement for parameters
f = @(t,x,u) [x(1).*(r(1) + u - b(1).*x(2)); ...
            x(2).*(-r(2) - b(2).*x(3) + c(2).*x(1));...
            x(3).*(-r(3) - b(3).*x(4) + c(3).*x(2));...
            x(4).*(-r(4) + c(4).*x(3))];
    
P = @(u) [(r(2)*c(4) + b(2)*r(4))/(c(4)*c(2)); ...
     (r(1) + u)/b(1); ...
     r(4)/c(4); ...
     (c(3)*(r(1) + u) - r(3)*b(1))/(b(1)*b(3))];

K = @(x,u) KK(x,u,P,b,c);
K1 = @(x) K(x,u_min);
K2 = @(x) K(x,u_max);

delta_t = 10^(-5);

%%
P_min = P(u_min);
P_max = P(u_max);
P_curr = P(u_min);
if x_0(1) > P_curr(1)
    u_curr = u_min;
elseif x_0(1) < P_curr(1)
    u_curr = u_max;
else disp('x_0 = P_1');
end
%x_0(3:4) = [P_min(3); 1.5];

f_s = @(t,x)f_synth(t,x,u_min,u_max,f,P,r,b);

options = odeset('Events',@(t,x)events_func(t,x,P_curr));

sol = ode45(@(t,x) f_s(t,x), t_0:delta_t:t_1, x_0, options);
solut = sol.y;      %trajectory
time = sol.x;       %time
mom_switch = numel(time);       %for display of switch moments
time_switch = time(end);
x_switch = solut(:,end);

num = 0;
while sol.ie
    %u_curr = u_max - u_curr + u_min;    
       
    if t_1-time(end)<delta_t
        break;
    end
    sol = ode45(@(t,x) f_s(t,x), time(end):delta_t:t_1, solut(:,end), options);
    
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

K1_lim = K1(P_min);
K2_lim = K2(P_max);

%%
figure
plot(time, K_ev1, 'r', time, K_ev2, 'g', time_switch, K_ev1(mom_switch), 'b*',...
    time_switch, K_ev2(mom_switch), 'b*', [min(time) max(time)], [K1_lim, K1_lim], 'm',...
    [min(time) max(time)], [K2_lim, K2_lim], 'm')
xlabel('t');
ylabel('K');
legend('K_1','K_2', 'switches');
grid minor

%%
P_min = P(u_min);
P_max = P(u_max);

figure
ax1 = subplot(3,2,1);
plot(time, solut(1,:), time_switch, x_switch(1,:),'go',...
    [min(time) max(time)], [P_min(1), P_min(1)], 'r');
xlabel('t');
ylabel('x_1');
legend('x_1(t)', 'switches');%,'P_1');
grid minor    

ax2 = subplot(3,2,2);
plot(time, solut(2,:), time_switch, x_switch(2,:),'go', ...
    [min(time) max(time)], ones(1,2)*P_min(2), 'r', ...
    [min(time) max(time)], ones(1,2)*P_max(2), 'r');
xlabel('t');
ylabel('x_2');
legend('x_2(t)', 'switches','P_2');
grid minor  

ax3 = subplot(3,2,3);
plot(time, solut(3,:), time_switch, x_switch(3,:),'go', ...
    [min(time) max(time)], ones(1,2)*P_min(3),'r');
xlabel('t');
ylabel('x_3');
legend('x_3(t)', 'switches','P_3');
grid minor  

ax4 = subplot(3,2,4);
plot(time, solut(4,:), time_switch, x_switch(4,:),'go', ...
    [min(time) max(time)], ones(1,2)*P_min(4), 'r', ...
    [min(time) max(time)], ones(1,2)*P_max(4), 'r');
xlabel('t');
ylabel('x_4');
legend('x_4(t)', 'switches','P_4');
grid minor  

ax5 = subplot(3,2,5);
plot(solut(2,mom_switch(end):end), solut(4,mom_switch(end):end), x_switch(2,:), x_switch(4,:), 'go',...
    [P_min(2) P_max(2)], [P_min(4) P_max(4)], 'r')

xlabel('x_2');
ylabel('x_4');
legend('x_4(x_2)', 'switches','P_2 P_4');
grid minor

linkaxes([ax4,ax3,ax2,ax1],'x'); 

%%
figure
hold on
plot3(time(mom_switch(end):end), solut(2,mom_switch(end):end), solut(4,mom_switch(end):end), 'b')

[T,X2] = meshgrid(time,P_min(2):P_max(2));

P4_fun = (c(3)*X2 - r(3))/b(3);
surf(T,X2,P4_fun, 'FaceAlpha',0.3)

xlabel('t');
ylabel('x_2');
zlabel('x_4');
hold off

