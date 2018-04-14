clear
clc
%% x1=p1, x3=p3
r = [1,1,1,1]';
b = [1,1,1,0]';
c = [0,1,1,1]';
u_min = [1,1.5];
u_max = [2,2];

x_0 = [2,5,1,5]';

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

 
V1 = P(u_min);
V2 = P([u_max(1),u_min(2)]);
V3 = P([u_min(1),u_max(2)]);
V4 = P(u_max); 
ok = V1(4) > 0 %requirement for parameters

K = @(x,u) KK(x,u,P,b,c);
K1 = @(x) K(x,u_min);
K2 = @(x) K(x,[u_max(1),u_min(2)]);
K3 = @(x) K(x,[u_min(1),u_max(2)]);
K4 = @(x) K(x,u_max);

delta_t = 10^(-5);

%%

P_curr = P(u_min);
f_s = @(t,x)f_synth(t,x,u_min,u_max,f,P,r,b,c);
f_s(0,x_0)
%%

options = odeset('Events',@(t,x)events_func(t,x,P_curr), 'MaxStep', 1e-2);

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
%for display
if isempty(mom_switch1)
    mom_switch1 = 1;
    time_switch1 = t_0;
    x_switch1 = x_0;
end
if isempty(mom_switch2)
    mom_switch2 = 1;
    time_switch2 = t_0;
    x_switch2 = x_0;
end
%%
len = numel(solut(1,:));
K_ev1 = zeros(1,len);
K_ev2 = zeros(1,len);
K_ev3 = zeros(1,len);
K_ev4 = zeros(1,len);

K1_lim = K1(V1);
K2_lim = K2(V2);
K3_lim = K3(V3);
K4_lim = K4(V4);

for i = 1:len
    K_ev1(i) = K1(solut(:,i));
    K_ev2(i) = K2(solut(:,i));
    K_ev3(i) = K3(solut(:,i));
    K_ev4(i) = K4(solut(:,i));
end

figure
hold on
plot(time, K_ev1, 'r', time(mom_switch1), K_ev1(mom_switch1), 'c*', ...
    time(mom_switch2), K_ev1(mom_switch2), 'm*');
plot(time, K_ev2, 'b', time(mom_switch1), K_ev2(mom_switch1), 'c*', ...
    time(mom_switch2), K_ev2(mom_switch2), 'm*');
plot(time, K_ev3, 'g', time(mom_switch1), K_ev3(mom_switch1), 'c*', ...
    time(mom_switch2), K_ev3(mom_switch2), 'm*');
plot(time, K_ev4, 'k', time(mom_switch1), K_ev4(mom_switch1), 'c*', ...
    time(mom_switch2), K_ev4(mom_switch2), 'm*');
se = [time(1) time(end)];
plot(se, [K1_lim K1_lim], 'r--', se, [K2_lim K2_lim], 'b--', ...
    se, [K3_lim K3_lim], 'g--', se, [K4_lim K4_lim], 'k--');
xlabel('t');
ylabel('K');
legend('K_1','sw1','sw2','K_2','sw1','sw2','K_3','sw1','sw2','K_4','sw1','sw2');
grid minor
hold off

%% d2

d2_min = (c(3)*solut(2,:) + u_min(2) - r(3))/(b(3));
d2_max = (c(3)*solut(2,:) + u_max(2) - r(3))/(b(3));

%%

figure
ax1 = subplot(3,2,1);
plot(time, solut(1,:), time_switch1, x_switch1(1,:),'co',...
    time_switch2, x_switch2(1,:),'mo', ...
    [min(time) max(time)], [V1(1), V1(1)], 'r');
xlabel('t');
ylabel('x_1');
legend('x_1(t)', 'switches1','switches2','P_1');
grid minor    

ax2 = subplot(3,2,2);
plot(time, solut(2,:), time_switch1, x_switch1(2,:),'co',...
    time_switch2, x_switch2(2,:),'mo', ...
    [min(time) max(time)], [V1(2) V1(2)], 'r', ...
    [min(time) max(time)], [V2(2) V2(2)], 'r');
xlabel('t');
ylabel('x_2');
legend('x_2(t)', 'switches1','switches2'),%'P_2');
grid minor  

ax3 = subplot(3,2,3);
plot(time, solut(3,:), time_switch1, x_switch1(3,:),'co', ...
    time_switch2, x_switch2(3,:),'mo',...
    [min(time) max(time)], [V1(3) V1(3)],'r');
xlabel('t');
ylabel('x_3');
legend('x_3(t)', 'switches1','switches2','P_3');
grid minor  

ax4 = subplot(3,2,4);
% plot(time, solut(4,:), time_switch1, x_switch1(4,:),'co',...
%     time_switch2, x_switch2(4,:),'mo',...
%     [min(time) max(time)], [V1(4) V1(4)], 'r', ...
%     [min(time) max(time)], [V3(4) V3(4)], 'r', ...
%     [min(time) max(time)], [V2(4) V2(4)], 'k', ...
%     [min(time) max(time)], [V4(4) V4(4)], 'k');

plot(time, solut(4,:), time_switch1, x_switch1(4,:),'co',...
     time_switch2, x_switch2(4,:),'mo', time, d2_min, 'r', time, d2_max, 'r'); 
xlabel('t');
ylabel('x_4');
legend('x_4(t)', 'switches1','switches2')%,'P_4');
grid minor  

ax5 = subplot(3,2,5);
plot(solut(2,:), solut(4,:), solut(2,end), solut(4,end), 'g*', ...
    [V1(2) V2(2)], [V1(4) V2(4)], 'r', [V2(2) V4(2)], [V2(4) V4(4)], 'r',...
    [V4(2) V3(2)], [V4(4) V3(4)], 'r', [V3(2) V1(2)], [V3(4) V1(4)], 'r');
xlabel('x_2');
ylabel('x_4');
legend('x_4(x_2)', 'end','P_4,P_2')
grid minor

linkaxes([ax4,ax3,ax2,ax1],'x'); 

%%
figure
plot(solut(1,1:mom_switch2(1)),solut(2,1:mom_switch2(1)), 'k', [V1(1) V1(1)],...
    [V1(2) V2(2)],'r', solut(1,1), solut(2,1), 'g*');
grid minor
%%

%mom_temp = min(mom_switch1(end),mom_switch2(end));
mom_temp = 1;

figure
hold on
plot3(time(mom_temp:end), solut(2,mom_temp:end), solut(4,mom_temp:end), 'b')

[T,X2] = meshgrid(time(1):time(end),V1(2):V2(2));

P4_fun1 = (c(3)*X2 + u_min(2) - r(3))/b(3);
P4_fun2 = (c(3)*X2 + u_max(2) - r(3))/b(3);

surf(T,X2,P4_fun1, 'FaceAlpha',0.1);
surf(T,X2,P4_fun2, 'FaceAlpha',0.1)

[T,X4] = meshgrid(time(1):time(end),V1(4):V3(4));
P2_fun1 = V1(2) + 0*X4;
surf(T,P2_fun1,X4, 'FaceAlpha',0.1);

[T,X4] = meshgrid(time(1):time(end),V2(4):V4(4));
P2_fun1 = V2(2) + 0*X4;
surf(T,P2_fun1,X4, 'FaceAlpha',0.1);

xlabel('t');
ylabel('x_2');
zlabel('x_4');
hold off

%%
figure
plot([V1(2) V2(2)], [V1(4) V2(4)], 'r', [V2(2) V4(2)], [V2(4) V4(4)], 'r',...
    [V4(2) V3(2)], [V4(4) V3(4)], 'r', [V3(2) V1(2)], [V3(4) V1(4)], 'r', x_0(2),x_0(4), 'g*');
grid minor