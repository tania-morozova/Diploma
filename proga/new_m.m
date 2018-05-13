clear
clc
%%
r = [1,1,1,1]';
b = [1,1,1,0]';
c = [0,1,1,1]';
u_min = [1,1.5];
u_max = [2,2];

x_0 = [3,5,0.5,5]';

t_0 = 0;
t_1 = 8;
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

delta_t = 10^(-7);
%%
st = linspace(-1, -8, 20);
time_stop = zeros(1,numel(st));
for i = 1:numel(st)
    eps = 10^st(i);
    i
    f_s = @(t,x)f_synth(t,x,u_min,u_max,f,P,r,b,c,eps);
    options = odeset('Events',@(t,x)events_func3(t,x,c,b,r,P,u_min,u_max,eps),...
        'MaxStep', 1e-4,'RelTol',1e-13);%,'AbsTol',1e-6);
    sol = ode45(@(t,x) f_s(t,x), t_0:delta_t:t_1, x_0, options);
    time_stop(i) = sol.x(end);
end

%%
figure
semilogx(st,time_stop,'r.-','LineWidth',2);
xlabel('\alpha = lg(\epsilon)','FontSize', 16, 'FontWeight','bold');
ylabel('?????','FontSize', 16, 'FontWeight','bold');
legend({'?????'},'FontSize',14)
grid on
%%
set(gcf,'PaperUnits','inches','Paperposition',[0 0 10 8],'PaperPositionMode','manual');%14x9 inches window
print('time_eps.eps','-depsc');
%%
eps = 10^(-30);
f_s = @(t,x)f_synth(t,x,u_min,u_max,f,P,r,b,c,eps);
options = odeset('Events',@(t,x)events_func3(t,x,c,b,r,P,u_min,u_max,eps), ...
    'MaxStep', 1e-4,'RelTol',1e-13);%,'AbsTol',1e-13);
sol = ode45(@(t,x) f_s(t,x), t_0:delta_t:t_1, x_0, options);
time_stop(i) = sol.x(end);

%% d2
solut = sol.y;
time = sol.x;
d2_min = (c(3)*solut(2,:) + u_min(2) - r(3))/(b(3));
d2_max = (c(3)*solut(2,:) + u_max(2) - r(3))/(b(3));

%% x1 x2 x3 x4
figure

ax1 = subplot(3,2,1);
plot(time, solut(1,:),'b', [min(time) max(time)], [V1(1), V1(1)], 'r--')%, ...
    %time_switch1, x_switch1(1,:),'co',time_switch2, x_switch2(1,:),'mo');
xlabel('t');
ylabel('x_1');
legend('x_1(t)','P_1')%, 'switches1','switches2','P_1');
grid minor    

ax2 = subplot(3,2,2);
plot(time, solut(2,:), 'b', ...
    [min(time) max(time)], [V1(2) V1(2)], 'r--', ...
    [min(time) max(time)], [V2(2) V2(2)], 'm--')%, ...
    %time_switch1, x_switch1(2,:),'co', time_switch2, x_switch2(2,:),'mo');
xlabel('t');
ylabel('x_2');
legend('x_2(t)','P_2^{min}', 'P_2^{max}')%, 'switches1','switches2'),%'P_2');
grid minor  


ax3 = subplot(3,2,3);
plot(time, solut(3,:), 'b', [min(time) max(time)], [V1(3) V1(3)],'r--')%, ...
    %time_switch1, x_switch1(3,:),'co', time_switch2, x_switch2(3,:),'mo');
xlabel('t');
ylabel('x_3');
legend('x_3(t)', 'P_3')%, 'switches1','switches2','P_3');
grid minor  

ax4 = subplot(3,2,4);
% plot(time, solut(4,:), time_switch1, x_switch1(4,:),'co',...
%     time_switch2, x_switch2(4,:),'mo',...
%     [min(time) max(time)], [V1(4) V1(4)], 'r', ...
%     [min(time) max(time)], [V3(4) V3(4)], 'r', ...
%     [min(time) max(time)], [V2(4) V2(4)], 'k', ...
%     [min(time) max(time)], [V4(4) V4(4)], 'k');

plot(time, solut(4,:), 'b', time, d2_min, 'r--', time, d2_max, 'm--')%,...
    %time_switch1, x_switch1(4,:),'co',time_switch2, x_switch2(4,:),'mo',); 
xlabel('t');
ylabel('x_4');
legend('x_4(t)','d_{min}(x_2)','d_{max}(x_2)')%, 'switches1','switches2')%,'P_4');
grid minor  

ax5 = subplot(3,2,5);
plot(solut(2,:), solut(4,:), 'b', solut(2,end), solut(4,end), 'g*', ...
    [V1(2) V2(2)], [V1(4) V2(4)], 'r', [V2(2) V4(2)], [V2(4) V4(4)], 'r',...
    [V4(2) V3(2)], [V4(4) V3(4)], 'r', [V3(2) V1(2)], [V3(4) V1(4)], 'r');
xlabel('x_2');
ylabel('x_4');
legend('x_4(x_2)', 'end','P_4,P_2')
grid minor

linkaxes([ax4,ax3,ax2,ax1],'x'); 