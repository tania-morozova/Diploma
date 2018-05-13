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

K = @(x,u) KK(x,u,P,b,c);
K1 = @(x) K(x,u_min);
K2 = @(x) K(x,[u_max(1),u_min(2)]);
K3 = @(x) K(x,[u_min(1),u_max(2)]);
K4 = @(x) K(x,u_max);

delta_t = 10^(-5);

%%
eps = 10^(-20);
P_curr = P(u_min);
f_s = @(t,x)f_synth(t,x,u_min,u_max,f,P,r,b,c,eps);
f_s(0,x_0);
%%

options = odeset('Events',@(t,x)events_func(t,x,P_curr), 'MaxStep', 1e-3);

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

%% d2

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


%% x1 x2 x3 x4 eps
figure

ax1 = subplot(2,2,1);
plot(time, solut(1,:),'b.-', [min(time) max(time)], [V1(1), V1(1)], 'r', ...
    [min(time) max(time)], [V1(1)-eps, V1(1)-eps], 'g--',...
    [min(time) max(time)], [V1(1)+eps, V1(1)+eps], 'g--')%, ...
    %time_switch1, x_switch1(1,:),'co',time_switch2, x_switch2(1,:),'mo');
xlabel('t');
ylabel('x_1');
legend('x_1(t)','P_1','P_1 - \epsilon','P_1 + \epsilon')%, 'switches1','switches2','P_1');
grid minor    

ax2 = subplot(2,2,2);
plot(time, solut(2,:), 'b.-', ...
    [min(time) max(time)], [V1(2) V1(2)], 'r', ...
    [min(time) max(time)], [V2(2) V2(2)], 'm', ...
    [min(time) max(time)], [V1(2)-eps V1(2)-eps], 'g--', ...
    [min(time) max(time)], [V2(2)+eps V2(2)+eps], 'g--')%, ...
    %time_switch1, x_switch1(2,:),'co', time_switch2, x_switch2(2,:),'mo');
xlabel('t');
ylabel('x_2');
legend('x_2(t)','P_2^{min}', 'P_2^{max}','P_2^{min}-\alpha_2\epsilon', 'P_2^{max}+\alpha_2\epsilon')%, 'switches1','switches2'),%'P_2');
grid minor  


ax3 = subplot(2,2,3);
plot(time, solut(3,:), 'b.-', [min(time) max(time)], [V1(3) V1(3)],'r', ...
    [min(time) max(time)], [V1(3)-eps, V1(3)-eps], 'g--',...
    [min(time) max(time)], [V1(3)+eps, V1(3)+eps], 'g--')%, ...
    %time_switch1, x_switch1(3,:),'co', time_switch2, x_switch2(3,:),'mo');
xlabel('t');
ylabel('x_3');
legend('x_3(t)', 'P_3','P_3 - \epsilon','P_3 + \epsilon','Location','northwest')%, 'switches1','switches2','P_3');
grid minor  

ax4 = subplot(2,2,4);
% plot(time, solut(4,:), time_switch1, x_switch1(4,:),'co',...
%     time_switch2, x_switch2(4,:),'mo',...
%     [min(time) max(time)], [V1(4) V1(4)], 'r', ...
%     [min(time) max(time)], [V3(4) V3(4)], 'r', ...
%     [min(time) max(time)], [V2(4) V2(4)], 'k', ...
%     [min(time) max(time)], [V4(4) V4(4)], 'k');

plot(time, solut(4,:), 'b.-', time, d2_min, 'r', time, d2_max, 'm', ...
    time, d2_min-eps, 'g--', time, d2_max+eps, 'g--')%,...
    %time_switch1, x_switch1(4,:),'co',time_switch2, x_switch2(4,:),'mo',); 
xlabel('t');
ylabel('x_4');
legend('x_4(t)','d_{min}(x_2)','d_{max}(x_2)','d_{min}(x_2)-\alpha_4\epsilon','d_{max}(x_2)+\alpha_4\epsilon', 'Location','southwest')%, 'switches1','switches2')%,'P_4');
grid minor  
%% x_1(t)
figure
hold on
plot(time, solut(1,:), time_switch1, x_switch1(1,:),'co',...
    time_switch2, x_switch2(1,:),'mo', ...
    [min(time) max(time)], [V1(1), V1(1)], 'r');
xlabel('t');
ylabel('x_1');
%legend('x_1(t)', 'switches1','switches2','P_1');
grid minor    
%% (x_2 x_4 t)
%mom_temp = min(mom_switch1(end),mom_switch2(end));
mom_temp = 1;
t_c = 2;

figure1 = figure
axes1 = axes('Parent',figure1);
hold on
plot3(time(mom_temp:end), solut(2,mom_temp:end), solut(4,mom_temp:end), 'b', 'LineWidth',2)

[T,X2] = meshgrid(time(mom_temp)-t_c:t_1,V1(2):V2(2));

P4_fun1 = (c(3)*X2 + u_min(2) - r(3))/b(3);
P4_fun2 = (c(3)*X2 + u_max(2) - r(3))/b(3);

surf(T,X2,P4_fun1, 'FaceAlpha',0.5);
surf(T,X2,P4_fun2, 'FaceAlpha',0.5)

[T,X4] = meshgrid(time(mom_temp)-t_c:t_1,V1(4):V3(4));
P2_fun1 = V1(2) + 0*X4;
surf(T,P2_fun1,X4, 'FaceAlpha',0.5);

[T,X4] = meshgrid(time(mom_temp)-t_c:t_1,V2(4):V4(4));
P2_fun1 = V2(2) + 0*X4;
surf(T,P2_fun1,X4, 'FaceAlpha',0.5);

xlabel('t');
ylabel('x_2');
zlabel('x_4');
view(axes1,[-24 14.8]);
hold off

%% K_i(x,u)
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
plot(time, K_ev1, 'r')%, time(mom_switch1))%, K_ev1(mom_switch1), 'c*', ...
    %time(mom_switch2), K_ev1(mom_switch2), 'm*');
plot(time, K_ev2, 'b')%, time(mom_switch1))%, K_ev2(mom_switch1), 'c*', ...
    %time(mom_switch2), K_ev2(mom_switch2), 'm*');
plot(time, K_ev3, 'g')%, time(mom_switch1))%, K_ev3(mom_switch1), 'c*', ...
    %time(mom_switch2), K_ev3(mom_switch2), 'm*');
plot(time, K_ev4, 'k')%, time(mom_switch1))%, K_ev4(mom_switch1), 'c*', ...
    %time(mom_switch2), K_ev4(mom_switch2), 'm*');
se = [time(1) time(end)];
plot(se, [K1_lim K1_lim], 'r--', se, [K2_lim K2_lim], 'b--', ...
    se, [K3_lim K3_lim], 'g--', se, [K4_lim K4_lim], 'k--');
xlabel('t');
ylabel('K');
legend('K_1','K_2','K_3','K_4','K_1^{min}','K_2^{min}','K_3^{min}','K_4^{min}');
grid minor
hold off


%% P2 P4
ax_l = 0.2
figure
hold on
plot([V1(2) V2(2)], [V1(4) V2(4)], 'r', 'LineWidth', 2)
plot([V2(2) V4(2)], [V2(4) V4(4)], 'r', 'LineWidth', 2)
plot([V4(2) V3(2)], [V4(4) V3(4)], 'r', 'LineWidth', 2)
plot([V3(2) V1(2)], [V3(4) V1(4)], 'r', 'LineWidth', 2);
axis([V2(1) - ax_l, V2(2) + ax_l,  V1(4) - ax_l, V4(4) + ax_l]) 
grid minor
xlabel('x_2','FontSize', 18, 'FontWeight','bold')
ylabel('x_4','FontSize', 18, 'FontWeight','bold')

%%
set(gcf,'PaperUnits','inches','Paperposition',[0 0 10 8],'PaperPositionMode','manual');%14x9 inches window
print('ord3.eps','-depsc');
%% K(x,u)

f_ut = @(t,x) f_u(t,x,u_min,u_max,P,r,b,c);

K_gen = zeros(1,numel(time));
for i = 1:numel(time)
    K_gen(i) = K(solut(:,i),f_ut_pr(f_ut(time(i),solut(:,i)),u_min,u_max));
end

figure
plot(time, K_gen, 'r');
grid minor

%%
l_p = 1.5;
figure

ax1 = subplot(2,2,1);
hold on
plot(time, solut(1,:),'b','LineWidth',2) 
plot([min(time) max(time)], [V1(1), V1(1)], 'r--','LineWidth',l_p)%, ...
    %time_switch1, x_switch1(1,:),'co',time_switch2, x_switch2(1,:),'mo');
xlabel('t','FontSize', 18, 'FontWeight','bold');
ylabel('x_1','FontSize', 18, 'FontWeight','bold');
legend({'x_1(t)','P_1'},'FontSize', 12,'Location','southeast')%, 'switches1','switches2','P_1');
grid minor    
hold off

ax2 = subplot(2,2,2);
hold on
plot(time, solut(2,:), 'b','LineWidth',2)
plot([min(time) max(time)], [V1(2) V1(2)], 'r--', 'LineWidth',l_p)
plot([min(time) max(time)], [V2(2) V2(2)], 'g--','LineWidth',l_p)%, ...
    %time_switch1, x_switch1(2,:),'co', time_switch2, x_switch2(2,:),'mo');
xlabel('t','FontSize', 18, 'FontWeight','bold');
ylabel('x_2','FontSize', 18, 'FontWeight','bold');
legend({'x_2(t)','P_2^{min}', 'P_2^{max}'},'FontSize', 12)%, 'switches1','switches2'),%'P_2');
grid minor  
hold off

ax3 = subplot(2,2,3);
hold on
plot(time, solut(3,:), 'b', 'LineWidth',2)
plot([min(time) max(time)], [V1(3) V1(3)],'r--','LineWidth',l_p)%, ...
    %time_switch1, x_switch1(3,:),'co', time_switch2, x_switch2(3,:),'mo');
xlabel('t','FontSize', 18, 'FontWeight','bold');
ylabel('x_3','FontSize', 18, 'FontWeight','bold');
legend({'x_3(t)', 'P_3'},'FontSize', 12)%, 'switches1','switches2','P_3');
grid minor  
hold off


ax4 = subplot(2,2,4);
hold on
% plot(time, solut(4,:), time_switch1, x_switch1(4,:),'co',...
%     time_switch2, x_switch2(4,:),'mo',...
%     [min(time) max(time)], [V1(4) V1(4)], 'r', ...
%     [min(time) max(time)], [V3(4) V3(4)], 'r', ...
%     [min(time) max(time)], [V2(4) V2(4)], 'k', ...
%     [min(time) max(time)], [V4(4) V4(4)], 'k');

plot(time, solut(4,:), 'b', 'LineWidth',2)
plot(time, d2_min, 'r--','LineWidth',l_p)
plot(time, d2_max, 'g--','LineWidth',l_p)%,...
    %time_switch1, x_switch1(4,:),'co',time_switch2, x_switch2(4,:),'mo',); 
xlabel('t','FontSize', 18, 'FontWeight','bold');
ylabel('x_4','FontSize', 18, 'FontWeight','bold');
legend({'x_4(t)','P_4^{min}(x_2)','P_4^{max}(x_2)'},'FontSize', 12)%, 'switches1','switches2')%,'P_4');
grid minor 
hold off