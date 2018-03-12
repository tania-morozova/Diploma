clear
clc
%%
%nobody extincts?
r = [1,1,1,1]';
b = [1,1,1,0]';
c = [0,1,1,1]';
u_min = 1;
u_max = 2;

x_0 = [10,1,1,1]';

t_0 = 0;
t_1 = 30;

%%

ok = u_min - (r(3)*b(1) - c(3)*r(1))/c(3) > 0
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

delta_t = 10^(-4);


%%

P_curr = P(u_min);
if x_0(1) > P_curr(1)
    u_curr = u_min;
elseif x_0(1) < P_curr(1)
    u_curr = u_max;
else disp('x_0 = P_1');
end

options = odeset('Events',@(t,x)events_func(t,x,P_curr));

sol = ode45(@(t,x) f(t,x,u_curr), t_0:delta_t:t_1, x_0, options);
solut = sol.y;
time = sol.x;
mom_switch = numel(time);
time_switch = time(end);
x_switch = solut(:,end);

num = 0;
while sol.ie
    u_curr = u_max - u_curr + u_min;    
        
    sol = ode45(@(t,x) f(t,x,u_curr), time(end):delta_t:t_1, solut(:,end), options);
    
    solut = [solut, sol.y];
    time = [time, sol.x];
    mom_switch = [mom_switch, mom_switch(end) + numel(sol.x)];
    time_switch = [time_switch, time(end)];
    x_switch = [x_switch, solut(:,end)];
    num = num+1;
end

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
legend('K_1','K_2');
grid minor

%%
figure
subplot(2,2,1)
plot(time, solut(1,:), time_switch, x_switch(1,:),'go');
xlabel('t');
ylabel('x_1');
grid minor    

subplot(2,2,2)
plot(time, solut(2,:), time_switch, x_switch(2,:),'go');
xlabel('t');
ylabel('x_2');
grid minor  

subplot(2,2,3)
plot(time, solut(3,:), time_switch, x_switch(3,:),'go');
xlabel('t');
ylabel('x_3');
grid minor  

subplot(2,2,4)
plot(time, solut(4,:), time_switch, x_switch(4,:),'go');
xlabel('t');
ylabel('x_4');
grid minor  
