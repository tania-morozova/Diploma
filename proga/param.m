%%
%periodical
r = [2,10,10,10]';
b = [1,10,0.1,0]';
c = [0,0.1,1,1]';
u_min = 10;
u_max = 20;

x_0 = [10,1,1,1]';

t_0 = 0;
t_1 = 10;

%% Ex1
r = [1,1,1,1]';
b = [1,1,1,0]';
c = [0,1,1,1]';
u_min = 1;
u_max = 2;

x_0 = [2,2.39,1,1.5]';

t_0 = 0;
t_1 = 15;

%% Ex2
%x3 not separated from P3 
r = [1,1,1,1]';
b = [1,1,1,0]';
c = [0,1,1,1]';
u_min = 1;
u_max = 2;

x_0 = [2,2.29,0.9,1]';

t_0 = 0;
t_1 = 50;

%%
%strange inf
r = [1,7,1,1]';
b = [0.1,10,0.1,0]';
c = [0,0.1,10,0.1]';
u_min = 1;
u_max = 2;

x_0 = [10,1,1,1]';

t_0 = 0;
t_1 = 10;

%%
%nobody extincts
r = [1,1,1,1]';
b = [1,1,1,0]';
c = [0,1,1,1]';
u_min = 1;
u_max = 2;

x_0 = [3,1,1,1]';

t_0 = 0;
t_1 = 100;

%%
r = [1,1,1,1]';
b = [1,1,1,0]';
c = [0,1,1,1]';
u_min = 1;
u_max = 2;

x_0 = [3,10,10,10]';

t_0 = 0;
t_1 = 10;