%%
%periodical
r = [2,10,10,10]';
b = [1,10,0.1,0]';
c = [0,0.1,5,1]';
u_min = [5, 1];
u_max = [6, 3];

x_0 = [10,1,1,1]';

t_0 = 0;
t_1 = 20;


%%
r = [1,7,1,1]';
b = [0.1,10,0.1,0]';
c = [0,0.1,10,0.1]';
u_min = [1, 1];
u_max = [2, 2];

x_0 = [10,1,1,1]';

t_0 = 0;
t_1 = 10;

%%
r = [1,2,1,2]';
b = [1,10,1,0]';
c = [0,1,10,2]';
u_min = [1, 1];
u_max = [2, 2];

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