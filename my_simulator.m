%% Load the parameter of the platform

clc
clear
load('CRAFT_configuration.mat', 'configuration');


% Cable exit points on the pulleys are denoted as A, 
% while cable anchor points on the platform are denoted as B.
A= configuration.A;
B= configuration.B;
frame_node = configuration.frame_nodes; %CDPR frame position 3*8
mp_node = configuration.mp_nodes; % moving platform position 4


%% define the parameter here
payload_mass = 10 ; %kg
g = 9.81; % N/kg
%desired position of the platform
pos = [1.5,2.0,1.5,0,0,0]';
tau_max = 100; % maximum tension
tau_min = 1; % minimum tension to avoid the deform of the cable 


p_d = [sum(frame_node(1,:))/8 ; sum(frame_node(2,:))/8 ;sum(frame_node(3,:))/8];
p_d = pos(1:3);
% To normalize
B_norm_x = sum(B(1,:))/8;
B_norm_y = sum(B(2,:))/8;
B_norm_z = sum(B(3,:))/8;
B_norm = [B_norm_x;B_norm_y;B_norm_z];

% B_p is when mass payload reach the desired position
% B_p = B - B_norm* ones(1,8);
B_p = B + p_d* ones(1,8);


% no need to normalize frame 
% mp_node_norm_x = sum(mp_node(1,:))/8;
% mp_node_norm_y = sum(mp_node(2,:))/8;
% mp_node_norm_z = sum(mp_node(3,:))/8;
% mp_node_norm = [mp_node_norm_x;mp_node_norm_y;mp_node_norm_z];
% mp_node = mp_node - mp_node_norm* ones(1,8);
mp_node = mp_node + p_d* ones(1,8);



plot3(frame_node(1,:),frame_node(2,:),frame_node(3,:),'ok');
hold on
plot3(A(1,:),A(2,:),A(3,:),'or')
hold on
plot3(B_p(1,:),B_p(2,:),B_p(3,:),'ob')

%% Plot the frame of the CRAFT

% eight vertice of the frame cube
verts = configuration.frame_nodes';


% faces connected by nodes
faces=[1 7 4 6;
       2 5 3 8;
       2 4 6 8;
       1 3 5 7;
       1 3 8 6;
       2 5 7 4];

% plot the cube
patch('Faces',faces,'Vertices',verts,'Facecolor','none','EdgeColor','k','LineWidth',2)
axis equal
axis([-0.5 4 -0.5 5 0 3])

hold

%% Plot the cube of the mass payload

% eight vertice of the frame cube
verts = mp_node';


% faces connected by nodes
faces=[1 7 4 6;
       2 5 3 8;
       2 4 6 8;
       1 3 5 7;
       1 3 8 6;
       2 5 7 4];

% plot the cube
patch('Faces',faces,'Vertices',verts,'Facecolor','none','EdgeColor','k')

hold on


%% connect the exit points to the anchor points
X = [A(1,:); B_p(1,:)];
Y = [A(2,:); B_p(2,:)];
Z = [A(3,:); B_p(3,:)];
plot3(X,Y,Z,'r','LineWidth',1.5);


%% calculate wrinch matrix
x1 = pos(4);
x2 = pos(5);
x3 = pos(6);

b_R_p = [cos(x3)*cos(x2),-sin(x3)*cos(x1)+cos(x3)*sin(x2)*sin(x1),sin(x3)*sin(x1)+cos(x3)*sin(x2)*cos(x1);
cos(x2)*sin(x3),cos(x1)*cos(x3)+sin(x1)*sin(x2)*sin(x3),-cos(x3)*sin(x1)+cos(x1)*sin(x2)*sin(x3);-sin(x2),cos(x2)*sin(x1),cos(x1)*cos(x2)];
l = A - B_p - b_R_p * B;
l_length = [norm(l(:,1)),norm(l(:,2)),norm(l(:,3)),norm(l(:,4)),norm(l(:,5)),norm(l(:,6)),norm(l(:,7)),norm(l(:,8))];
u = l./l_length;

M1 = cross(b_R_p*B(:,1),u(:,1));
M2 = cross(b_R_p*B(:,2),u(:,2));
M3 = cross(b_R_p*B(:,3),u(:,3));
M4 = cross(b_R_p*B(:,4),u(:,4));
M5 = cross(b_R_p*B(:,5),u(:,5));
M6 = cross(b_R_p*B(:,6),u(:,6));
M7 = cross(b_R_p*B(:,7),u(:,7));
M8 = cross(b_R_p*B(:,8),u(:,8));
G = [0,0,-payload_mass * 9.81,0,0,0]';
M =[M1 M2 M3 M4 M5 M6 M7 M8];
W = [u;M];

%% tension distribution
N = null(W);
x0 = W\-G;
lb = tau_min*ones(8,1) - x0;
ub = tau_max*ones(8,1) - x0;
fun = @(x) norm(x0+ N*[x(1);x(2)]);
cond_A = [N;-N];
cond_b = [ub; -lb];
lamda = fmincon(fun,-20*ones(2,1),cond_A,cond_b);
Tau = x0+ N*lamda;
fprintf('The desired cable tension distribution is: [');
fprintf('%g, ', Tau(1:end-1));
fprintf('%g]\n', Tau(end));
lamda
