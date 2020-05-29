% Pendulum with Disk, take 2

blankSlate

%%%% Define variables %%%

% positions:
syms x0_1 y0_1 o0_1 x0_2 y0_2 o0_2 real

q = [x0_1 y0_1 o0_1 x0_2 y0_2 o0_2]'; % generalized coordinates

% velocities and time:
syms dx0_1 dy0_1 do0_1 dx0_2 dy0_2 do0_2 t real

dq = [dx0_1 dy0_1 do0_1 dx0_2 dy0_2 do0_2]'; % generalized velocities

% Inertial parameters and link dimensions
syms m1 m2 I1 I2 l1 lc1 real

M = diag([m1 m1 I1 m2 m2 I2]); % Mass matrix

% External forces and moments
syms F1x F1y M1 F2x F2y M2 g real;

Qe = [F1x, F1y - m1*g, M1-M2, F2x, F2y - m2*g, M2]';

%%%% Define constraint vector

C = sym('c',[4, 1]);

rotmat = @(x) [cos(x), - sin(x); sin(x), cos(x)]; % rotation matrix from link frame to lab frame

r0_1 = q(1:2); % link 1 center of mass position in lab frame
r0_2 = q(4:5); % link 2 center of mass position in lab frame

R0_1 = rotmat(o0_1); % Rotation matrix from link 1 frame to lab frame;
R0_2 = rotmat(o0_2); % Rotation matrix from link 2 frame to lab frame;

r1_0 = [-lc1; 0]; % Position of origin (pin joint on ground) in link 1 frame
r1_A = [l1-lc1; 0]; % Position of pin joint A in link 1 frame
r2_A = [0; 0]; % Position of pin joint A in link 2 frame

C(1:2) = r0_1 + R0_1*r1_0; % Pin joint at [0;0];
C(3:4) = r0_1 + R0_1*r1_A - (r0_2 - R0_2*r2_A); % Pin joint at A

% Define Constraint jacobian

Cq = jacobian(C,q);
Cqt = diff(Cq,t); 
Ct = diff(C,t);

% Compute Quadratic Velocity Vector (?)
Qd = -( jacobian(Cq*dq,q)*dq + 2*Cqt*dq + Ct );

% Compute augmented matrix

n = size(Cq,1);
A = [M Cq';Cq zeros(n)];

% Define RHS vector:

Q = [Qe;Qd];

%% Set parameters

mset = [1 1];
Iset = [0.5 1];
lset = 2;
lcset = 1;
gset = 9.8;

oi_set = [-pi/2 0];
doi_set = [0 0];

Torque_fun = @(t) [0,sin(t)*sum(mset)*gset];
Force_fun = @(t) [0,0,0,0];

%% Set numerical LHS, RHS functions

Aset = subs(A,[m1 m2 I1 I2 l1 lc1],[mset Iset lset lcset]);
Qset = subs(Q,[g,m1,m2,l1,lc1],[gset, mset, lset lcset]);

matlabFunction(Aset,Qset,'vars',{[q;dq]',[M1,M2],[F1x,F1y,F2x,F2y]},'file','AQ_PendDisk_num','sparse',true);

%% Solve C to get initial conditions

Cset = subs(C,[o0_1 o0_2 l1 lc1],[oi_set lset lcset]);
q_initial = solve(Cset == 0);

Cq_initial_set = subs(Cq,[o0_1 o0_2 l1 lc1],[oi_set lset lcset]);
Cq_initial_set = double(Cq_initial_set);

Cqdq_initial_set = subs(Cq_initial_set*dq,[do0_1 do0_2],doi_set); % first derivative, as there is no explicit dependence on time (otherwise would have to add diff(C,t))

dq_initial = solve( Cqdq_initial_set == 0);

% Compile initial conditions
q_initial = [q_initial.x0_1 q_initial.y0_1 oi_set(1) ...
    q_initial.x0_2 q_initial.y0_2 oi_set(2)];
q_initial = double(q_initial);

dq_initial = double([dq_initial.dx0_1 dq_initial.dy0_1 doi_set(1) ...
    dq_initial.dx0_2 dq_initial.dy0_2 doi_set(2)]);

%% Simulate over time
t0 = 0;
tf = 2*pi;
n = 2;
tspan = linspace(t0,tf,n);

opts = odeset('RelTol',1e-10,'abstol',1e-10,'normcontrol','on','refine',100);
[t,qdq] = ode23(@(t,z) PendDiskODEfun(t,z,Torque_fun(t),Force_fun(t)),tspan,[q_initial,dq_initial]); 

q_num = qdq(:,1:6);
dq_num = qdq(:,7:12);

close all
plot(t,qdq(:,[3,6])*180/pi)
legend('\theta_1','\theta_2')

%% Get reaction forces

n = length(t);
dq_ddq_L = zeros(n,2*6+4);

for i = 1:length(t)
   dq_ddq_L(i,:) = PendDiskODEfun(t(i),qdq(i,:),Torque_fun(t(i)),Force_fun(t(i)),true); 
end

% Define the jacobian of C as a function

Cq_num = subs(Cq,[l1 lc1],[lset lcset]);
Cq_fun = matlabFunction(Cq_num,'vars',{q'},'sparse',true);

ddq_num = dq_ddq_L(:,7:12);
Lambda = dq_ddq_L(:,13:16);

Qc_num = NaN(n,6);
for i = 1:n
    Qc_num(i,:) = -Cq_fun(q_num(i,:))'*Lambda(i,:)';
end
% Qc_num gives resultant constraint forces and moments.
% To resolve components due to individual constraints, we need the
% component of the constraint jacobian due to those constraints
% e.g. ground contact is C(1:2)

% Cq_ground = jacobian(C(1:2),q(1:3)); % this will be equal to Cq(1:2,1:3);
% disp(Cq(1:2,1:3) - Cq_ground); % = [0,0,0;0,0,0];

% Since we already have the Jacobian of the constraint vector due to ground
% contact, we need to identify the associated ground contact lagrange
% multipliers (which should be Lambda(1:2));
% The resultants of the ground contact reaction forces are then
%   Qc_ground = - Cq_ground'*Lambda_ground
%             = - Cq(1:2,1:3)'*Lambda(1:2)

Qc_ground_res = NaN(n,3);

for i = 1:n
    Cq_ground = Cq_fun(q_num(i,:));
    Cq_ground = Cq_ground(1:2,1:3);
    Qc_ground_res(i,1:3) = (- Cq_ground'*Lambda(i,1:2)')'; % Note that this is NOT equal to Qc_num(i,1:3)
end

% Now, we identify an equipollent system at ground contact. The only thing
% that may change is the moment
F_ground_res = [Qc_ground_res(:,1:2),zeros(n,1)]; % 3D forces due to ground contact

r0_1to0 = zeros(n,3);
for i = 1:n
    r0_1to0(i,1:2) = double(subs(R0_1*r1_0,[o0_1,lc1],[q_num(i,3),lcset]));
end
Torque0_1react = -cross(r0_1to0',F_ground_res')';

Qc_ground = double(Qc_ground_res + Torque0_1react);

disp(['Maximum torque at origin = ',num2str(max(Qc_ground(:,3)),2)]) % should be ~ zero


%% Animation

FR = 30; % Frames per second; always played back at 30 FPS
filename = 'PendDisk';
saveanimation = true;
tQ = linspace(t0,tf,round(FR*diff([t0 tf])+1));

qQ = interp1(t,q_num(:,1:6),tQ);

Qc_ground_Q = interp1(t,Qc_ground,tQ)/(sum(mset)*gset);

cs = @(a) [cos(a),sin(a)];

O_1 = qQ(:,1:2) - lcset*cs(qQ(:,3));
A_1 = qQ(:,1:2) + lcset*cs(qQ(:,3));
A_2 = qQ(:,4:5);

rG = sqrt(Iset(2)/mset(2)); % radius of gyration of the disk

disk_line = A_2 + rG*cs(qQ(:,6));

lmax = lset + rG;

close all
figure
writerObj = VideoWriter(filename,'MPEG-4');
writerObj.FrameRate = 30;
if saveanimation
    open(writerObj);
end
for i = 1:length(tQ)
    cla
    % plot static points
    plot(0,0,'ko','markersize',10)
    hold on
    
    % plot link 1
    plot([O_1(i,1) A_1(i,1)],[O_1(i,2) A_1(i,2)],'rs-','linewidth',2)
    
    % plot Disk
    [x_circ,y_circ] = circle_xy(A_2(i,1),A_2(i,2),rG);
    plot(x_circ,y_circ,'b-','linewidth',2)
    
    % plot disk reference line
    plot([A_2(i,1),disk_line(i,1)],[A_2(i,2),disk_line(i,2)],'b-','linewidth',1)
    
    % plot COMs
    plot(qQ(i,[1,4]),qQ(i,[2,5]),'kx')
    
    % plot reaction forces
    
    quiver(0,0,Qc_ground_Q(i,1),Qc_ground_Q(i,2))
    
    text(0,1,['t = ',num2str(tQ(i),'%.2f')],'units','normalized','verticalalignment','top')
    xlim(1.5*lmax*[-1 1])
    ylim(1.5*lmax*[-1 1])
    drawnow
    if saveanimation
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
        
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 1
            imwrite(imind,cm,[filename,'.gif'],'gif', 'Loopcount',inf,'DelayTime',1/FR);
        else
            imwrite(imind,cm,[filename,'.gif'],'gif','WriteMode','append','DelayTime',1/FR);
        end
    else
        pause(1/28)
    end
end
close(writerObj)

