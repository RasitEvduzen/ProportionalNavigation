clc; clear all; close all;
% True Proportional Navigation - 3D TPN
% Guidance: a = N * Vc * (Omega x r_hat)
% Written By: Rasit Evduzen
% 28-Apr-2026

%% Parameters
N       = 4;
Ts      = 1e-2;
Tmax    = 15;
Rmin    = 10;
Amax    = 30 * 9.81;
StepVis = 50;

%% Initial Conditions (3D)
Rm0 = [0;    0;   0];    % Missile initial position [m]
Vm0 = [800;  0;   0];    % Missile initial velocity [m/s], speed = 800 m/s along +X
Rt0 = [6000; 1500; 500]; % Shahed initial position  [m]
Vt0 = [-180; 0;  50];    % Shahed initial velocity  [m/s], head-on + slight ascent

%% Load STL Geometry
MissileMesh = stlread('missile.STL');
ShahedMesh  = stlread('shahed.stl');

%% Simulate
[tm, tt, t_h, Omega_h, Vc_h, R_h, acmd_h, miss, hit] = ...
    Sim3DTPN(Rm0, Vm0, Rt0, Vt0, N, Ts, Tmax, Rmin, Amax);

fprintf('Miss Distance: %.3f m  |  Hit: %d\n', miss, hit);

%% Animate
NoD = length(t_h);
figure('units','normalized','position',[0 0 1 1],'color','w')
for i = 1:NoD
    if mod(i, StepVis) == 0 || i == NoD
        clf
        DrawFrame(tm, tt, t_h, Omega_h, Vc_h, R_h, acmd_h, i, N, MissileMesh, ShahedMesh)
        drawnow
    end
end

%% Utility Functions

function [traj_m, traj_t, t_h, Omega_h, Vc_h, R_h, acmd_h, miss, hit] = ...
    Sim3DTPN(Rm0, Vm0, Rt0, Vt0, N, Ts, Tmax, Rmin, Amax)

Rm = Rm0;
Vm = Vm0;
Rt = Rt0;
Vt = Vt0;
t   = 0;
hit = 0;

Nmax    = round(Tmax/Ts) + 2;
traj_m  = zeros(Nmax, 3);
traj_t  = zeros(Nmax, 3);
t_h     = zeros(Nmax, 1);
Omega_h = zeros(Nmax, 1);
Vc_h    = zeros(Nmax, 1);
R_h     = zeros(Nmax, 1);
acmd_h  = zeros(Nmax, 1);

idx = 1;
traj_m(1,:) = Rm';
traj_t(1,:) = Rt';
R_h(1)      = norm(Rt - Rm);

while t < Tmax

    % LOS vector and range
    r     = Rt - Rm;               % LOS vector: target pos - missile pos
    R     = norm(r);               % scalar range
    r_hat = r / R;                 % unit LOS vector

    % Relative velocity
    V_rel = Vt - Vm;               % relative velocity vector

    % Closing velocity: positive when range is decreasing
    Vc = -(V_rel' * r_hat);

    % Termination
    if R < Rmin
        hit = 1;
        break;
    end
    if Vc < 0 && t > 0.5
        break;
    end

    % LOS angular rate vector:  Omega = (r x V_rel) / R^2
    Omega = cross(r, V_rel) / R^2;

    % TPN guidance law: a = N * Vc * (Omega x r_hat)
    % acceleration is perpendicular to LOS (not missile velocity)
    % theoretically optimal for constant velocity targets
    a_cmd = N * Vc * cross(Omega, r_hat);

    % Saturate to kinematic limit (direction preserved)
    a_norm = norm(a_cmd);
    if a_norm > Amax
        a_cmd = a_cmd * (Amax / a_norm);
    end

    % RK4 integration
    Sm = RK4([Rm; Vm], a_cmd,   Ts);
    St = RK4([Rt; Vt], [0;0;0], Ts);
    Rm = Sm(1:3);
    Vm = Sm(4:6);
    Rt = St(1:3);
    Vt = St(4:6);
    t  = t + Ts;

    idx = idx + 1;
    traj_m(idx,:) = Rm';
    traj_t(idx,:) = Rt';
    t_h(idx)      = t;
    Omega_h(idx)  = norm(Omega) * 180/pi;
    Vc_h(idx)     = Vc;
    R_h(idx)      = R;
    acmd_h(idx)   = a_norm / 9.81;
end

traj_m  = traj_m(1:idx,:);
traj_t  = traj_t(1:idx,:);
t_h     = t_h(1:idx);
Omega_h = Omega_h(1:idx);
Vc_h    = Vc_h(1:idx);
R_h     = R_h(1:idx);
acmd_h  = acmd_h(1:idx);
miss    = norm(Rt - Rm);
end

function x_next = RK4(x, a, Ts)
k1 = Dyn3D(x,            a);
k2 = Dyn3D(x + Ts/2*k1,  a);
k3 = Dyn3D(x + Ts/2*k2,  a);
k4 = Dyn3D(x + Ts*k3,    a);
x_next = x + (Ts/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function dx = Dyn3D(x, a)
dx = [x(4); x(5); x(6); a(1); a(2); a(3)];
end

function DrawFrame(tm, tt, t_h, Omega_h, Vc_h, R_h, acmd_h, k, N, meshM, meshT)
mx = tm(k,1); my = tm(k,2); mz = tm(k,3);
tx = tt(k,1); ty = tt(k,2); tz = tt(k,3);

Vmd = [1 0 0];
Vtd = [-1 0 0];
if k > 1
    dM = tm(k,:) - tm(k-1,:);
    dT = tt(k,:) - tt(k-1,:);
    if norm(dM) > 0, Vmd = dM / norm(dM); end
    if norm(dT) > 0, Vtd = dT / norm(dT); end
end

subplot(4,3,[1 2 4 5 7 8 10 11]), hold on, grid on, box on, view(35,20)
plot3(tm(1:k,1), tm(1:k,2), tm(1:k,3), 'b-', 'LineWidth', 2)
plot3(tm(1,1), tm(1,2), tm(1,3), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b')
PlotBody([mx;my;mz], Vmd, meshM, [0.15 0.40 0.85],  90,   0)
PlotBody([tx;ty;tz], Vtd, meshT, [0.85 0.30 0.10],   0, 180)
axis equal
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]')
title(sprintf('TPN 3D  (N=%d)  |  R=%.0fm  |  t=%.2fs', N, R_h(k), t_h(k)), 'FontSize', 11)
legend({'Missile'}, 'Location', 'northeast', 'FontSize', 9)

subplot(4,3,3), plot(t_h(1:k), R_h(1:k), 'k-', 'LineWidth', 2), grid on, hold on
xline(t_h(k), 'b--', 'LineWidth', 1)
xlabel('$t$ [s]', 'Interpreter', 'latex')
ylabel('$R$ [m]', 'Interpreter', 'latex')
title('Range $R(t)$', 'Interpreter', 'latex')
xlim([0 t_h(end)]); ylim([0 max(R_h)*1.05])

subplot(4,3,6), plot(t_h(1:k), Vc_h(1:k), 'Color', [0.35 0.7 0.15], 'LineWidth', 2), grid on, hold on
xline(t_h(k), 'b--', 'LineWidth', 1)
xlabel('$t$ [s]', 'Interpreter', 'latex')
ylabel('$V_c$ [m/s]', 'Interpreter', 'latex')
title('Closing Velocity $V_c(t)$', 'Interpreter', 'latex')
xlim([0 t_h(end)])

subplot(4,3,9), plot(t_h(1:k), Omega_h(1:k), 'Color', [0.9 0.4 0], 'LineWidth', 2), grid on, hold on
yline(0, 'k--', 'LineWidth', 1.5)
xline(t_h(k), 'b--', 'LineWidth', 1)
xlabel('$t$ [s]', 'Interpreter', 'latex')
ylabel('$|\Omega|$ [deg/s]', 'Interpreter', 'latex')
title('$|\Omega|(t) \rightarrow 0$', 'Interpreter', 'latex')
xlim([0 t_h(end)])

subplot(4,3,12)
area(t_h(1:k), acmd_h(1:k), 'FaceColor', [0.1 0.7 0.35], 'FaceAlpha', 0.3, ...
     'EdgeColor', [0.1 0.7 0.35], 'LineWidth', 1.5)
grid on; hold on
yline(30, 'r--', 'LineWidth', 1.2)
yline( 0, 'k-',  'LineWidth', 0.8)
xline(t_h(k), 'b--', 'LineWidth', 1)
xlabel('$t$ [s]', 'Interpreter', 'latex')
ylabel('$|a_c|$ [g]', 'Interpreter', 'latex')
title('$\mathbf{a} = N \cdot V_c \cdot (\Omega \times \hat{r})$', 'Interpreter', 'latex')
xlim([0 t_h(end)])

sgtitle(sprintf('True Proportional Navigation  $|$  $N = %d$', N), ...
        'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold')
end

function PlotBody(pos, vel_dir, mesh, facecolor, preRotY, preRotX)
yaw  = atan2d(vel_dir(2), vel_dir(1));
elev = atan2d(vel_dir(3), sqrt(vel_dir(1)^2 + vel_dir(2)^2));
R    = rotz(yaw) * roty(-elev) * roty(preRotY) * rotx(preRotX);
T    = [R, pos; 0 0 0 1];
pts  = TransformGeo(mesh.Points, T);
trisurf(mesh.ConnectivityList, pts(:,1), pts(:,2), pts(:,3), ...
        'FaceColor', facecolor, 'EdgeColor', 'none')
camlight; lighting gouraud
end

function transPoints = TransformGeo(points, Tmatrix)
[sze, ~] = size(points);
points(:, 4) = 1;
transPoints  = zeros(sze, 4);
for i = 1:sze
    transPoints(i,:) = Tmatrix * points(i,:)';
end
end