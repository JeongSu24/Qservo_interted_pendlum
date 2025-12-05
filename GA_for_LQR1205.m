%% GA-based LQR Optimization for Rotary Inverted Pendulum (Balance Only)
%  - 목표: 손으로 inverted 근처로 올려놨을 때 (몇 도 오차) LQR이 잘 버티도록 하는 K 찾기
%  - 실제 하드웨어 상황(전압 제한, 속도 추정기 LPF)을 GA 시뮬레이션에 반영

clear; clc; close all;

fprintf('================================================\n');
fprintf('   GA-LQR Balance Optimization Start\n');
fprintf('================================================\n');

%% 1. GA 설정 ---------------------------------------------------------------
% 유전자: [q_theta, q_alpha, q_d_theta, q_d_alpha, r_val]
nVars = 5;

lb = [1,    10,   0,   0,   0.1];
ub = [500, 1000,  5,   5,   10];

options = gaoptimset( ...
    'PopulationSize', 100, ...   % 필요하면 줄여도 됨
    'Generations',    50, ...
    'Display',        'iter', ...
    'PlotFcns',       @gaplotbestf ...
);

%% 2. GA 실행 ---------------------------------------------------------------
fprintf('최적화를 시작합니다...\n');
[x_best, fval] = ga(@lqr_cost_standup, nVars, [], [], [], [], lb, ub, [], options);

%% 3. 최적 값 출력 ----------------------------------------------------------
fprintf('\n================================================\n');
fprintf('   Optimization Completed!\n');
fprintf('================================================\n');

fprintf('Best Q and R found:\n');
fprintf('   q_theta     : %.4f\n', x_best(1));
fprintf('   q_alpha     : %.4f\n', x_best(2));
fprintf('   q_d_theta   : %.4f\n', x_best(3));
fprintf('   q_d_alpha   : %.4f\n', x_best(4));
fprintf('   R           : %.4f\n', x_best(5));

%% 4. 최종 K 계산 -----------------------------------------------------------
[A, B] = get_sys();

Qf = diag([x_best(1), x_best(2), x_best(3), x_best(4)]);
Rf = x_best(5);

K = lqr(A, B, Qf, Rf);

fprintf('\n최적의 LQR Gain K:\n');
disp(K);

save('K.mat','K');   % Simulink에서 load('K.mat')로 읽으면 됨
fprintf('K 저장 완료! Simulink에서 사용 가능합니다.\n');

%% ========================================================================
% Local Functions
% ========================================================================

%% 시스템 방정식 -----------------------------------------------------------
function [A, B] = get_sys()
    % Q-servo / QUBE pendulum 파라미터 (연구실 세팅에 맞게 필요시 수정)
    Rm = 8.4;  km = 0.042;
    mr = 0.095; r  = 0.085;
    mp = 0.024; Lp = 0.129;
    g  = 9.81;  l  = Lp/2;

    Jr = mr*r^2/3;
    Jp = mp*Lp^2/3;

    M  = [Jr,         mp*l*r;
          mp*l*r,     Jp      ];

    C  = [km^2/Rm, 0;
          0,       0];

    G  = [  0,      0;
           0,   -mp*g*l];

    invM = inv(M);

    A = zeros(4,4);
    A(1,3) = 1; A(2,4) = 1;
    A(3:4,1:2) = -invM*G;
    A(3:4,3:4) = -invM*C;

    B = [0; 0; invM*[km/Rm; 0]];
end

%% 비용 함수 (이름은 standup이지만, 현재는 "balance 성능" 평가) -------------
function J = lqr_cost_standup(param)
    [A, B] = get_sys();

    Q = diag(param(1:4));
    R = param(5);

    % LQR 해가 안 나오면 큰 penalty
    try
        K = lqr(A,B,Q,R);
    catch
        J = 1e6;
        return;
    end

    % ------------ 시뮬레이션 설정 ------------
    dt = 0.002;          % 시뮬레이션 time step [s]
    T  = 2.5;            % 전체 시뮬레이션 시간 [s]
    N  = round(T/dt);

    Vmax = 10;           % 실제 Q-servo 전압 제한 (±10V 가정, 다르면 수정)

    % 속도 추정기(LPF) 설정: G(s) = wc / (s + wc)
    wc = 50;             % rad/s  (Simulink Transfer Fcn의 50/(s+50)과 통일)
    a  = exp(-wc*dt);    % 이산화했을 때 계수 (대략적인 근사)

    % ------------ 초기 상태 ------------
    % 손으로 inverted 근처로 올려놨다고 가정 (예: 8도 정도 오차)
    theta0 = 0;
    alpha0 = deg2rad(8);   % pendulum angle error
    x_true = [theta0; alpha0; 0; 0];

    % 측정 상태 (각도는 바로, 속도는 LPF로 추정)
    theta_dot_meas = 0;
    alpha_dot_meas = 0;
    x_meas = [x_true(1); x_true(2); theta_dot_meas; alpha_dot_meas];

    x_ref = [0; 0; 0; 0];

    % ------------ 비용 계산 ------------
    J = 0;

    for k = 1:N
        % LQR 입력: 측정 상태 기반
        u_raw = -K * (x_meas - x_ref);

        % 전압 제한
        u = max(-Vmax, min(Vmax, u_raw));

        % 실제 시스템 상태 업데이트 (이상적인 연속 모델)
        xdot_true = A * x_true + B * u;
        x_true    = x_true + xdot_true * dt;

        % 실제 각도 / 속도
        theta_true     = x_true(1);
        alpha_true     = x_true(2);
        theta_dot_true = x_true(3);
        alpha_dot_true = x_true(4);

        % ---- 속도 추정기 (Simulink의 angle→Derivative→40/(s+40) 근사) ----
        % 간단한 1차 IIR: v_meas(k+1) = a*v_meas(k) + (1-a)*v_true(k)
        theta_dot_meas = a*theta_dot_meas + (1-a)*theta_dot_true;
        alpha_dot_meas = a*alpha_dot_meas + (1-a)*alpha_dot_true;

        % 측정 상태 갱신
        x_meas = [theta_true; alpha_true; theta_dot_meas; alpha_dot_meas];

        % ---- 비용 누적 (각도 오차 + 입력 크기) ----
        e_theta = theta_true;
        e_alpha = alpha_true;

        w_alpha = 20;
        w_theta = 5;
        w_u     = 0.2;

        J = J + ( ...
            w_alpha*abs(e_alpha) + ...
            w_theta*abs(e_theta) + ...
            w_u    *abs(u) ) * dt;
    end

    % 마지막에 각도가 0 근처에 있지 않으면 큰 penalty
    J = J + 1000*abs(x_true(2)) + 200*abs(x_true(1));
end

%% (참고) Settling time 계산용 함수 - 이건 Simulink MATLAB Function 블록용 예시
%  GA 최적화에는 쓰지 않고, 나중에 Simulink에서 따로 사용하면 됨.
function t_set = settling_time(alpha, t)
persistent settled
persistent t_last

if isempty(settled)
    settled = false;
    t_last = 0;
end

% 목표: alpha = 0 rad, ±2도 이내에 들어간 최초 시간 기억
if abs(alpha) < deg2rad(2)  % ±2°
    if ~settled
        settled = true;
        t_last  = t;
    end
end

t_set = t_last;
end
