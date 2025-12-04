%% QUBE-Servo 2: 150-Degree Step Response Optimization
clear; clc; close all;

% =========================================================
% 1. 파라미터 정의 (QUBE-Servo 2 정밀 데이터)
% =========================================================
val_Rm = 8.4; val_km = 0.042;
val_mr = 0.095; val_r = 0.085;
val_mp = 0.024; val_Lp = 0.129;
val_g  = 9.81;
val_l  = val_Lp / 2;
val_Jr = val_mr * val_r^2 / 3;
val_Jp = val_mp * val_Lp^2 / 3;

% =========================================================
% 2. 시스템 모델링 (State-Space)
% =========================================================
M_num = [val_Jr, val_mp*val_l*val_r; val_mp*val_l*val_r, val_Jp];
C_num = [val_km^2/val_Rm, 0; 0, 0];
G_num = [0, 0; 0, -val_mp*val_g*val_l]; % Inverted Control (-)

A_ss = zeros(4,4);
A_ss(1,3) = 1; A_ss(2,4) = 1;
invM = inv(M_num);
A_ss(3:4, 1:2) = -invM * G_num;
A_ss(3:4, 3:4) = -invM * C_num;
B_ss = [0; 0; invM * [val_km/val_Rm; 0]];

% =========================================================
% 3. LQR 튜닝 (가장 중요한 부분!)
% =========================================================
% 목표: 150도 회전 시 가장 빠른 Settling Time 달성
% 전략: Theta 가중치를 높여서 빠르게 가되, Alpha 가중치로 중심을 잡음

% [튜닝 가이드]
% q_theta: Arm의 목표 도달 속도 (높을수록 빠름)
% q_alpha: Pendulum 중심 잡기 (너무 낮으면 흔들려서 수렴 늦어짐)
% r_val  : 입력 전압 사용량 (낮을수록 과격하게 움직임 -> 전압 포화 주의)

q_theta = 60;   % Arm 오차 가중치
q_alpha = 50;   % Pendulum 각도 가중치 (넘어지지 않게 꽉 잡아야 함)
r_val   = 0.6;  % 입력 비용 (작을수록 반응 빠름)

Q = diag([q_theta, q_alpha, 0.1, 0.1]); 
R = r_val;
K = lqr(A_ss, B_ss, Q, R);

% =========================================================
% 4. 시뮬레이션: Step Response (0 -> 150 deg)
% =========================================================
% 목표 상태 정의 [Theta_target, Alpha_target, 0, 0]
target_deg = 150;
target_rad = target_deg * pi / 180;
x_ref = [target_rad; 0; 0; 0]; % 목표: Arm 150도, Pendulum 0도(수직)

% 폐루프 시스템 구성: dx = (A - BK)x + BK*x_ref
% (Error Dynamics: e = x - x_ref 를 이용한 등가 모델)
sys_cl = ss(A_ss - B_ss*K, B_ss*K, eye(4), 0);

t = 0:0.001:3; % 3초간 시뮬레이션
[y, t, x] = lsim(sys_cl, repmat(x_ref', length(t), 1), t, [0;0;0;0]);

% 입력 전압 역계산 (u = -K(x - x_ref))
u_history = zeros(length(t), 1);
for i = 1:length(t)
    error_state = x(i,:)' - x_ref;
    u_history(i) = -K * error_state;
end

% =========================================================
% 5. 2% Settling Time 자동 계산
% =========================================================
theta_res = y(:,1); % Arm Angle Response
final_val = target_rad;
tolerance = 0.02 * final_val; % 2% 허용 오차 범위

% 범위 안에 들어온 마지막 순간이 아니라, "들어와서 안 나가는" 시점을 찾음
settled_indices = find(abs(theta_res - final_val) > tolerance);
if isempty(settled_indices)
    ts = 0; % 처음부터 범위 내 (그럴 리 없지만)
else
    last_unsettled_idx = settled_indices(end);
    if last_unsettled_idx == length(t)
        ts = NaN; % 아직 수렴 안 함
    else
        ts = t(last_unsettled_idx + 1);
    end
end

% =========================================================
% 6. 결과 출력 및 시각화
% =========================================================
max_u = max(abs(u_history));
VOLT_LIMIT = 10;

fprintf('\n================================================\n');
fprintf('        Simulation Result (Target: 150 deg)\n');
fprintf('================================================\n');
fprintf('1. 2%% Settling Time :  [ %.4f seconds ]\n', ts);
fprintf('2. Max Voltage       :  [ %.2f V ]\n', max_u);
fprintf('------------------------------------------------\n');

if isnan(ts)
    fprintf(' [실패] 3초 안에 수렴하지 못했습니다. Q를 높이세요.\n');
elseif max_u > VOLT_LIMIT
    fprintf(' [경고] 전압 제한(10V)을 초과했습니다! (현실성 없음)\n');
    fprintf(' -> R 값을 키우거나 q_theta를 줄여야 합니다.\n');
else
    fprintf(' [성공] 전압 범위 내에서 안정적으로 수렴했습니다.\n');
    fprintf(' -> 더 줄이고 싶다면 R을 %.2f -> %.2f 로 낮춰보세요.\n', R, R*0.8);
end
fprintf('================================================\n');

figure('Name', '150 deg Step Response', 'Color', 'w', 'Position', [200 100 800 600]);

subplot(2,1,1);
plot(t, y(:,1)*180/pi, 'b-', 'LineWidth', 2); hold on;
yline(target_deg, 'k--', 'Target 150\circ');
yline(target_deg * 1.02, 'g:', '2% Upper');
yline(target_deg * 0.98, 'g:', '2% Lower');
xline(ts, 'r-', ['Ts = ' num2str(ts, '%.3f') 's']);
grid on;
ylabel('Arm Angle (deg)');
title(['Step Response: 0 \rightarrow 150 deg (Settling Time: ' num2str(ts, '%.3f') 's)']);
legend('Response', 'Target', '2% Bound');

subplot(2,1,2);
plot(t, u_history, 'm-', 'LineWidth', 1.5); hold on;
yline(VOLT_LIMIT, 'r--'); yline(-VOLT_LIMIT, 'r--');
grid on;
xlabel('Time (s)'); ylabel('Voltage (V)');
title(['Control Input (Max: ' num2str(max_u, '%.1f') 'V)']);