%% QUBE-Servo 2: 150-Degree Step Response Optimization (Saturation Added)
clear; clc; close all;

% 1. 파라미터 정의
val_Rm = 8.4; val_km = 0.042;
val_mr = 0.095; val_r = 0.085;
val_mp = 0.024; val_Lp = 0.129;
val_g  = 9.81;
val_l  = val_Lp / 2;
val_Jr = val_mr * val_r^2 / 3;
val_Jp = val_mp * val_Lp^2 / 3;

% 2. 시스템 모델링
M_num = [val_Jr, val_mp*val_l*val_r; val_mp*val_l*val_r, val_Jp];
C_num = [val_km^2/val_Rm, 0; 0, 0];
G_num = [0, 0; 0, -val_mp*val_g*val_l]; 

invM = inv(M_num);
A_ss = zeros(4,4);
A_ss(1,3) = 1; A_ss(2,4) = 1;
A_ss(3:4, 1:2) = -invM * G_num;
A_ss(3:4, 3:4) = -invM * C_num;
B_ss = [0; 0; invM * [val_km/val_Rm; 0]];

% =========================================================
% 3. LQR 튜닝 (최적화 제안 값)
% =========================================================
% [전략]
% 150도(2.6rad) 점프 시 전압 15V를 넘지 않으려면 K(1)이 약 5.7 이하여야 함.
% R을 높여서 K(1)을 억제하거나, Q_theta를 적절히 타협해야 함.

% 추천 세팅 (Speed & Safety Balance)
q_theta = 25;    % 너무 크면 초기 전압 포화됨 (2.6 rad * K(1) <= 15V 맞추기)
q_alpha = 10;    % 펜듈럼 균형
q_d_theta = 2;   % 암 속도 댐핑 (오버슈트 억제 -> 정착시간 단축)
q_d_alpha = 0.1; % 펜듈럼 속도 댐핑

r_val   = 3.0;   % 전압을 15V 안쪽으로 누르기 위해 R을 키움

Q = diag([q_theta, q_alpha, q_d_theta, q_d_alpha]); 
R = r_val;

K = lqr(A_ss, B_ss, Q, R);

fprintf('Calculated K(1) (Theta Gain): %.4f\n', K(1));
fprintf('Max Initial Voltage Estimate: %.2f V\n', abs(K(1) * (150*pi/180)));

% =========================================================
% 4. 시뮬레이션: 전압 포화(Saturation) 고려
% =========================================================
target_deg = 150;
target_rad = target_deg * pi / 180;
x_ref = [target_rad; 0; 0; 0]; 

dt = 0.001;
t = 0:dt:3;
x_sim = zeros(length(t), 4); % 상태 저장 변수
u_history = zeros(length(t), 1);

VOLT_LIMIT = 15; % 모터 물리적 한계

% 초기 상태
x_curr = [0; 0; 0; 0]; 
x_sim(1,:) = x_curr';

for i = 1:length(t)-1
    % 1. 에러 계산 및 전압 산출
    error_state = x_curr - x_ref;
    u_raw = -K * error_state;
    
    % 2. [핵심] 전압 포화 (Saturation) 적용
    % 실제 모터는 15V 이상을 낼 수 없음 -> 클램핑
    u_clamped = max(-VOLT_LIMIT, min(VOLT_LIMIT, u_raw));
    u_history(i) = u_clamped;
    
    % 3. 상태 업데이트 (Euler Method for simplicity)
    % dx = Ax + Bu
    dx = A_ss * x_curr + B_ss * u_clamped;
    x_next = x_curr + dx * dt;
    
    % 다음 스텝 준비
    x_sim(i+1, :) = x_next';
    x_curr = x_next;
end
u_history(end) = u_history(end-1); % 마지막 값 채움

% =========================================================
% 5. 2% Settling Time 계산
% =========================================================
y_theta = x_sim(:,1);
final_val = target_rad;
tolerance = 0.02 * final_val; 

% 범위 밖으로 나간 마지막 인덱스 찾기
unsettled_idx = find(abs(y_theta - final_val) > tolerance, 1, 'last');

if isempty(unsettled_idx)
    ts = 0;
else
    ts = t(unsettled_idx + 1);
end

% =========================================================
% 6. 결과 출력
% =========================================================
max_u_sim = max(abs(u_history));

fprintf('\n================================================\n');
fprintf('   Real Simulation Result (Saturation Applied)\n');
fprintf('================================================\n');
fprintf('1. 2%% Settling Time :  [ %.4f seconds ]\n', ts);
fprintf('2. Max Voltage Used  :  [ %.2f V ] (Limit: %dV)\n', max_u_sim, VOLT_LIMIT);
fprintf('------------------------------------------------\n');

% 그래프
figure('Name', 'Realistic Step Response', 'Color', 'w');
subplot(2,1,1);
plot(t, y_theta*180/pi, 'b-', 'LineWidth', 2); hold on;
yline(target_deg, 'k--');
yline(target_deg * 1.02, 'g:');
yline(target_deg * 0.98, 'g:');
xline(ts, 'r-', ['Ts = ' num2str(ts, '%.3f') 's']);
grid on; ylabel('Arm Angle (deg)');
title('Step Response with Voltage Limit');

subplot(2,1,2);
plot(t, u_history, 'm-', 'LineWidth', 1.5); hold on;
yline(VOLT_LIMIT, 'r--'); yline(-VOLT_LIMIT, 'r--');
grid on; xlabel('Time (s)'); ylabel('Voltage (V)');
title('Actual Input Voltage (Clamped)');