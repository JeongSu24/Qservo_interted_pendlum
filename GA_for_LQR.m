%% GA-based LQR Optimization for Rotary Inverted Pendulum (Integrated)
% 작성자: 제어공학 연구실 (Prof. Gemini & Student)
% 기능: 유전 알고리즘을 이용해 Q, R을 튜닝하고 최적의 K를 Simulink로 내보냄
clear; clc; close all;

fprintf('================================================\n');
fprintf('   Genetic Algorithm LQR Tuning Start\n');
fprintf('================================================\n');

%% 1. GA 설정
% 유전자: [q_theta, q_alpha, q_d_theta, q_d_alpha, r_val]
nVars = 5; 

% 탐색 범위 (Lower Bound, Upper Bound)
% 전략: q_theta(반응속도), q_alpha(균형)은 넓게, 속도항은 좁게, R은 적절히
lb = [1,   10,  0,   0,   0.1]; 
ub = [500, 1000, 10,  5,   10];  

% GA 옵션 설정
options = optimoptions('ga', ...
    'PopulationSize', 50, ...        % 개체 수 (많을수록 전역 최적해 찾을 확률 높음)
    'MaxGenerations', 30, ...        % 세대 수
    'Display', 'iter', ...           % 진행 상황 출력
    'PlotFcn', @gaplotbestf, ...     % 수렴 그래프 표시
    'UseParallel', false);           % 병렬 연산 (필요시 true)

%% 2. GA 실행 (Optimization)
% 비용 함수 핸들(@cost_func)을 넘겨줌
fprintf('최적화를 시작합니다. 잠시만 기다리게...\n');
[x_best, fval] = ga(@lqr_cost_function, nVars, [], [], [], [], lb, ub, [], options);

%% 3. 최적 결과 출력
fprintf('\n================================================\n');
fprintf('   Optimization Completed!\n');
fprintf('================================================\n');
fprintf('Best Parameters Found:\n');
fprintf(' - q_theta (Arm Pos)    : %.4f\n', x_best(1));
fprintf(' - q_alpha (Pend Angle) : %.4f\n', x_best(2));
fprintf(' - q_d_theta (Arm Vel)  : %.4f\n', x_best(3));
fprintf(' - q_d_alpha (Pend Vel) : %.4f\n', x_best(4));
fprintf(' - R (Input Cost)       : %.4f\n', x_best(5));
fprintf('------------------------------------------------\n');
fprintf('Best Cost (Score)      : %.4f\n', fval);

%% 4. 최종 검증 및 시각화 (Verification)
% 최적 파라미터로 정밀 시뮬레이션을 수행하고 그래프를 그립니다.
verify_and_plot(x_best);

%% 5. Simulink 연동을 위한 K값 확정 (Export)
% 최적의 파라미터로 K를 다시 계산해서 워크스페이스에 남깁니다.
[A_ss, B_ss] = get_sys_model(); % 시스템 모델 불러오기
Q_final = diag([x_best(1), x_best(2), x_best(3), x_best(4)]);
R_final = x_best(5);

K = lqr(A_ss, B_ss, Q_final, R_final);

fprintf('\n================================================\n');
fprintf('   Export to Simulink\n');
fprintf('================================================\n');
fprintf('최적화된 Gain K가 생성되었습니다:\n');
disp(K);
fprintf('이제 Simulink 파일에서 [Run] 버튼을 누르면 이 K값이 적용됩니다.\n');


%% ========================================================================
%  Local Functions (내부 함수 정의)
%  ========================================================================

% [함수 1] 시스템 모델링 (중복 방지를 위해 별도 함수로 분리)
function [A, B] = get_sys_model()
    % 파라미터 정의
    val_Rm = 8.4; val_km = 0.042;
    val_mr = 0.095; val_r = 0.085;
    val_mp = 0.024; val_Lp = 0.129;
    val_g  = 9.81; val_l  = val_Lp / 2;
    val_Jr = val_mr * val_r^2 / 3;
    val_Jp = val_mp * val_Lp^2 / 3;

    % 운동 방정식 계수 행렬
    M_num = [val_Jr, val_mp*val_l*val_r; val_mp*val_l*val_r, val_Jp];
    C_num = [val_km^2/val_Rm, 0; 0, 0];
    G_num = [0, 0; 0, -val_mp*val_g*val_l]; % Inverted Pendulum (-)

    invM = inv(M_num);
    
    % State Space 행렬 구성
    A = zeros(4,4);
    A(1,3) = 1; A(2,4) = 1;
    A(3:4, 1:2) = -invM * G_num;
    A(3:4, 3:4) = -invM * C_num;
    
    B = [0; 0; invM * [val_km/val_Rm; 0]];
end

% [함수 2] 비용 함수 (GA가 호출함)
function cost = lqr_cost_function(params)
    [A_ss, B_ss] = get_sys_model(); % 모델 가져오기
    
    % 파라미터 언패킹
    Q = diag([params(1), params(2), params(3), params(4)]);
    R = params(5);
    
    % LQR 계산 (실패 시 페널티)
    try
        K = lqr(A_ss, B_ss, Q, R);
    catch
        cost = 1e6; return; 
    end

    % 시뮬레이션 설정
    dt = 0.005; % GA 속도를 위해 dt는 5ms
    t_end = 2.0; 
    target_rad = 150 * pi / 180;
    
    x_curr = [0; 0; 0; 0];
    x_ref = [target_rad; 0; 0; 0];
    
    max_u = 0;
    hist_theta = zeros(floor(t_end/dt), 1);
    
    % 비선형 시뮬레이션 루프
    for i = 1:length(hist_theta)
        u_raw = -K * (x_curr - x_ref);
        % 전압 포화 (15V)
        u_clamped = max(-15, min(15, u_raw));
        if abs(u_raw) > max_u, max_u = abs(u_raw); end
        
        dx = A_ss * x_curr + B_ss * u_clamped;
        x_curr = x_curr + dx * dt;
        hist_theta(i) = x_curr(1);
    end

    % 점수 계산 (Cost Calculation)
    % 1. 정착 시간 (Settling Time)
    tol = 0.02 * target_rad;
    err = abs(hist_theta - target_rad);
    unsettled_idx = find(err > tol, 1, 'last');
    
    if isempty(unsettled_idx), ts = 0;
    else, ts = unsettled_idx * dt;
    end
    
    % 2. 각종 페널티
    penalty_unstable = 0;
    if ts >= t_end - dt, penalty_unstable = 2000; end % 수렴 안 하면 큰 벌점
    
    penalty_voltage = 0;
    if max_u > 15, penalty_voltage = (max_u - 15) * 50; end % 전압 초과 벌점 강화
    
    max_overshoot = max(hist_theta);
    penalty_overshoot = 0;
    if max_overshoot > target_rad * 1.05 % 5% 오버슈트 허용
        penalty_overshoot = (max_overshoot - target_rad) * 1000; 
    end

    % 최종 비용 (Settling Time을 줄이는 게 1순위)
    cost = (ts * 100) + penalty_voltage + penalty_unstable + penalty_overshoot;
end

% [함수 3] 결과 검증 및 그래프 (최종 확인용)
function verify_and_plot(params)
    [A_ss, B_ss] = get_sys_model();
    Q = diag(params(1:4)); R = params(5);
    K = lqr(A_ss, B_ss, Q, R);
    
    dt = 0.001; % 검증 때는 정밀하게 (1ms)
    t = 0:dt:3;
    target_deg = 150;
    target_rad = target_deg * pi / 180;
    
    x_curr = [0; 0; 0; 0];
    x_ref = [target_rad; 0; 0; 0];
    
    x_hist = zeros(length(t), 4);
    u_hist = zeros(length(t), 1);
    
    for i = 1:length(t)
        u_raw = -K * (x_curr - x_ref);
        u_clamped = max(-15, min(15, u_raw));
        
        u_hist(i) = u_clamped;
        dx = A_ss * x_curr + B_ss * u_clamped;
        x_curr = x_curr + dx * dt;
        x_hist(i,:) = x_curr';
    end
    
    % 그래프 그리기
    figure('Name', 'Optimized LQR Result', 'Color', 'w', 'Position', [100 100 800 600]);
    
    subplot(2,1,1);
    plot(t, x_hist(:,1)*180/pi, 'b-', 'LineWidth', 2); hold on;
    plot(t, x_hist(:,2)*180/pi, 'r:', 'LineWidth', 1.5);
    yline(target_deg, 'k--');
    legend('Arm (\theta)', 'Pendulum (\alpha)', 'Target');
    title('Step Response (0 \rightarrow 150 deg)');
    ylabel('Angle (deg)'); grid on;
    
    subplot(2,1,2);
    plot(t, u_hist, 'm-', 'LineWidth', 1.5);
    yline(15, 'r--'); yline(-15, 'r--');
    xlabel('Time (s)'); ylabel('Voltage (V)');
    title(['Control Input (R = ' num2str(R, '%.3f') ')']);
    grid on;
end
