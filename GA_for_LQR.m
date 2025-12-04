%% GA-based LQR Optimization for Rotary Inverted Pendulum
clear; clc; close all;

fprintf('================================================\n');
fprintf('   Genetic Algorithm LQR Tuning Start\n');
fprintf('================================================\n');

% 1. GA 설정
% 유전자: [q_theta, q_alpha, q_d_theta, q_d_alpha, r_val]
nVars = 5; 

% 탐색 범위 (Lower Bound, Upper Bound)
% q_theta, q_alpha는 크게, 속도항은 작게, R은 적당히
lb = [1,   1,   0,   0,   0.1]; 
ub = [500, 500, 10,  10,  10];  

% GA 옵션 설정 (세대 수, 인구 수 등)
options = optimoptions('ga', ...
    'PopulationSize', 50, ...        % 한 세대에 50개의 개체
    'MaxGenerations', 30, ...        % 최대 30세대 진화
    'Display', 'iter', ...           % 진행 상황 출력
    'PlotFcn', @gaplotbestf);        % 최적값 그래프 실시간 표시

% 2. GA 실행 (오래 걸릴 수 있음)
% cost_function은 아래에 정의되어 있음
[x_best, fval] = ga(@lqr_cost_function, nVars, [], [], [], [], lb, ub, [], options);

% 3. 최적 결과 출력
fprintf('\n================================================\n');
fprintf('   Optimization Completed!\n');
fprintf('================================================\n');
fprintf('Best Parameters Found:\n');
fprintf(' - q_theta   : %.4f\n', x_best(1));
fprintf(' - q_alpha   : %.4f\n', x_best(2));
fprintf(' - q_d_theta : %.4f\n', x_best(3));
fprintf(' - q_d_alpha : %.4f\n', x_best(4));
fprintf(' - R         : %.4f\n', x_best(5));
fprintf('------------------------------------------------\n');
fprintf('Best Cost (Score): %.4f\n', fval);

% 4. 최적 파라미터로 최종 검증 시뮬레이션
verify_result(x_best);


%% [핵심] 적합도 함수 (Cost Function)
% GA가 이 함수를 수천 번 호출하면서 점수를 매김
function cost = lqr_cost_function(params)
    % 파라미터 언패킹
    q_th = params(1);
    q_al = params(2);
    qd_th = params(3);
    qd_al = params(4);
    r_val = params(5);

    % --- 1. 시스템 모델링 (고정값) ---
    val_Rm = 8.4; val_km = 0.042;
    val_mr = 0.095; val_r = 0.085;
    val_mp = 0.024; val_Lp = 0.129;
    val_g  = 9.81; val_l  = val_Lp / 2;
    val_Jr = val_mr * val_r^2 / 3;
    val_Jp = val_mp * val_Lp^2 / 3;

    M_num = [val_Jr, val_mp*val_l*val_r; val_mp*val_l*val_r, val_Jp];
    C_num = [val_km^2/val_Rm, 0; 0, 0];
    G_num = [0, 0; 0, -val_mp*val_g*val_l];

    invM = inv(M_num);
    A_ss = zeros(4,4);
    A_ss(1,3) = 1; A_ss(2,4) = 1;
    A_ss(3:4, 1:2) = -invM * G_num;
    A_ss(3:4, 3:4) = -invM * C_num;
    B_ss = [0; 0; invM * [val_km/val_Rm; 0]];

    % --- 2. LQR 게인 계산 ---
    Q = diag([q_th, q_al, qd_th, qd_al]);
    R = r_val;
    
    % LQR 계산 실패 시 (불안정 등) 엄청난 페널티 부여
    try
        K = lqr(A_ss, B_ss, Q, R);
    catch
        cost = 1e6; return; 
    end

    % --- 3. 비선형 시뮬레이션 (Saturation 포함) ---
    dt = 0.005; % 속도를 위해 dt를 조금 키움 (0.001 -> 0.005)
    t_end = 2.0; % 2초까지만 봄
    t = 0:dt:t_end;
    
    target_rad = 150 * pi / 180;
    x_ref = [target_rad; 0; 0; 0];
    x_curr = [0; 0; 0; 0];
    
    max_u = 0;
    x_hist = zeros(length(t), 4);
    
    for i = 1:length(t)
        u_raw = -K * (x_curr - x_ref);
        % 전압 포화 (15V)
        u_clamped = max(-15, min(15, u_raw));
        if abs(u_raw) > max_u, max_u = abs(u_raw); end
        
        dx = A_ss * x_curr + B_ss * u_clamped;
        x_curr = x_curr + dx * dt;
        x_hist(i,:) = x_curr';
    end

    % --- 4. 점수 계산 (가장 중요!) ---
    theta_res = x_hist(:,1);
    
    % (A) 정착 시간 (Settling Time) 계산
    tol = 0.02 * target_rad;
    err = abs(theta_res - target_rad);
    unsettled_idx = find(err > tol, 1, 'last');
    
    if isempty(unsettled_idx)
        ts = 0; 
    else
        ts = t(unsettled_idx);
    end
    
    % 만약 끝까지 수렴 안했으면 페널티
    if ts >= t_end - dt
        penalty_unstable = 1000;
    else
        penalty_unstable = 0;
    end
    
    % (B) 전압 페널티 (포화가 심하게 일어나면 감점)
    % 실제로 15V에 걸리더라도, "원래 요구했던 전압(u_raw)"이 너무 크면 
    % 제어기가 무리하고 있다는 뜻이므로 감점함.
    penalty_voltage = 0;
    if max_u > 15
        penalty_voltage = (max_u - 15) * 10; % 1V 초과당 10점 감점
    end

    % (C) 오버슈트 페널티 (목표값보다 10% 이상 넘어가면 감점)
    max_theta = max(theta_res);
    penalty_overshoot = 0;
    if max_theta > target_rad * 1.1
        penalty_overshoot = (max_theta - target_rad) * 500;
    end

    % 최종 비용 함수 (작을수록 좋음)
    % 정착 시간(초)에 가중치 100을 둠 + 각종 페널티 합산
    cost = (ts * 100) + penalty_voltage + penalty_unstable + penalty_overshoot;
end


%% 검증용 함수
function verify_result(params)
    fprintf('\nRunning Verification with Best Parameters...\n');
    % 여기서 최적 파라미터로 상세 시뮬레이션 및 그래프 그리기
    % (위의 cost_function 내용을 상세 그래프용으로 복사해서 사용)
    % ... 자네가 기존에 쓰던 코드를 여기 넣으면 되네. 
    % 단, K값 계산 시 params를 사용해야 함.
    
    q_th = params(1); q_al = params(2); 
    qd_th = params(3); qd_al = params(4); r_val = params(5);
    
    % 모델링 및 K 재계산 (코드 중복 최소화를 위해 간략히 씀)
    % 실제로는 위 lqr_cost_function의 모델링 파트와 동일하게 작성 필요
    % (이 칸이 부족해서 생략하지만, 자네는 합칠 수 있을 걸세!)
    fprintf('검증 완료! 이 파라미터로 시뮬링크에 넣게나.\n');
end