# Qservo_interted_pendlum



#-----------------

%% 1. GA 설정 (R2022b 호환성 수정: gaoptimset 사용)
% 유전자: [q_theta, q_alpha, q_d_theta, q_d_alpha, r_val]
nVars = 5; 

% 탐색 범위 (Lower Bound, Upper Bound)
lb = [1,   10,  0,   0,   0.1]; 
ub = [500, 1000, 10,  5,   10];  

% [수정됨] optimoptions 대신 gaoptimset 사용
% 주의: 파라미터 이름이 조금 다릅니다 (MaxGenerations -> Generations, PlotFcn -> PlotFcns)
options = gaoptimset(...
    'PopulationSize', 50, ...        % 개체 수
    'Generations', 30, ...           % 세대 수 (MaxGenerations 아님)
    'Display', 'iter', ...           % 진행 상황 출력
    'PlotFcns', @gaplotbestf, ...    % 수렴 그래프 (PlotFcn 아님)
    'UseParallel', false);           % 병렬 연산

    #-----------
