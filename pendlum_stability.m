%% Rotary Inverted Pendulum: Pole & Zero Count Check
clear; clc; close all;

% 1. 파라미터 정의
syms Jr Jp mr mp r Lp l g km Rm s real
% Jr: Arm Inertia, Jp: Pendulum Inertia
% mr: Arm Mass, mp: Pendulum Mass
% r: Arm Length, Lp: Pendulum Total Length, l: Pendulum CoM Length

% 2. 파라미터 값 정의 (Table 2.2: QUBE-Servo 2 Parameters 기반)
val_Rm = 8.4;             % Terminal resistance [Ohm]
val_km = 0.042;           % Back-emf constant / Torque constant [V/(rad/s) or N.m/A]
val_mr = 0.095;           % Rotary arm mass [kg]
val_r  = 0.085;           % Rotary arm length (pivot to tip) [m]
val_mp = 0.024;           % Pendulum link mass [kg]
val_Lp = 0.129;           % Pendulum link total length [m]
val_g  = 9.81;            % Gravity [m/s^2]

% 3. 유도된 파라미터 계산 (Equation of Motion 이미지 기반)
% (1) l: 펜듈럼 무게중심 거리 (균일한 막대 가정 시 전체 길이의 절반)
val_l  = val_Lp / 2;      
% (2) Jr: 회전 팔 관성 모멘트 (공식: mr * r^2 / 3)
val_Jr = val_mr * val_r^2 / 3; 
% (3) Jp: 펜듈럼 관성 모멘트 (공식: mp * Lp^2 / 3)
val_Jp = val_mp * val_Lp^2 / 3;

% 값과 심볼 매핑
vals = {Jr, Jp, mr, mp, r, Lp, l, g, km, Rm};
nums = {val_Jr, val_Jp, val_mr, val_mp, val_r, val_Lp, val_l, val_g, val_km, val_Rm};

% 4. 운동 방정식 행렬 (제공된 식: Linearized Equation 기반)
% 식 1: Jr*theta_dd + mp*l*r*alpha_dd = tau - br*theta_d (br 무시)
%       tau = (km/Rm)*(Vm - km*theta_d)
%       -> (Jr*s^2 + (km^2/Rm)*s)*Theta + (mp*l*r*s^2)*Alpha = (km/Rm)*Vm

M11 = Jr*s^2 + (km^2/Rm)*s;
M12 = mp*l*r*s^2;

% 식 2: Jp*alpha_dd + mp*l*r*theta_dd + mp*g*l*alpha = -bp*alpha_d (bp 무시)
%       -> (mp*l*r*s^2)*Theta + (Jp*s^2 + mp*g*l)*Alpha = 0
% [주의!] 이미지 식에 따라 mp*g*l 부호를 (+)로 설정함. (Pendant type)
% 만약 역진자(Inverted) 제어라면 'M22 = Jp*s^2 - mp*g*l'로 변경해야 함.

M21 = mp*l*r*s^2;
M22 = Jp*s^2 - mp*g*l;  % 이미지 식 그대로 적용 (+)

A_mat = [M11, M12; M21, M22];
B_vec = [(km/Rm); 0];

% 5. 전달함수 유도 및 변환
X = A_mat \ B_vec;
TF_Theta = simplify(X(1)); % Theta (Arm)
TF_Alpha = simplify(X(2)); % Alpha (Pendulum)

[n_th, d_th] = numden(TF_Theta);
sys_theta = tf(double(sym2poly(subs(n_th, vals, nums))), ...
               double(sym2poly(subs(d_th, vals, nums))));

[n_al, d_al] = numden(TF_Alpha);
sys_alpha = tf(double(sym2poly(subs(n_al, vals, nums))), ...
               double(sym2poly(subs(d_al, vals, nums))));

% 6. 결과 확인 (Pole/Zero)
fprintf('\n=== QUBE-Servo 2 Model Parameters ===\n');
fprintf('Calculated Jr: %.3e kg.m^2\n', val_Jr);
fprintf('Calculated Jp: %.3e kg.m^2\n', val_Jp);
fprintf('Calculated l (CoM): %.4f m\n', val_l);
fprintf('--------------------------------------\n');
fprintf('System Poles (Stability Check):\n');
pole(sys_alpha)

% 7. 시각화
figure('Color','w');
pzmap(sys_alpha);
title('Pole-Zero Map (QUBE-Servo 2 Model)');
grid on;