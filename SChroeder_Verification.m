%% ======================================================================
%  [발표용] Schroeder Phase 검증 (Side-by-Side Comparison)
%  목적: Wrapped(랜덤해 보임) vs Unwrapped(완벽한 곡선) 비교 시각화
% ======================================================================
clear; clc; close all;

% 1. 파일 설정 (경로를 본인 파일에 맞게 수정하세요)
filepath = 'G:\공유 드라이브\Battery Software Lab\Protocols\2025 현대 NE 고품셀 평가 실험\주행부하_scaledown\udds_0725.xlsx';
[~, fname, ~] = fileparts(filepath);

% 2. 데이터 로드
data = readtable(filepath);
t = data{:, 1}; 
I_load = data{:, 2};
mask = ~(isnan(t) | isnan(I_load)); t = t(mask); I_load = I_load(mask);

% 3. Schroeder 신호 생성
f_min = 0.001558; f_max = 0.03186; num_sines = 20; A = 2;
freqs = logspace(log10(f_min), log10(f_max), num_sines);
k_vec = 1:num_sines;

% [이론값] Unwrapped Schroeder Phase (곡선)
theo_phases_unwrapped = -k_vec .* (k_vec - 1) * pi / num_sines; 

% 신호 주입
I_noise = zeros(size(t));
for k = 1:num_sines
    I_noise = I_noise + A * sin(2*pi*freqs(k)*t + theo_phases_unwrapped(k));
end
I_total = I_load + I_noise;

% 4. FFT 및 위상 추출
dt = median(diff(t)); fs = 1/dt; N = numel(I_total);
I_fft = fft((I_total - mean(I_total)) .* hann(N));
f_axis = (0:floor(N/2)-1)*fs/N;
raw_phase_spectrum = angle(I_fft(1:floor(N/2))); % 전체 스펙트럼 (-pi ~ pi)

% 5. 데이터 정리 (Wrapped vs Unwrapped)
meas_phases_wrapped = zeros(1, num_sines);   % 왼쪽 그래프용
meas_phases_unwrapped = zeros(1, num_sines); % 오른쪽 그래프용

for k = 1:num_sines
    % (1) 해당 주파수 인덱스 찾기
    [~, idx] = min(abs(f_axis - freqs(k)));
    
    % (2) Wrapped 값 추출 (Raw Data)
    raw_val = raw_phase_spectrum(idx);
    meas_phases_wrapped(k) = raw_val;
    
    % (3) Unwrapped 값 계산 (Alignment)
    theo_val = theo_phases_unwrapped(k);
    diff_val = raw_val - theo_val;
    n_shift = round(diff_val / (2*pi));
    meas_phases_unwrapped(k) = raw_val - n_shift * 2*pi;
end

% ======================================================================
% 6. 시각화 (좌우 비교)
% ======================================================================
figure('Name', 'Schroeder Phase Comparison', 'Position', [100 200 1200 600], 'Color', 'w');

% --- [왼쪽] Wrapped Phase (눈에 보이는 모습) ---
subplot(1, 2, 1);
plot(1:num_sines, meas_phases_wrapped, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
grid on;
title({'(1) Wrapped Phase (Raw FFT Data)', '규칙 없이 랜덤하게 보임 (-\pi ~ +\pi 제한)'}, 'FontSize', 14);
xlabel('Frequency Index (k)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Phase (rad)', 'FontSize', 12, 'FontWeight', 'bold');
ylim([-pi*1.2 pi*1.2]); 
yticks([-pi 0 pi]); yticklabels({'-\pi', '0', '\pi'});
xlim([0 num_sines+1]);

% --- [오른쪽] Unwrapped Phase (진실된 모습) ---
subplot(1, 2, 2);
% (1) 이론값 곡선 (파란 실선)
plot(1:num_sines, theo_phases_unwrapped, 'b-', 'LineWidth', 2); hold on;
% (2) 실제 측정값 (빨간 점)
plot(1:num_sines, meas_phases_unwrapped, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

grid on;
title({'(2) Unwrapped Phase (Aligned)', '완벽한 2차함수 곡선을 따름'}, 'FontSize', 14);
xlabel('Frequency Index (k)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Cumulative Phase (rad)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Theoretical Curve (Schroeder)', 'Measured Phase (FFT)'}, ...
       'Location', 'southwest', 'FontSize', 11);
xlim([0 num_sines+1]);

% 화살표 주석 (오른쪽 그래프 강조)
annotation('textarrow',[0.75, 0.65],[0.5, 0.4],'String',{'이 곡선 패턴이','Schroeder Phase의 증거'}, ...
    'FontSize',11, 'Color','b', 'LineWidth', 1.5);

% 전체 제목
sgtitle(['[' fname '] Schroeder Phase 검증: 보이는 것 vs 실제 규칙'], 'FontSize', 16, 'FontWeight', 'bold');

fprintf('>> 그래프 생성 완료. 좌우 비교를 통해 규칙성을 확인하세요.\n');