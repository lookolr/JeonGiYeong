clear; clc; close all;

%% 1. 파라미터 설정
% --------------------------------------------------------
target_mean_freq = 0.00979;   % [Hz] 고정된 중심 주파수 (Mean)
base_f_min       = 0.0001;    % [Hz] 기준이 되는 원래 최소 주파수
base_f_max       = 100;       % [Hz] 기준이 되는 원래 최대 주파수

% [핵심] 확장 비율 설정 (1.0 = 유지, 1.2 = 20% 확장, 1.5 = 50% 확장)
expansion_ratio  = 1.0;       

n_sines          = 20;        % Sine wave 개수
dt               = 1;     % 샘플링 시간
t_end            = 6000;      % 총 시간 (데이터 개수: 6,000,001개)

% 저장 경로
save_dir = 'G:\공유 드라이브\Battery Software Group (2025)\Internship\25년_동계인턴\전기영\Figure\6주차 FFT\99_Sets';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% 2. 확장된 주파수 범위 계산
% --------------------------------------------------------
% (1) 원래 범위의 로그 반지름(반폭) 계산
original_radius_log = log10(base_f_max) - log10(target_mean_freq);

% (2) 확장된 반지름 계산 (비율 적용)
expanded_radius_log = original_radius_log * expansion_ratio;

% (3) 새로운 최소/최대 주파수 도출 (중심 기준 대칭 확장)
f_min_new = target_mean_freq / 10^(expanded_radius_log);
f_max_new = target_mean_freq * 10^(expanded_radius_log);

% (4) 주파수 벡터 생성
freq_vec = logspace(log10(f_min_new), log10(f_max_new), n_sines)';

%% 3. 분포 생성 (Sigma 재계산)
% --------------------------------------------------------
% 확장된 범위 안에 예쁜 종 모양이 꽉 차도록 Sigma를 범위에 맞춰 조정
current_radius_log = (log10(f_max_new) - log10(f_min_new)) / 2;
calculated_sigma   = current_radius_log / 2.5; 

% 가우시안 가중치 계산
mu_log = log10(target_mean_freq);
x_log  = log10(freq_vec);
weights = exp( - (x_log - mu_log).^2 / (2 * calculated_sigma^2) );
amplitudes = weights / max(weights); 

%% 4. 신호 합성
% --------------------------------------------------------
k_idx = (1:n_sines)';
% Schroeder Phase (PAPR 최소화)
phases = - (k_idx .* (k_idx - 1) * pi) / n_sines;

t_vec = (0:dt:t_end)';
I_multi = zeros(size(t_vec));

fprintf('--- 범위 확장 시뮬레이션 ---\n');
fprintf(' >> 고정 중심(Mean): %.5f Hz\n', target_mean_freq);
fprintf(' >> 기존 범위: %.5f ~ %.5f Hz\n', base_f_min, base_f_max);
fprintf(' >> 확장 범위(x%.1f): %.5f ~ %.5f Hz\n', expansion_ratio, f_min_new, f_max_new);
fprintf(' >> 데이터 개수: %d 개 (Excel 저장 불가 -> MAT 저장)\n', length(t_vec));

% 벡터 연산으로 속도 최적화 (For loop 대신)
for i = 1:n_sines
    I_multi = I_multi + amplitudes(i) * cos(2*pi*freq_vec(i)*t_vec + phases(i));
end

%% 5. 시각화
% --------------------------------------------------------
hFig = figure('Position', [100, 100, 1000, 800]);

% [Top] Frequency Domain
subplot(2, 1, 1);
stem(freq_vec, amplitudes, 'filled', 'Color', 'r', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(freq_vec, amplitudes, 'r-', 'LineWidth', 2);

% 기준선들 표시
xline(target_mean_freq, 'b--', 'LineWidth', 2, 'Label', 'Center (Fixed)');
xline(base_f_min, 'k:', 'LineWidth', 1.5, 'Label', 'Original Min');
xline(base_f_max, 'k:', 'LineWidth', 1.5, 'Label', 'Original Max');

% 새로 확장된 범위 표시
xline(f_min_new, 'g-', 'LineWidth', 1, 'Label', 'New Min');
xline(f_max_new, 'g-', 'LineWidth', 1, 'Label', 'New Max');

set(gca, 'XScale', 'log'); 
grid on;
xlabel('Frequency [Hz] (Log Scale)');
ylabel('Amplitude (Weight)');
title(sprintf('Expanded Range (Ratio: %.1f, Mean Fixed)', expansion_ratio));
xlim([f_min_new*0.5, f_max_new*2]); % 시야를 조금 더 넓게
ylim([0, 1.3]);

% [Bottom] Time Domain
subplot(2, 1, 2);
plot(t_vec, I_multi, 'k');
grid on;
xlabel('Time [s]');
ylabel('Current [A]');
title('Generated Current Profile');
xlim([0, 100]); % 앞부분 100초만 확대해서 보여줌 (전체는 너무 빽빽함)

%% 6. 저장 (수정됨: XLSX -> MAT)
% --------------------------------------------------------
file_tag = sprintf('MultiSine_Expanded_Ratio%.1f', expansion_ratio);

% [중요] .mat 파일로 저장 (대용량 데이터 처리용)
full_mat_path = fullfile(save_dir, [file_tag '.mat']);

% 필요한 변수들만 골라서 저장
save(full_mat_path, 't_vec', 'I_multi', 'freq_vec', 'amplitudes', 'dt', 't_end');

fprintf(' >> 데이터 저장 완료 (MAT): %s\n', full_mat_path);

% 그림 저장
full_fig_path = fullfile(save_dir, [file_tag '.png']);
saveas(hFig, full_fig_path);
fprintf(' >> 그림 저장 완료 (PNG): %s\n', full_fig_path);