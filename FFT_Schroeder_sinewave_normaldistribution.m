clear; clc; close all;

%% 1. 설정
% --------------------------------------------------------
% skew_mode 설정
%  0 : 정규분포 (Skewness = 0)
% -5 : Negative Skew (데이터는 고주파 밀집, Mean은 저주파 쪽으로 이동)
%  5 : Positive Skew (데이터는 저주파 밀집, Mean은 고주파 쪽으로 이동)
alpha = 5; 

f_min = 0.0002; 
f_max = 100;    
n_sines = 50; % 분포를 부드럽게 보기 위해 개수를 늘림
dt = 0.0;
t_end = 5000;

% 저장 경로 (기존 경로 유지)
save_dir = 'G:\공유 드라이브\Battery Software Group (2025)\Internship\25년_동계인턴\전기영\Figure\6주차 FFT\Skew_Normal_Analysis';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% 2. Skew Normal Distribution 기반 가중치 생성
% --------------------------------------------------------
% 주파수를 로그 스케일로 변환하여 계산
freq_vec = logspace(log10(f_min), log10(f_max), n_sines)';
x = log10(freq_vec);

% 분포의 파라미터 (로그 스케일 기준)
mu = (log10(f_min) + log10(f_max)) / 2; % 중심 (0.01Hz 근처)
sigma = 0.5; % 분포의 너비

% Skew Normal PDF 수식: 2 * phi(x) * Phi(alpha * x)
% phi: 표준정규분포 PDF, Phi: 표준정규분포 CDF
z = (x - mu) / sigma;
pdf_normal = (1/(sigma*sqrt(2*pi))) * exp(-z.^2 / 2);
cdf_skew = 0.5 * (1 + erf(alpha * z / sqrt(2)));
weights = 2 * pdf_normal .* cdf_skew;

% 진폭 정규화 (최대값 1A)
amplitudes = weights / max(weights);

% 평균 주파수 계산 (Log-mean)
log_mean = sum(x .* weights) / sum(weights);
mean_freq = 10^log_mean;

%% 3. 신호 합성 (Schroeder Phase)
% --------------------------------------------------------
k_idx = (1:n_sines)';
phases = - (k_idx .* (k_idx - 1) * pi) / n_sines;
t_vec = (0:dt:t_end)';
I_multi = zeros(size(t_vec));

for i = 1:n_sines
    I_multi = I_multi + amplitudes(i) * cos(2*pi*freq_vec(i)*t_vec + phases(i));
end

%% 4. 시각화
% --------------------------------------------------------
hFig = figure('Position', [100, 100, 1000, 800]);

subplot(2, 1, 1);
stem(freq_vec, amplitudes, 'Marker', 'none'); hold on;
plot(freq_vec, amplitudes, 'LineWidth', 2);
xline(mean_freq, 'r--', 'LineWidth', 2, 'Label', sprintf('Mean: %.4f Hz', mean_freq));
set(gca, 'XScale', 'log'); grid on;
xlabel('Frequency [Hz]'); ylabel('Normalized Amplitude');
title(sprintf('Skew Normal Distribution (Alpha = %d)', alpha));

subplot(2, 1, 2);
plot(t_vec, I_multi); grid on;
xlabel('Time [s]'); ylabel('Current [A]');
title('Time Domain Profile');
xlim([0, 1000]);

fprintf('--- 분석 결과 ---\n');
fprintf(' >> 설정 왜도(Alpha): %d\n', alpha);
fprintf(' >> 계산된 평균 주파수: %.5f Hz\n', mean_freq);