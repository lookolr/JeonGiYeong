%% ======================================================================
%  Standalone Multi-Sine Generator (Skewed Freq Distribution)
%  목적: 주파수 분포를 f_min 또는 f_max 쪽으로 치우치게 생성
% ======================================================================
clear; clc; close all;

% [사용자 설정]
seed_number = 1;        

% -------------------------------------------------------------
% [1] 시뮬레이션 및 파형 기본 설정
% -------------------------------------------------------------
dt          = 1;        % 샘플링 간격 [s]
t_end       = 4000;     % 총 시간 [s]
A_individual = 2;       % 개별 사인파 진폭 [A]
num_sines   = 20;       % 사인파 개수

% [2] 주파수 범위 설정
f_min = 0.000318;      % 최저 주파수
f_max = 0.0386;        % 최고 주파수

% -------------------------------------------------------------
% [3] 주파수 분포 치우침(Skew) 설정 (⭐ 여기가 핵심입니다)
% -------------------------------------------------------------
% skew_power = 1.0  -> 기존과 동일 (Log-scale 등간격)
% skew_power < 1.0  -> f_max (고주파) 쪽에 촘촘하게 쏠림 (예: 0.4 추천)
% skew_power > 1.0  -> f_min (저주파) 쪽에 촘촘하게 쏠림 (예: 2.5 추천)
skew_power = 0.4;   % <--- 이 값을 바꿔가며 테스트해보세요!

% 1. 저장 경로 설정
save_dir = 'G:\공유 드라이브\Battery Software Group (2025)\Internship\25년_동계인턴\전기영\Figure\4주차 FFT\FFT_Schroeder sine wave';

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
fprintf('저장 경로: %s\n', save_dir);

% 2. 시간 벡터 생성
t_vec = (0:dt:t_end)'; 

% 3. 주파수 및 위상 생성 (Skew 적용 + Schroeder Phase)
% (1) Skew된 주파수 벡터 생성 로직
%  - 단계 1: 0부터 1까지 선형 벡터 생성
norm_idx = linspace(0, 1, num_sines);

%  - 단계 2: 거듭제곱을 이용해 간격 왜곡 (Skewing)
skewed_idx = norm_idx .^ skew_power;

%  - 단계 3: Log 스케일 범위에 매핑
log_f_min = log10(f_min);
log_f_max = log10(f_max);
log_freqs = log_f_min + (log_f_max - log_f_min) * skewed_idx;

%  - 단계 4: 최종 주파수 변환
freqs = 10 .^ log_freqs;

% (2) Schroeder Phase 생성
phases = zeros(1, num_sines);
for k = 1:num_sines
    phases(k) = - (k * (k - 1) / num_sines) * pi;
end

% 4. 합성 (Superposition)
I_vec = zeros(size(t_vec));
fprintf('>> Sine Wave 합성 중... (Skew Power: %.2f)\n', skew_power);
for k = 1:num_sines
    I_vec = I_vec + A_individual * sin(2*pi*freqs(k)*t_vec + phases(k));
end

%% (C) 엑셀 저장 (덮어쓰기 방지 포함)
export_tbl = table(t_vec, I_vec, 'VariableNames', {'Time', 'Current'});

% 파일명에 Skew 정보 추가하여 구분
save_xlsx_name = sprintf('MultiSine_Schroeder_Skew%.1f_dt%d.xlsx', skew_power, dt);
full_xlsx_path = fullfile(save_dir, save_xlsx_name);

if exist(full_xlsx_path, 'file')
    try
        delete(full_xlsx_path);
        fprintf('  >> [알림] 기존 파일 삭제 후 새로 생성합니다.\n');
    catch
        fprintf('  >> [경고] 파일 삭제 실패 (사용 중일 수 있음).\n');
    end
end

writetable(export_tbl, full_xlsx_path);
fprintf('  >> Data Saved: %s\n', save_xlsx_name);

%% (D) FFT 처리
fs = 1/dt;
N = numel(I_vec); 
win = hann(N,'periodic'); % Hanning Window
I_fft = fft((I_vec - mean(I_vec)) .* win);
halfIdx = 1:floor(N/2); 
f = (halfIdx-1)*fs/N;
mag = abs(I_fft(halfIdx))/(sum(win)/2);
phs_rad = angle(I_fft(halfIdx));

%% (E) 피크 탐색
[pks, locs] = findpeaks(mag, f, 'MinPeakProminence', 0.05*max(mag), 'SortStr','descend');
nShow = min(20, numel(pks));
if nShow > 0, idx=1:nShow; locs_col=locs(idx); pks_col=pks(idx); else, locs_col=[]; pks_col=[]; end

%% (F) 시각화 및 저장
fig_title = sprintf('Multi-Sine (Skew Power: %.1f)', skew_power);
f_handle = figure('Name', fig_title, 'Position',[100 50 800 900]);

% 1. Time Domain
subplot(3,1,1);
plot(t_vec, I_vec, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.0); 
title(sprintf('Time Domain: Max Amp %.2f A', max(abs(I_vec)))); 
xlabel('Time (s)'); ylabel('Current (A)');
grid on; xlim([0 t_end]);

CF = max(abs(I_vec)) / rms(I_vec);
text(0.02*t_end, max(I_vec)*0.8, sprintf('Crest Factor: %.2f', CF), 'FontSize', 10, 'BackgroundColor', 'w');

% 2. Magnitude (여기서 주파수 쏠림 확인 가능)
subplot(3,1,2);
plot(f, mag, 'k'); hold on;
if nShow>0, stem(locs_col, pks_col, 'r', 'filled'); end
% 주입된 주파수 위치 표시 (파란 점선)
for k=1:num_sines, xline(freqs(k), ':b', 'Alpha', 0.4); end
title(sprintf('FFT Magnitude (Skew: %.1f -> Freq Distribution Check)', skew_power)); 
grid on; xlim([0 f_max*1.5]);

% 3. Phase
subplot(3,1,3);
plot(f, phs_rad, '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 2); hold on;
inj_indices = zeros(1, num_sines);
for k = 1:num_sines, [~, m] = min(abs(f - freqs(k))); inj_indices(k) = m; end
plot(f(inj_indices), phs_rad(inj_indices), 'ro', 'LineWidth', 1.5, 'MarkerSize', 6);

theory_phase = zeros(1, num_sines);
for k=1:num_sines
    theory_phase(k) = angle(exp(1i * (- (k * (k - 1) / num_sines) * pi)));
end
plot(freqs, unwrap(theory_phase), 'b--', 'LineWidth', 1.0); 
title('FFT Phase'); 
grid on; xlim([0 f_max*1.5]); ylim([-pi pi]);
yticks([-pi 0 pi]); yticklabels({'-\pi', '0', '\pi'});

% Figure 저장
save_fig_name = sprintf('MultiSine_Skew%.1f_FFT.png', skew_power);
saveas(f_handle, fullfile(save_dir, save_fig_name));
fprintf('  >> Figure Saved: %s\n', save_fig_name);
fprintf('\n======= Skewed Sine Wave 생성 완료 =======\n');