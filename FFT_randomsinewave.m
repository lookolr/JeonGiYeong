%% ======================================================================
%  Standalone Multi-Sine Generator (No File Reading)
%  목적: 주행부하 파일 없이, 20개의 Sine Wave를 직접 합성하여 생성
% ======================================================================
clear; clc; close all;

% [사용자 설정]
seed_number = 1;        % 시드 번호 (이걸 바꾸면 위상이 바뀜)
dt          = 1;        % 샘플링 간격 [s]
t_end       = 2000;     % 총 시간 [s] (최저 주파수 주기를 충분히 포함하도록 설정)
A_individual = 2;       % 개별 사인파 진폭 [A]
num_sines   = 20;       % 사인파 개수

% [주파수 설정]
f_min = 0.0003183;       % 최저 주파수 (tau=500s)
f_max = 0.0386;        % 최고 주파수 (tau=4.12s)

% 1. 저장 경로 설정
save_dir = sprintf('G:\\공유 드라이브\\Battery Software Group (2025)\\Internship\\25년_동계인턴\\전기영\\Figure\\4주차 FFT\\FFT_random sine wave\\Seed %d', seed_number);
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
fprintf('저장 경로: %s\n', save_dir);

% 2. 시간 벡터 생성 (파일에서 읽는 대신 직접 생성)
t_vec = (0:dt:t_end)'; 

% 3. 위상 및 주파수 생성
rng(seed_number); % 시드 고정
freqs = logspace(log10(f_min), log10(f_max), num_sines);
phases = 2 * pi * rand(1, num_sines); % Random Phase (0~2pi)

% 4. 합성 (Superposition)
I_vec = zeros(size(t_vec));
fprintf('>> Sine Wave 합성 중... (Seed: %d)\n', seed_number);

for k = 1:num_sines
    % I = A * sin(2*pi*f*t + phase)
    I_vec = I_vec + A_individual * sin(2*pi*freqs(k)*t_vec + phases(k));
end

%% (C) 엑셀 저장
export_tbl = table(t_vec, I_vec, 'VariableNames', {'Time', 'Current'});

% 파일명 저장 (이제 특정 주행모드 이름 대신 Generic한 이름을 사용)
save_xlsx_name = sprintf('MultiSine_Only_Seed%d_1s.xlsx', seed_number);
full_xlsx_path = fullfile(save_dir, save_xlsx_name);
writetable(export_tbl, full_xlsx_path);
fprintf('  >> Data Saved: %s\n', save_xlsx_name);

%% (D) FFT 처리
fs = 1/dt;
N = numel(I_vec); 
win = hann(N,'periodic');

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
fig_title = sprintf('Multi-Sine Generation (Seed: %d)', seed_number);
f_handle = figure('Name', fig_title, 'Position',[100 50 800 900]);

% 1. Time Domain
subplot(3,1,1);
plot(t_vec, I_vec, 'Color', [0.85 0.32 0.09], 'LineWidth', 1.0); 
title(sprintf('Time Domain (Seed %d): Max Amp %.2f A', seed_number, max(abs(I_vec)))); 
xlabel('Time (s)'); ylabel('Current (A)');
grid on; xlim([0 t_end]);

% 2. Magnitude
subplot(3,1,2);
plot(f, mag, 'k'); hold on;
if nShow>0, stem(locs_col, pks_col, 'r', 'filled'); end
for k=1:num_sines, xline(freqs(k), ':b', 'Alpha', 0.4); end
title('FFT Magnitude'); grid on; xlim([0 f_max*1.5]);

% 3. Phase
subplot(3,1,3);
plot(f, phs_rad, '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 2); hold on;
% 주입한 주파수 위치의 위상만 빨간 점으로 표시
inj_indices = zeros(1, num_sines);
for k = 1:num_sines, [~, m] = min(abs(f - freqs(k))); inj_indices(k) = m; end
plot(f(inj_indices), phs_rad(inj_indices), 'ro', 'LineWidth', 1.5, 'MarkerSize', 6);
title(sprintf('FFT Phase (Seed %d)', seed_number)); 
grid on; xlim([0 f_max*1.5]); ylim([-pi pi]);
yticks([-pi 0 pi]); yticklabels({'-\pi', '0', '\pi'});

% Figure 저장
save_fig_name = sprintf('MultiSine_Only_Seed%d_FFT.png', seed_number);
saveas(f_handle, fullfile(save_dir, save_fig_name));
fprintf('  >> Figure Saved: %s\n', save_fig_name);

fprintf('\n======= 생성 완료 =======\n');