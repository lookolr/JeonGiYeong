%% ======================================================================
%  Standalone Multi-Sine Generator (Schroeder Phase)
%  [수정완료] 엑셀 파일 초기화(삭제 후 저장) 기능 추가 버전
% ======================================================================
clear; clc; close all;

% [사용자 설정]
seed_number = 1;        

% -------------------------------------------------------------
% [설정 Tip]
% 정밀한 시뮬레이션을 원하시면 dt = 0.1 추천
% 빠른 확인을 원하시면 dt = 1 사용 (현재 설정)
% -------------------------------------------------------------
dt          = 0.01;        % 샘플링 간격 [s]
t_end       = 8000;     % 총 시간 [s] (tau=500s 커버를 위해 넉넉하게)
A_individual = 2;       % 개별 사인파 진폭 [A]
num_sines   = 20;       % 사인파 개수

% [주파수 설정]
f_min = 0.00018;      % 최저 주파수 -> tau2 = 159.15s (또는 500s) 대응
f_max = 0.53;        % 최고 주파수 -> tau1 = 1.59s (또는 4.12s) 대응

% 1. 저장 경로 설정
save_dir = sprintf('G:\\공유 드라이브\\Battery Software Group (2025)\\Internship\\25년_동계인턴\\전기영\\Figure\\4주차 FFT\\FFT_Schroeder sine wave');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
fprintf('저장 경로: %s\n', save_dir);

% 2. 시간 벡터 생성
t_vec = (0:dt:t_end)'; 

% 3. 주파수 및 위상 생성 (Schroeder Phase 적용)
freqs = logspace(log10(f_min), log10(f_max), num_sines);

% Schroeder Phase 공식: phi_k = - (k * (k-1) / K) * pi
phases = zeros(1, num_sines);
for k = 1:num_sines
    phases(k) = - (k * (k - 1) / num_sines) * pi;
end

% 4. 합성 (Superposition)
I_vec = zeros(size(t_vec));
fprintf('>> Sine Wave 합성 중... (Schroeder Phase)\n');
for k = 1:num_sines
    I_vec = I_vec + A_individual * sin(2*pi*freqs(k)*t_vec + phases(k));
end

%% (C) 엑셀 저장 (문제 해결 파트)
export_tbl = table(t_vec, I_vec, 'VariableNames', {'Time', 'Current'});

% 파일명에 'Schroeder' 명시
save_xlsx_name = sprintf('MultiSine_Schroeder_1s.xlsx');
full_xlsx_path = fullfile(save_dir, save_xlsx_name);

% [핵심 수정] 기존 파일이 존재하면 삭제하여 '유령 데이터' 방지
if exist(full_xlsx_path, 'file')
    try
        delete(full_xlsx_path);
        fprintf('  >> [알림] 기존 파일 삭제 후 새로 생성합니다.\n');
    catch
        fprintf('  >> [경고] 엑셀 파일이 열려있어 삭제에 실패했을 수 있습니다. 파일을 닫아주세요.\n');
    end
end

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
fig_title = sprintf('Multi-Sine Generation (Schroeder Phase)');
f_handle = figure('Name', fig_title, 'Position',[100 50 800 900]);

% 1. Time Domain
subplot(3,1,1);
plot(t_vec, I_vec, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.0); 
title(sprintf('Time Domain (Schroeder): Max Amp %.2f A', max(abs(I_vec)))); 
xlabel('Time (s)'); ylabel('Current (A)');
grid on; xlim([0 t_end]);

% 텍스트로 파고율(Crest Factor) 표시
CF = max(abs(I_vec)) / rms(I_vec);
text(0.02*t_end, max(I_vec)*0.8, sprintf('Crest Factor: %.2f', CF), 'FontSize', 10, 'BackgroundColor', 'w');

% 2. Magnitude
subplot(3,1,2);
plot(f, mag, 'k'); hold on;
if nShow>0, stem(locs_col, pks_col, 'r', 'filled'); end
for k=1:num_sines, xline(freqs(k), ':b', 'Alpha', 0.4); end
title('FFT Magnitude'); grid on; xlim([0 f_max*1.5]);

% 3. Phase (Schroeder 곡선 확인)
subplot(3,1,3);
% 배경: 전체 FFT 위상
plot(f, phs_rad, '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 2); hold on;
% 주입한 주파수 위치의 위상 (빨간 점)
inj_indices = zeros(1, num_sines);
for k = 1:num_sines, [~, m] = min(abs(f - freqs(k))); inj_indices(k) = m; end
plot(f(inj_indices), phs_rad(inj_indices), 'ro', 'LineWidth', 1.5, 'MarkerSize', 6);

% Schroeder 이론 곡선 (파란 선) - 시각적 비교용
theory_phase = zeros(1, num_sines);
for k=1:num_sines
    theory_phase(k) = angle(exp(1i * (- (k * (k - 1) / num_sines) * pi)));
end
% 이론값은 점선으로 연결하여 경향성 확인
plot(freqs, unwrap(theory_phase), 'b--', 'LineWidth', 1.0); 
title('FFT Phase (Red: Measured, Blue: Schroeder Curve)'); 
grid on; xlim([0 f_max*1.5]); ylim([-pi pi]);
yticks([-pi 0 pi]); yticklabels({'-\pi', '0', '\pi'});

% Figure 저장
save_fig_name = sprintf('MultiSine_Schroeder_FFT.png');
saveas(f_handle, fullfile(save_dir, save_fig_name));
fprintf('  >> Figure Saved: %s\n', save_fig_name);
fprintf('\n======= Schroeder Phase 생성 완료 =======\n');