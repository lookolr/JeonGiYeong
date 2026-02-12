%% ======================================================================
%  Multi-Sine Superposition & FFT Analysis (Final + Excel Save)
%  목적: 주행 부하 + Random Phase Multi-Sine 합성
%  출력 1: Figure (Time / Mag / Phase) - 파일당 1개
%  출력 2: Excel 데이터 (Time-Current) - 1초 간격 리샘플링
% ======================================================================
clear; clc; close all;

% 0) 저장 경로 설정
save_dir = 'G:\공유 드라이브\Battery Software Group (2025)\Internship\25년_동계인턴\전기영\Figure\4주차 FFT\FFT_random Phase';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
fprintf('저장 경로: %s\n', save_dir);

% 0-1) 파일 목록
driving_files = {
    'G:\공유 드라이브\Battery Software Lab\Protocols\2025 현대 NE 고품셀 평가 실험\주행부하_scaledown\udds_0725.xlsx'
    'G:\공유 드라이브\Battery Software Lab\Protocols\2025 현대 NE 고품셀 평가 실험\주행부하_scaledown\us06_0725.xlsx'
    'G:\공유 드라이브\Battery Software Lab\Protocols\2025 현대 NE 고품셀 평가 실험\주행부하_scaledown\BSL_CITY1_0725.xlsx'
    'G:\공유 드라이브\Battery Software Lab\Protocols\2025 현대 NE 고품셀 평가 실험\주행부하_scaledown\BSL_HW1_0725.xlsx'
};

% 파라미터
f_min = 0.001558; f_max = 0.03186; num_sines = 20; A_individual = 2;

for fileIdx = 1:numel(driving_files)
    % (A) 데이터 읽기
    filename = driving_files{fileIdx};
    [~, fname_only, ~] = fileparts(filename);
    data = readtable(filename);
    t_vec = data{:, 1}; I_vec = data{:, 2};
    mask = ~(isnan(t_vec) | isnan(I_vec));
    t_vec = t_vec(mask); I_vec = I_vec(mask);
    I_original = I_vec;

    % (B) 합성 (Random Phase)
    freqs = logspace(log10(f_min), log10(f_max), num_sines);
    phases = 2 * pi * rand(1, num_sines); % Random Phase
    I_noise = zeros(size(t_vec));
    for k = 1:num_sines
        I_noise = I_noise + A_individual * sin(2*pi*freqs(k)*t_vec + phases(k));
    end
    I_vec = I_original + I_noise; % 최종 합성 전류

    % (C) [복구됨] 데이터 리샘플링 (1s 간격) 및 엑셀 저장 -----------
    t_1s = (0 : 1 : floor(max(t_vec)))'; 
    I_1s = interp1(t_vec, I_vec, t_1s, 'linear'); % 선형 보간
    export_tbl = table(t_1s, I_1s, 'VariableNames', {'Time', 'Current'});
    
    save_xlsx_name = sprintf('%s_MultiSine_Random_1s.xlsx', fname_only);
    full_xlsx_path = fullfile(save_dir, save_xlsx_name);
    writetable(export_tbl, full_xlsx_path);
    fprintf('  >> Data Saved: %s\n', save_xlsx_name);
    % -------------------------------------------------------------

    % (D) FFT 처리
    dt = median(diff(t_vec)); if dt==0||isnan(dt), dt=1; end; fs = 1/dt;
    N = numel(I_vec); win = hann(N,'periodic');
    I_fft = fft((I_vec - mean(I_vec)) .* win);
    halfIdx = 1:floor(N/2); f = (halfIdx-1)*fs/N;
    mag = abs(I_fft(halfIdx))/(sum(win)/2);
    phs_rad = angle(I_fft(halfIdx));

    % (E) 피크 탐색
    [pks, locs] = findpeaks(mag, f, 'MinPeakProminence', 0.05*max(mag), 'SortStr','descend');
    peak_num = 10; nShow = min(peak_num, numel(pks));
    if nShow > 0
        idx=1:nShow; locs_col=locs(idx); pks_col=pks(idx);
    else
        locs_col=[]; pks_col=[];
    end

    % (F) 시각화 및 저장
    fig_title = sprintf('[%d] %s (Final + Data)', fileIdx, fname_only);
    f_handle = figure('Name', fig_title, 'Position',[100 50 800 900]);

    % 1. Time (전체 시간 출력)
    subplot(3,1,1);
    plot(t_vec, I_vec, 'Color', [0.85 0.32 0.09], 'LineWidth', 0.8); hold on;
    plot(t_vec, I_original, 'Color', [0 0.44 0.74], 'LineWidth', 1.2);
    legend('Superimposed', 'Original'); title('Time Domain'); grid on; 
    xlim([0 max(t_vec)]); % 전체 시간

    % 2. Mag
    subplot(3,1,2);
    plot(f, mag, 'k'); hold on;
    if nShow>0, stem(locs_col, pks_col, 'r', 'filled'); end
    for k=1:num_sines, xline(freqs(k), ':b', 'Alpha', 0.4); end
    title('FFT Magnitude'); grid on; xlim([0 f_max*1.5]);

    % 3. Phase
    subplot(3,1,3);
    plot(f, phs_rad, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 2); hold on;
    inj_indices = zeros(1, num_sines);
    for k = 1:num_sines, [~, m] = min(abs(f - freqs(k))); inj_indices(k) = m; end
    plot(f(inj_indices), phs_rad(inj_indices), 'ro', 'LineWidth', 1.5, 'MarkerSize', 5);
    title('FFT Phase'); grid on; xlim([0 f_max*1.5]); ylim([-pi pi]);

    % Figure 저장
    save_fig_name = sprintf('%s_Final_PhaseCheck.png', fname_only);
    saveas(f_handle, fullfile(save_dir, save_fig_name));
    fprintf('  >> Figure Saved: %s\n', save_fig_name);
end
fprintf('모든 분석 및 저장 완료!\n');