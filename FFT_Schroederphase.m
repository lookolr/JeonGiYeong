%% ======================================================================
%  드라이빙 부하 + Multi-Sine (Schroeder) → FFT & Phase Check
%  ● Freq Range: 0.001558 Hz ~ 0.03186 Hz (Log Scale 20개)
%  ● Amplitude : 개별 2A
%  ● Phase     : Schroeder Phase (Deterministic & Low Crest Factor)
%  ● Output    : Figure (PNG) & Data (Excel, 1s sampling)
%  ● Check     : Phase Spectrum 추가 (위상 분포 확인용)
% ======================================================================
clear; clc; close all;

% 0) 기본 저장 경로 설정 ---------------------------------------------------
base_dir = 'G:\공유 드라이브\Battery Software Group (2025)\Internship\25년_동계인턴\전기영\Figure\4주차 FFT\FFT_Schroeder_Phase_Check';
fig_dir = fullfile(base_dir, 'Figure');
data_dir = fullfile(base_dir, 'Synthesized_Data');

if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

fprintf('결과 저장 경로:\n - Figure: %s\n - Data  : %s\n', fig_dir, data_dir);
fprintf('적용된 위상: Schroeder Phase (Phase Spectrum 확인 가능)\n');

% 0-1) 분석할 파일 목록 ---------------------------------------------------
driving_files = {
    'G:\공유 드라이브\Battery Software Lab\Protocols\2025 현대 NE 고품셀 평가 실험\주행부하_scaledown\udds_0725.xlsx'
    'G:\공유 드라이브\Battery Software Lab\Protocols\2025 현대 NE 고품셀 평가 실험\주행부하_scaledown\us06_0725.xlsx'
    'G:\공유 드라이브\Battery Software Lab\Protocols\2025 현대 NE 고품셀 평가 실험\주행부하_scaledown\BSL_CITY1_0725.xlsx'
    'G:\공유 드라이브\Battery Software Lab\Protocols\2025 현대 NE 고품셀 평가 실험\주행부하_scaledown\BSL_HW1_0725.xlsx'
};

% ── Multi-Sine 파라미터 설정 ──────────────────────────────────────────
f_start = 0.001558;          % [Hz]
f_end   = 0.03186;           % [Hz]
num_sines = 20;              % 사인파 개수
A_each  = 2;                 % [A]

freqs_add = logspace(log10(f_start), log10(f_end), num_sines);
all_peak_tbl = struct;       

for fileIdx = 1:numel(driving_files)
    
    %% (A) 데이터 읽기 ---------------------------------------------------
    filename = driving_files{fileIdx};
    [~, fname_only, ~] = fileparts(filename);
    data = readtable(filename); 
    try
        t_vec = data{:, 1}; I_vec = data{:, 2};
    catch
        error('데이터 읽기 실패: %s', filename);
    end
    mask  = ~(isnan(t_vec) | isnan(I_vec));
    t_vec = t_vec(mask); I_vec = I_vec(mask);
    
    %% (B) Multi-Sine 합성 (Schroeder Phase 적용) -----------------------
    % 공식: phase[k] = -k * (k-1) * pi / N
    phases = zeros(1, num_sines);
    for k = 1:num_sines
        phases(k) = -k * (k - 1) * pi / num_sines;
    end
    
    I_noise = zeros(size(t_vec));
    for k = 1:num_sines
        f_k = freqs_add(k); 
        p_k = phases(k);
        % 개별 사인파 누적
        I_noise = I_noise + A_each * sin(2*pi*f_k * t_vec + p_k);
    end
    
    I_vec_added = I_vec + I_noise; 
    
    %% (C) 데이터 리샘플링 (1s 간격) 및 엑셀 저장 -------------------------
    t_1s = (0 : 1 : floor(max(t_vec)))'; 
    I_1s = interp1(t_vec, I_vec_added, t_1s, 'linear');
    export_tbl = table(t_1s, I_1s, 'VariableNames', {'Time', 'Current'});
    
    save_xlsx_name = sprintf('%s_MultiSine_Schroeder_1s.xlsx', fname_only);
    full_xlsx_path = fullfile(data_dir, save_xlsx_name);
    
    writetable(export_tbl, full_xlsx_path);
    fprintf('  >> Excel Saved: %s\n', save_xlsx_name);

    %% (D) FFT 및 위상(Phase) 계산 ---------------------------------------
    dt = median(diff(t_vec)); if dt==0||isnan(dt), dt=1; end; fs = 1/dt;
    N = numel(I_vec_added); win = hann(N,'periodic');
    
    % FFT 수행
    I_fft = fft((I_vec_added - mean(I_vec_added)) .* win);
    
    halfIdx = 1:floor(N/2); 
    f       = (halfIdx-1)*fs/N;
    
    % 1. 진폭 (Magnitude)
    mag     = abs(I_fft(halfIdx))/(sum(win)/2);
    
    % 2. 위상 (Phase) - [중요]
    phs_rad = angle(I_fft(halfIdx));
    
    %% (E) 피크 탐색 -----------------------------------------------------
    [pks, locs] = findpeaks(mag, f, 'MinPeakProminence', 0.01*max(mag), 'SortStr','descend');
    peak_num = 25; nShow = min(peak_num, numel(pks));
    if nShow == 0
        peak_tbl = table(); locs_col=[]; pks_col=[]; tau_col=[]; magdb_col=[];
    else
        idx=1:nShow; locs_col=locs(idx); pks_col=pks(idx);
        safe_locs=locs_col; safe_locs(safe_locs==0)=NaN;
        tau_col=1./(2*pi*safe_locs); magdb_col=20*log10(pks_col);
        peak_tbl=table(locs_col(:), tau_col(:), pks_col(:), magdb_col(:), ...
            'VariableNames', {'Freq_Hz','Tau_s','Amplitude','Mag_dB'});
    end
    all_peak_tbl.(sprintf('file%d',fileIdx)) = peak_tbl;
    
    %% (F) 시각화 (시간 / 진폭 / 위상) ------------------------------------
    fig_title = sprintf('[%d] %s (Schroeder Phase)', fileIdx, fname_only);
    f_handle = figure('Name', fig_title, 'Position',[100 50 800 900]); % 세로로 길게
       
    % (F-1) 시간 영역
    subplot(3,1,1);
    plot(t_vec, I_vec_added, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 0.8); grid on;
    xlabel('Time (s)'); ylabel('Current (A)');
    title(sprintf('Time-domain: %s (Schroeder)', fname_only), 'Interpreter', 'none');
    
    % (F-2) 주파수 영역 (Magnitude)
    subplot(3,1,2);
    plot(f, mag, 'Color', [0.2 0.2 0.2], 'LineWidth', 1.0); hold on;
    if nShow > 0, stem(locs_col, pks_col, 'r', 'filled', 'LineWidth', 1.2); end
    
    % 주입한 주파수 위치 표시
    for k = 1:num_sines, xline(freqs_add(k), ':b', 'Alpha', 0.4); end
    
    legend('Spectrum', 'Detected Peaks', 'Injected Freqs');
    grid on; xlabel('Frequency (Hz)'); ylabel('|I(f)|');
    title(sprintf('FFT Magnitude (Range: %.4f ~ %.4f Hz)', f_start, f_end));
    xlim([0 f_end * 1.5]); 
    
    % (F-3) 주파수 영역 (Phase) - [추가됨]
    subplot(3,1,3);
    
    % 배경 노이즈(주행부하)의 위상: 회색 점으로 표시
    plot(f, phs_rad, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 2); hold on;
    
    % 우리가 주입한 주파수(Schroeder Phase) 위치의 위상: 빨간 동그라미로 강조
    % (정확한 인덱스를 찾아서 표시)
    inj_indices = zeros(1, num_sines);
    for k = 1:num_sines
        [~, min_idx] = min(abs(f - freqs_add(k)));
        inj_indices(k) = min_idx;
    end
    plot(f(inj_indices), phs_rad(inj_indices), 'ro', 'LineWidth', 1.5, 'MarkerSize', 5);
    
    grid on; xlabel('Frequency (Hz)'); ylabel('Phase (rad)');
    ylim([-pi pi]); 
    yticks([-pi -pi/2 0 pi/2 pi]); yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
    title('FFT Phase Spectrum (Grey: Noise, Red: Schroeder Signal)');
    xlim([0 f_end * 1.5]); % 관심 영역 확대
    
    % === [저장] ===
    save_filename = sprintf('%s_MultiSine_Schroeder_FFT_Phase.png', fname_only);
    full_save_path = fullfile(fig_dir, save_filename);
    
    saveas(f_handle, full_save_path);
    fprintf('  >> Figure Saved: %s\n', save_filename);
end

fprintf('\n======= Multi-Sine (Schroeder) 분석 완료 =======\n');
fprintf('Data saved in: %s\n', data_dir);