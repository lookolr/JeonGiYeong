%% ======================================================================
%  드라이빙 부하(Current) 파일  →  FFT(Hann 창)  →  진폭 & 위상(Phase) 분석
%  + Phase Spectrum (위상) 그래프 추가
% ======================================================================
clear; clc; close all;

% 0) 저장 경로 설정 --------------------------------------------------------
save_dir = 'G:\공유 드라이브\Battery Software Group (2025)\Internship\25년_동계인턴\전기영\Figure\4주차 FFT\Phase_Check';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
    fprintf('저장 폴더가 없어 새로 생성했습니다: %s\n', save_dir);
else
    fprintf('저장 경로 확인됨: %s\n', save_dir);
end

% 0-1) 분석할 파일 목록 ---------------------------------------------------
driving_files = {
    'G:\공유 드라이브\Battery Software Lab\Protocols\2025 현대 NE 고품셀 평가 실험\주행부하_scaledown\udds_0725.xlsx'
    'G:\공유 드라이브\Battery Software Lab\Protocols\2025 현대 NE 고품셀 평가 실험\주행부하_scaledown\us06_0725.xlsx'
    'G:\공유 드라이브\Battery Software Lab\Protocols\2025 현대 NE 고품셀 평가 실험\주행부하_scaledown\BSL_CITY1_0725.xlsx'
    'G:\공유 드라이브\Battery Software Lab\Protocols\2025 현대 NE 고품셀 평가 실험\주행부하_scaledown\BSL_HW1_0725.xlsx'
};

% Pulse (필요시 사용)
t_end = 300; dt = 0.1;
pulse.t = (0:dt:t_end)'; pulse.I = zeros(size(pulse.t));
pulse.I(pulse.t>=10 & pulse.t<=20) = 1;
% driving_files = [ driving_files; {pulse} ]; % 펄스 포함하려면 주석 해제

all_peak_tbl = struct; 

for fileIdx = 1:numel(driving_files)
    
    %% (A) 데이터 읽기 ---------------------------------------------------
    item = driving_files{fileIdx};
    
    if isstruct(item)
        t_vec = item.t(:); I_vec = item.I(:);
        file_label = 'Pulse'; 
    else
        filename = item;
        data = readtable(filename);
        try
            t_vec = data{:, 1}; I_vec = data{:, 2};
        catch
            error('파일 읽기 오류: %s', filename);
        end
        file_label = filename;  
    end
    
    % NaN 제거
    mask  = ~(isnan(t_vec) | isnan(I_vec));
    t_vec = t_vec(mask);
    I_vec = I_vec(mask);
    
    %% (B) 샘플링 정보 ---------------------------------------------------
    dt = median(diff(t_vec));
    if dt == 0 || isnan(dt), dt = 1; end
    fs = 1/dt;
    
    %% (C) FFT (진폭 & 위상 계산) ----------------------------------------
    N       = numel(I_vec);
    win     = hann(N,'periodic');
    
    % FFT 계산
    I_fft   = fft((I_vec - mean(I_vec)) .* win);
    
    halfIdx = 1:floor(N/2);
    f       = (halfIdx-1) * fs / N; 
    
    % 1. 진폭 (Magnitude)
    mag     = abs(I_fft(halfIdx)) / (sum(win)/2);
    
    % 2. 위상 (Phase) - [중요] 여기서 위상 추출
    % angle 함수는 라디안(-pi ~ pi) 값을 반환합니다.
    phs_rad = angle(I_fft(halfIdx)); 
    
    %% (D) 피크 탐색 -----------------------------------------------------
    [pks, locs] = findpeaks(mag, f, 'MinPeakProminence', 0.05*max(mag), 'SortStr','descend'); 
    peak_num = 10; nShow = min(peak_num, numel(pks)); 
    
    if nShow == 0
        peak_tbl = table();
        locs_col = []; pks_col = []; 
    else
        idx       = 1:nShow;
        locs_col  = locs(idx); locs_col = locs_col(:);      
        pks_col   = pks(idx);  pks_col  = pks_col(:);       
        tau_col   = 1./(2*pi*locs_col);
        magdb_col = 20*log10(pks_col);
        
        peak_tbl  = table(locs_col, tau_col, pks_col, magdb_col, ...
            'VariableNames', {'Freq_Hz','Tau_s','Amplitude','Mag_dB'});
    end
    all_peak_tbl.(sprintf('file%d',fileIdx)) = peak_tbl;
    
    %% (E) 시각화 (3단 구성: 시간 / 진폭 / 위상) ---------------------------
    [~, fname_only, ext] = fileparts(file_label);
    if isstruct(item), fname_only = 'Pulse'; ext=''; end
    
    fig_title = sprintf('[%d] %s Analysis', fileIdx, fname_only);
    
    % Figure 크기를 좀 더 세로로 길게 설정
    f_handle = figure('Name', fig_title, 'Position', [100 50 800 800]);
    
    % 1. 시간 영역
    subplot(3,1,1);
    plot(t_vec, I_vec, 'LineWidth', 1.0); grid on;
    xlabel('Time (s)'); ylabel('Current (A)');
    title(['(1) Time-domain: ' fname_only], 'Interpreter', 'none');
    
    % 2. 주파수 영역 (Magnitude)
    subplot(3,1,2);
    plot(f, mag, 'LineWidth', 1.0); hold on;
    if nShow > 0
        stem(locs_col, pks_col, 'r', 'filled', 'LineWidth', 1.3);
        legend('Spectrum', 'Detected Peaks');
    end
    grid on; xlabel('Frequency (Hz)'); ylabel('|I(f)| (Amp)');
    title('(2) FFT Magnitude Spectrum');
    
    % 3. 주파수 영역 (Phase) - [추가됨]
    subplot(3,1,3);
    % 위상은 점들이 많으므로 점(scatter) 형식이나 얇은 선으로 표현하는 게 잘 보임
    plot(f, phs_rad, '.', 'Color', [0.4 0.4 0.4], 'MarkerSize', 3); 
    grid on; 
    xlabel('Frequency (Hz)'); ylabel('Phase (rad)');
    ylim([-pi pi]); % 위상은 -3.14 ~ +3.14 사이
    yticks([-pi -pi/2 0 pi/2 pi]);
    yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
    title('(3) FFT Phase Spectrum (Randomness Check)');
    
    % === [저장] ===
    save_filename = sprintf('%s_PhaseCheck.png', fname_only);
    full_save_path = fullfile(save_dir, save_filename);
    saveas(f_handle, full_save_path); 
    fprintf('  >> Figure Saved: %s\n', save_filename);
    
end

fprintf('\n======= 분석 완료 =======\n');