clear; clc; close all;

%% 1. 시뮬레이션 통합 설정 (최적화 버전)
% --------------------------------------------------------
% [1] 고정 파라미터
f_min_fixed  = 0.0002;  % [Hz] 
f_max_fixed  = 100;     % [Hz] 
n_sines      = 20;      % Sine wave 개수
dt           = 0.0025;  
t_end        = 5000;    

% [2] 변수 파라미터 설정 (총 99개 조합)
% (A) Mean (중심 주파수): 0.001 ~ 0.1 (로그 스케일 11개)
target_means = logspace(log10(0.001), log10(0.1), 11);

% (B) Sigma (퍼짐 정도): 0.1 ~ 1.0 (선형 스케일 9개)
target_sigmas = linspace(0.1, 1.0, 9);

% [3] 저장 경로 설정
save_dir = 'G:\공유 드라이브\Battery Software Group (2025)\Internship\25년_동계인턴\전기영\Figure\6주차 FFT\99_Sets_Optimization';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% 2. 데이터 생성 및 저장 (.mat 방식)
% --------------------------------------------------------
fprintf('============ 데이터 생성 시작 (총 %d 세트) ============\n', length(target_means) * length(target_sigmas));
fprintf(' >> 설정: dt=%.4fs (100Hz 대응), t_end=%ds\n', dt, t_end);

% 공통 벡터 미리 생성 (속도 향상)
t_vec = (0:dt:t_end)'; % 2,000,001 행 (약 200만개)
freq_vec = logspace(log10(f_min_fixed), log10(f_max_fixed), n_sines)';
x_log = log10(freq_vec);

total_count = 0;
tic; % 시간 측정 시작

% (1) Mean 반복
for m_idx = 1:length(target_means)
    mu_freq = target_means(m_idx); 
    mu_log  = log10(mu_freq);      
    
    % (2) Sigma 반복
    for s_idx = 1:length(target_sigmas)
        sigma_val = target_sigmas(s_idx); 
        
        total_count = total_count + 1;
        
        % ----------------------------------------------------
        % A. 가중치(Weights) 계산 (Pure Gaussian)
        % ----------------------------------------------------
        z = (x_log - mu_log) / sigma_val;
        weights = exp(-z.^2 / 2); 
        amplitudes = weights / max(weights); % 정규화
        
        % ----------------------------------------------------
        % B. 신호 합성 (Schroeder Phase)
        % ----------------------------------------------------
        k_idx = (1:n_sines)';
        phases = - (k_idx .* (k_idx - 1) * pi) / n_sines;
        
        I_multi = zeros(size(t_vec));
        
        % 20개 사인파 합성
        for i = 1:n_sines
            I_multi = I_multi + amplitudes(i) * cos(2*pi*freq_vec(i)*t_vec + phases(i));
        end
        
        % ----------------------------------------------------
        % C. .mat 파일 저장 (속도/용량 최적)
        % ----------------------------------------------------
        file_tag = sprintf('MultiSine_Mean%.4f_Sigma%.1f', mu_freq, sigma_val);
        full_path_mat = fullfile(save_dir, [file_tag '.mat']);
        
        % 필요한 변수만 골라서 저장
        save(full_path_mat, 't_vec', 'I_multi', 'freq_vec', 'amplitudes', 'mu_freq', 'sigma_val');
        
        % 10개마다 진행 상황 출력
        if mod(total_count, 10) == 0 || total_count == 99
            elapsed_time = toc;
            fprintf('[%d/99] 저장완료: Mean=%.4f, Sigma=%.1f (경과시간: %.1f초)\n', ...
                total_count, mu_freq, sigma_val, elapsed_time);
        end
        
    end % Sigma Loop
end % Mean Loop

fprintf('\n============ 모든 작업 완료 ============\n');
fprintf('파일 저장 경로: %s\n', save_dir);