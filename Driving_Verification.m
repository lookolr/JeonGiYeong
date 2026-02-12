%% ======================================================================
%  주행부하 위상 펼쳐보기 (Unwrapped Phase Check)
%  목적: 주행부하의 위상이 "매끄러운 곡선"인지 "랜덤한 꺾은선"인지 확인
% ======================================================================
clear; clc; close all;
% 1. 파일 하나 선택 (예: UDDS) - 경로 수정 필요
filepath = 'G:\공유 드라이브\Battery Software Lab\Protocols\2025 현대 NE 고품셀 평가 실험\주행부하_scaledown\BSL_HW1_0725.xlsx';
[~, fname, ~] = fileparts(filepath);

% 2. 데이터 로드
data = readtable(filepath);
t = data{:, 1}; 
I_load = data{:, 2};
mask = ~(isnan(t) | isnan(I_load)); t = t(mask); I_load = I_load(mask);

% 3. FFT 수행
dt = median(diff(t)); fs = 1/dt; N = numel(I_load);
I_fft = fft((I_load - mean(I_load)) .* hann(N)); % Hann 창 적용

% 주파수 절반만 취급
halfIdx = 1:floor(N/2);
f_axis = (halfIdx-1)*fs/N;
fft_complex = I_fft(halfIdx);

% 4. 위상 추출 및 Unwrapping (핵심)
% (1) Wrapped Phase: -pi ~ +pi 사이로 접힌 값 (기존에 보던 것)
phase_wrapped = angle(fft_complex); 

% (2) Unwrapped Phase: 2pi 경계를 넘으면 계속 이어붙여서 펴주는 함수
phase_unwrapped = unwrap(phase_wrapped);

% ======================================================================
% 5. 시각화 (비교)
% ======================================================================
figure('Name', 'Driving Cycle Unwrapped Phase', 'Position', [100 100 1000 500], 'Color', 'w');

% 왼쪽: 기존에 보던 것 (Wrapped)
subplot(1,2,1);
plot(f_axis, phase_wrapped, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 2);
title({['(1) Wrapped Phase (' fname ')'], '-pi ~ +pi 제한됨'}, 'Interpreter', 'none');
xlabel('Frequency (Hz)'); ylabel('Phase (rad)');
ylim([-pi pi]); grid on;
yticks([-pi 0 pi]); yticklabels({'-\pi', '0', '\pi'});

% 오른쪽: 쭉 펴본 것 (Unwrapped)
subplot(1,2,2);
plot(f_axis, phase_unwrapped, 'b-', 'LineWidth', 0.5);
title({['(2) Unwrapped Phase (' fname ')'], '제한을 풀고 연결한 모습'}, 'Interpreter', 'none');
xlabel('Frequency (Hz)'); ylabel('Cumulative Phase (rad)');
grid on;

% 설명 텍스트 추가
sgtitle('주행부하는 "규칙(공식)"이 없으므로 Unwrapped 그래프도 불규칙하게 나타납니다', 'FontSize', 12, 'FontWeight', 'bold');

fprintf('>> 분석 완료. 오른쪽 그래프가 매끄러운 곡선이 아니라 "불규칙한 선"인지 확인하세요.\n');