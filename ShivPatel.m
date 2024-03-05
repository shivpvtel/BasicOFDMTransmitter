
% Shiv Patel
% CMPEN 462 Mini project - Basic OFDM Transmitter

% For Simplicity, I commented out the code for the graphs of the other 
% Modulations Schemes and only left 64QAM. 

% For Number 7, you can uncomment line 133 to display the buffer for 64QAM

%%%%%%%%%% NUMBER 1 %%%%%%%%%%%%%
phrase = 'WirelessCommunicationSystemsandSecurityShivPatel';
ascii_values = uint8(phrase);

num_subcarriers = 2048; 
symbols_per_subframe = 14; 
BPS64QAM = 6;
num_bits = num_subcarriers * symbols_per_subframe * BPS64QAM;
num_repeats = ceil(num_bits / (length(ascii_values) * 8));
repeated_ascii = repmat(ascii_values, 1, num_repeats);
binary_string = dec2bin(repeated_ascii, 8);
bit_stream = reshape(binary_string', 1, []) - '0';
bit_stream = bit_stream(1:num_bits);
disp(num_repeats);
disp(num_bits);

fs = 30.72e6; 
CP_first_duration = 5.2e-6; 
CP_remaining_duration = 4.7e-6; 

%%%%%%%%%% NUMBER 2/3 %%%%%%%%%%%%%

parallel_BPSK = cell(1, symbols_per_subframe);
parallel_QPSK = cell(1, symbols_per_subframe);
parallel_QAM16 = cell(1, symbols_per_subframe);
parallel_QAM64 = cell(1, symbols_per_subframe);


for j = 1:symbols_per_subframe
    idx_start = (j-1) * num_subcarriers + 1;
    idx_end = j*num_subcarriers;
    
    if idx_end > length(bit_stream)
        break; 
    end
    
    BPSK_modulation = pskmod(bit_stream(idx_start:idx_end), 2);
    QPSK_modulation = pskmod(bit_stream(idx_start:idx_end), 4);
    QAM16_modulation = qammod(bit_stream(idx_start:idx_end), 16, 'UnitAveragePower', true);
    QAM64_modulation = qammod(bit_stream(idx_start:idx_end), 64, 'UnitAveragePower', true);
    
    parallel_BPSK{j} = BPSK_modulation;
    parallel_QPSK{j} = QPSK_modulation;
    parallel_QAM16{j} = QAM16_modulation;
    parallel_QAM64{j} = QAM64_modulation;
end

OFDM_symbols_time_domain = zeros(num_subcarriers, symbols_per_subframe);
OFDM_symbols_time_domain2 = zeros(num_subcarriers, symbols_per_subframe);
OFDM_symbols_time_domain3 = zeros(num_subcarriers, symbols_per_subframe);
OFDM_symbols_time_domain4 = zeros(num_subcarriers, symbols_per_subframe);

%%%%%%%%%% NUMBER 5 %%%%%%%%%%%%%

for k = 1:symbols_per_subframe
    parallel_data = parallel_BPSK{k};
    OFDM_symbol_time_domain = ifft(parallel_data, num_subcarriers);
    OFDM_symbols_time_domain(:, k) = OFDM_symbol_time_domain;

    parallel_data2 = parallel_QPSK{k};
    OFDM_symbol_time_domain2 = ifft(parallel_data2, num_subcarriers);
    OFDM_symbols_time_domain2(:, k) = OFDM_symbol_time_domain2;

    parallel_data3 = parallel_QAM16{k};
    OFDM_symbol_time_domain3 = ifft(parallel_data3, num_subcarriers);
    OFDM_symbols_time_domain3(:, k) = OFDM_symbol_time_domain3;

    parallel_data4 = parallel_QAM64{k};
    OFDM_symbol_time_domain4 = ifft(parallel_data4, num_subcarriers);
    OFDM_symbols_time_domain4(:, k) = OFDM_symbol_time_domain4;
end

%%%%%%%%%% NUMBER 6 %%%%%%%%%%%%%

CP_first_samples = round(CP_first_duration * fs);
CP_remaining_samples = round(CP_remaining_duration * fs);

num_OFDM_symbols = size(OFDM_symbols_time_domain, 2);
num_OFDM_symbols2 = size(OFDM_symbols_time_domain2, 2);
num_OFDM_symbols3 = size(OFDM_symbols_time_domain3, 2);
num_OFDM_symbols4 = size(OFDM_symbols_time_domain4, 2);

OFDM_symbols_with_CP = zeros(num_subcarriers + CP_first_samples, num_OFDM_symbols); 
OFDM_symbols_with_CP2 = zeros(num_subcarriers + CP_first_samples, num_OFDM_symbols2); 
OFDM_symbols_with_CP3 = zeros(num_subcarriers + CP_first_samples, num_OFDM_symbols3); 
OFDM_symbols_with_CP4 = zeros(num_subcarriers + CP_first_samples, num_OFDM_symbols4); 

OFDM_symbols_with_CP(:, 1) = [OFDM_symbols_time_domain(end-CP_first_samples+1:end, 1); OFDM_symbols_time_domain(:, 1)];
OFDM_symbols_with_CP2(:, 1) = [OFDM_symbols_time_domain2(end-CP_first_samples+1:end, 1); OFDM_symbols_time_domain2(:, 1)];
OFDM_symbols_with_CP3(:, 1) = [OFDM_symbols_time_domain3(end-CP_first_samples+1:end, 1); OFDM_symbols_time_domain3(:, 1)];
OFDM_symbols_with_CP4(:, 1) = [OFDM_symbols_time_domain4(end-CP_first_samples+1:end, 1); OFDM_symbols_time_domain4(:, 1)];


for k = 2:num_OFDM_symbols
    OFDM_symbol_with_CP = [OFDM_symbols_time_domain(end-CP_remaining_samples+1:end, k); OFDM_symbols_time_domain(:, k)];
    OFDM_symbols_with_CP(1:length(OFDM_symbol_with_CP), k) = OFDM_symbol_with_CP;
end

for k = 2:num_OFDM_symbols2
    OFDM_symbol_with_CP2 = [OFDM_symbols_time_domain2(end-CP_remaining_samples+1:end, k); OFDM_symbols_time_domain2(:, k)];
    OFDM_symbols_with_CP2(1:length(OFDM_symbol_with_CP2), k) = OFDM_symbol_with_CP2;
end

for k = 2:num_OFDM_symbols3
    OFDM_symbol_with_CP3 = [OFDM_symbols_time_domain3(end-CP_remaining_samples+1:end, k); OFDM_symbols_time_domain3(:, k)];
    OFDM_symbols_with_CP3(1:length(OFDM_symbol_with_CP3), k) = OFDM_symbol_with_CP3;
end

for k = 2:num_OFDM_symbols4
    OFDM_symbol_with_CP4 = [OFDM_symbols_time_domain4(end-CP_remaining_samples+1:end, k); OFDM_symbols_time_domain4(:, k)];
    OFDM_symbols_with_CP4(1:length(OFDM_symbol_with_CP4), k) = OFDM_symbol_with_CP4;
end

%%%%%%%%%% NUMBER 7 %%%%%%%%%%%%%

output_buffer = OFDM_symbols_with_CP(:);
output_buffer2 = OFDM_symbols_with_CP2(:);
output_buffer3 = OFDM_symbols_with_CP3(:);
output_buffer4 = OFDM_symbols_with_CP4(:);
% disp(output_buffer);
% disp(output_buffer2);
% disp(output_buffer3);

%disp(output_buffer4);

%%%%%%%%%% NUMBER 8 %%%%%%%%%%%%%

% %BPSK
% figure;
% subplot(2, 1, 1); 
% plot(real(OFDM_symbols_time_domain(:, 1)));
% title('BPSK Modulation Scheme - Time Domain Signal Without CP');
% xlabel('Sample Index');
% ylabel('Amplitude');
% grid on;
% 
% subplot(2, 1, 2); 
% plot(real(OFDM_symbols_with_CP(:, 1)));
% title('BPSK Modulation Scheme - Time Domain Signal With CP');
% xlabel('Sample Index');
% ylabel('Amplitude');
% grid on;
% 
% 
% %QPSK
% figure;
% subplot(2, 1, 1); 
% plot(real(OFDM_symbols_time_domain2(:, 1)));
% title('QPSK Modulation Scheme - Time Domain Signal Without CP');
% xlabel('Sample Index');
% ylabel('Amplitude');
% grid on;
% 
% subplot(2, 1, 2); 
% plot(real(OFDM_symbols_with_CP2(:, 1)));
% title('QPSK Modulation Scheme - Time Domain Signal With CP');
% xlabel('Sample Index');
% ylabel('Amplitude');
% grid on;
% 
% %16QAM
% figure;
% subplot(2, 1, 1); 
% plot(real(OFDM_symbols_time_domain3(:, 1)));
% title('16QAM Modulation Scheme - Time Domain Signal Without CP');
% xlabel('Sample Index');
% ylabel('Amplitude');
% grid on;
% 
% subplot(2, 1, 2); 
% plot(real(OFDM_symbols_with_CP3(:, 1)));
% title('16QAM Modulation Scheme - Time Domain Signal With CP');
% xlabel('Sample Index');
% ylabel('Amplitude');
% grid on;
% 


%64QAM
figure;
subplot(2, 1, 1); 
plot(real(OFDM_symbols_time_domain4(:, 1)));
title('64QAM Modulation Scheme - Time Domain Signal Without CP');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

subplot(2, 1, 2); 
plot(real(OFDM_symbols_with_CP4(:, 1)));
title('64QAM Modulation Scheme - Time Domain Signal With CP');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

%%%%%%%%%% NUMBER 10 %%%%%%%%%%%%%
fc = 2.4e9; 
t = (0:length(output_buffer)-1)'/fs; 
carrier_signal = cos(2*pi*fc*t);
upconverted_signal = real(output_buffer) .* carrier_signal;
upconverted_signal2 = real(output_buffer2) .* carrier_signal;
upconverted_signal3 = real(output_buffer3) .* carrier_signal;
upconverted_signal4 = real(output_buffer4) .* carrier_signal;
% figure;
% plot(t, upconverted_signal);
% title('BPSK Up-Converted Signal');
% xlabel('Time (s)');
% ylabel('Amplitude');
% figure;
% plot(t, upconverted_signal2);
% title('QPSK Up-Converted Signal');
% xlabel('Time (s)');
% ylabel('Amplitude');
% figure;
% plot(t, upconverted_signal3);
% title('16QAM Up-Converted Signal');
% xlabel('Time (s)');
% ylabel('Amplitude');
figure;
plot(t, upconverted_signal4);
title('64QAM Up-Converted Signal');
xlabel('Time (s)');
ylabel('Amplitude');
