% --- FULL WIRELESS SIMULATOR (FINAL WORKING VERSION) ---
clear; close all; clc;

figure('Position', [100 100 850 650]); hold on;

% 1) AWGN
cfg = get_default_config();
cfg.channel = 'awgn'; 
cfg.snr_end = 12;
[s1, b1] = run_simulation(cfg);
semilogy(s1, b1, 'o-', 'LineWidth', 1.4, 'DisplayName', 'AWGN');

% 2) Rayleigh, no diversity
cfg = get_default_config();
cfg.channel = 'rayleigh'; 
cfg.n_branches = 1;
[s2, b2] = run_simulation(cfg);
semilogy(s2, b2, 's-', 'LineWidth', 1.4, 'DisplayName', 'Rayleigh (no diversity)');

% 3) Rayleigh + MRC (Perfect CSI)
cfg = get_default_config();
cfg.n_branches = 2; 
cfg.combining = 'MRC'; 
cfg.estimation = 'perfect';
[s3, b3] = run_simulation(cfg);
semilogy(s3, b3, '^-', 'LineWidth', 1.4, 'DisplayName', 'Rayleigh + MRC (perfect CSI)');

% 4) Rayleigh + MRC (MMSE Estimation - NOVEL)
cfg = get_default_config();
cfg.n_branches = 2; 
cfg.combining = 'MRC'; 
cfg.estimation = 'mmse';
[s4, b4] = run_simulation(cfg);
semilogy(s4, b4, 'd-', 'LineWidth', 1.4, 'DisplayName', 'Rayleigh + MRC (MMSE est.)');

% 5) Rayleigh + MRC (Correlated Antennas - NOVEL)
cfg = get_default_config();
cfg.n_branches = 2; 
cfg.combining = 'MRC'; 
cfg.estimation = 'perfect';
cfg.rho = 0.8; 
[s5, b5] = run_simulation(cfg);
semilogy(s5, b5, 'v--', 'LineWidth', 1.4, 'DisplayName', 'Rayleigh + MRC (Correlated \rho=0.8)');

grid on; set(gca, 'YScale', 'log');
xlabel('SNR (dB)'); ylabel('BER');
title('BPSK: AWGN vs Rayleigh vs MRC (Novel Analysis)');
legend('Location', 'southwest');
ylim([1e-5 1]);

% Save the plot automatically to your folder
saveas(gcf, 'ber_curves.png');
fprintf('Simulation complete. Plot saved as ber_curves.png\n');


% ==========================================
% LOCAL FUNCTIONS
% ==========================================

function cfg = get_default_config()
    cfg.mod = 'BPSK';
    cfg.channel = 'rayleigh';
    cfg.combining = 'MRC';
    cfg.estimation = 'perfect';
    cfg.n_branches = 2;
    cfg.rho = 0.0;
    cfg.snr_start = 0;
    cfg.snr_step = 2;
    cfg.snr_end = 20;
    cfg.n_bits = 100000;
    cfg.seed = 42;
end

function [snrs, bers] = run_simulation(cfg)
    rng_stream = RandStream('mt19937ar', 'Seed', cfg.seed);
    snrs = cfg.snr_start : cfg.snr_step : cfg.snr_end;
    bers = zeros(size(snrs));

    for k = 1:numel(snrs)
        bers(k) = run_point(cfg, snrs(k), rng_stream);
        fprintf('  SNR = %3d dB   BER = %.5e\n', snrs(k), bers(k));
    end
end

function ber = run_point(cfg, snr_db, rng_stream)
    bits = randi(rng_stream, [0 1], cfg.n_bits, 1);
    symbols = 2*bits - 1; % BPSK
    N = numel(symbols);

    if strcmpi(cfg.channel, 'awgn')
        snr_lin = 10^(snr_db/10);
        noise = sqrt(1/(2*snr_lin)) * (randn(rng_stream, N, 1) + 1j*randn(rng_stream, N, 1));
        y_eq = symbols + noise;
    else
        snr_lin = 10^(snr_db/10);
        noise_var = 1 / snr_lin;
        
        % --- BLOCK FADING GENERATION ---
        block_size = 100; % Channel stays constant for 100 symbols
        num_blocks = ceil(N / block_size);
        
        if cfg.n_branches == 1
            h_block = (randn(rng_stream, 1, num_blocks) + 1j*randn(rng_stream, 1, num_blocks)) / sqrt(2);
        else
            R = toeplitz(cfg.rho .^ (0:cfg.n_branches-1));
            L = chol(R, 'lower');
            raw = (randn(rng_stream, cfg.n_branches, num_blocks) + 1j*randn(rng_stream, cfg.n_branches, num_blocks)) / sqrt(2);
            h_block = L * raw;
        end
        
        % Repeat the channel for every symbol in the block
        h = repelem(h_block, 1, block_size);
        h = h(:, 1:N); % Trim to exact length

        noise = sqrt(noise_var/2) * (randn(rng_stream, cfg.n_branches, N) + 1j*randn(rng_stream, cfg.n_branches, N));
        y = h .* symbols.' + noise;

        % Channel Estimation
        if strcmpi(cfg.estimation, 'perfect')
            h_hat = h;
        else
            h_hat = zeros(size(h));
            pilot_spacing = block_size; % 1 pilot per block
            pidx = 1:pilot_spacing:N;
            
            for i = 1:cfg.n_branches
                y_p = y(i, pidx).';
                x_p = symbols(pidx);
                h_ls = y_p ./ x_p;
                
                if strcmpi(cfg.estimation, 'mmse')
                    h_ls = h_ls * (1 / (1 + noise_var));
                end
                
                % Zero-order hold (repeat the estimate for the block)
                h_est_block = repelem(h_ls, pilot_spacing);
                h_hat(i, :) = h_est_block(1:N);
            end
        end

        % MRC Combining
        w = conj(h_hat);
        num = sum(w .* y, 1);
        den = sum(abs(h_hat).^2, 1);
        y_eq = num ./ (den + 1e-12);
    end

    % Demodulate
    bits_hat = double(real(y_eq(:)) >= 0);
    ber = mean(bits ~= bits_hat);
end