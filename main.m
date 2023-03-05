clc; clearvars; close all;
cvx_setup;

%% -------------------------------------------------------------------------
% 1. Create array objects.
% -------------------------------------------------------------------------
% planar array size
M = 8; % number of rows
N = 8; % number of columns

% create transmit and receive arrays (custom array object)
atx = array.create(M,N); % transmit array
arx = array.create(M,N); % receive array

% number of transmit and receive antennas
Nt = M * N;
Nr = M * N;

%% -------------------------------------------------------------------------
% 2. Set phase and amplitude control resolution.
% -------------------------------------------------------------------------
bits_phase = 6; % phase shifter resolution (Inf for infinite resolution)
bits_amp = 6; % attenuator resolution (Inf for infinite resolution)
dB_step = 0.5; % attenuation (in dB) per LSB (log-stepped attenuators)

%% -------------------------------------------------------------------------
% 3. Define TX and RX coverage regions.
% -------------------------------------------------------------------------
% coverage regions span azimuth and elevation
az_deg = [-56:8:56].';
el_deg = [-8:8:8].';

% assume TX and RX coverage regions to be the same
tx_az = az_deg;
tx_el = el_deg;
rx_az = az_deg;
rx_el = el_deg;

% full TX coverage region in az-el
num_tx_az = length(tx_az);
num_tx_el = length(tx_el);
tx_azel = [repmat(tx_az,num_tx_el,1) repelem(tx_el,num_tx_az,1)];

% full RX coverage region in az-el
num_rx_az = length(rx_az);
num_rx_el = length(rx_el);
rx_azel = [repmat(rx_az,num_rx_el,1) repelem(rx_el,num_rx_az,1)];

% size of codebooks we will create (number of beams)
Mt = length(tx_azel(:,1));
Mr = length(rx_azel(:,1));
disp(['Number of TX beams: ' num2str(Mt)]);
disp(['Number of RX beams: ' num2str(Mr)]);

%% -------------------------------------------------------------------------
% 4. Baseline codebooks are conjugate beamforming (CBF) codebooks.
% -------------------------------------------------------------------------
% get array response vectors
Atx = atx.get_array_response(tx_azel(:,1)*pi/180,tx_azel(:,2)*pi/180);
Arx = arx.get_array_response(rx_azel(:,1)*pi/180,rx_azel(:,2)*pi/180);

% conjugate beamforming codebooks
F_cbf = Atx;
W_cbf = Arx;

% quantize phase and amplitude of beams
F_cbf = quantize_codebook(F_cbf,bits_phase,bits_amp,dB_step);
W_cbf = quantize_codebook(W_cbf,bits_phase,bits_amp,dB_step);

% plot the transmit beams steering in the azimuth plane
idx_azimuth = find(tx_azel(:,2) == 0);
atx.show_beam_codebook_azimuth(F_cbf(:,idx_azimuth));
rlim([0,Nt]);

%% -------------------------------------------------------------------------
% 5. Create channel matrix (the estimate of the channel).
% -------------------------------------------------------------------------
channel_type = 'spherical-wave';
if strcmpi(channel_type,'spherical-wave')
    % spherical-wave channel
    H = generate_spherical_wave_channel(atx,arx);
else
    % Rayleigh channel
    H = cgauss_rv(0,1,Nr,Nt);
end

%% -------------------------------------------------------------------------
% 6. Run LoneSTAR.
% -------------------------------------------------------------------------
% set coverage variance (sigma^2) (in dB)
sigma_sq_tx_dB = -22;
sigma_sq_rx_dB = -22;

% set channel estimation error variance (epsilon^2) (in dB)
eps_sq_dB = -Inf; % set to -Inf for no channel estimation error

% precompute constraints
sigma_sq_tx = 10.^(sigma_sq_tx_dB./10);
sigma_sq_rx = 10.^(sigma_sq_rx_dB./10);
rhs_constraint_f = (sigma_sq_tx * Nt^2 * Mt);
rhs_constraint_w = (sigma_sq_rx * Nr^2 * Mr);
eps_sq = 10.^(eps_sq_dB./10);

% initialize codebooks
F = F_cbf;
W = W_cbf;

% print objective before LoneSTAR
E_init = W' * H * F;
disp(['Before LoneSTAR: ' num2str(10*log10(norm(E_init,'fro')^2)) ' dB']);

WH = W' * H ./ sqrt(Nt*Nr);

if eps_sq_dB == -Inf
    cvx_begin quiet
        cvx_precision low
        variable F(Nt,Mt) complex
            minimize( sum_square_abs(vec(WH * F)) )
        subject to
            sum_square_abs(Nt - diag(Atx' * F)) <= rhs_constraint_f;
            abs(vec(F)) <= 1;
    cvx_end
else
    eps_sq_W = norm(W,'fro')^2 * eps_sq ./ (Nt*Nr);
    cvx_begin quiet
        cvx_precision low
        variable F(Nt,Mt) complex
            minimize( sum_square_abs(vec(WH * F)) + eps_sq_W * sum_square_abs(vec(F)) )
        subject to
            sum_square_abs(Nt - diag(Atx' * F)) <= rhs_constraint_f;
            abs(F(:)) <= 1;
    cvx_end
end

% make sure no NaN in F
if sum(sum(isnan(F)))
    input('F contains NaN');
end
    
% project TX codebook
F = quantize_codebook(F,bits_phase,bits_amp,dB_step);
    
% solve for RX codebook
HF = H * F ./ sqrt(Nt*Nr);
if eps_sq_dB == -Inf
    cvx_begin quiet
        cvx_precision low
        variable W(Nr,Mr) complex
            minimize( sum_square_abs(vec(HF' * W)) )
        subject to
            sum_square_abs(Nr - diag(Arx' * W)) <= rhs_constraint_w;
            abs(vec(W)) <= 1;
    cvx_end
else
    eps_sq_F = norm(F,'fro')^2 * eps_sq ./ (Nt*Nr);
    cvx_begin quiet
        cvx_precision low
        variable W(Nr,Mr) complex
            minimize( sum_square_abs(vec(HF' * W)) + eps_sq_F * sum_square_abs(vec(W)) )
        subject to
            sum_square_abs(Nr - diag(Arx' * W)) <= rhs_constraint_w;
            abs(vec(W)) <= 1;
    cvx_end
end
    
% make sure no NaN in W
if sum(sum(isnan(W)))
    error('W contains NaN');
end

% project RX codebook
W = quantize_codebook(W,bits_phase,bits_amp,dB_step);

% print objective of design
E_post = W' * H * F;
disp(['After LoneSTAR: ' num2str(10*log10(norm(E_post,'fro')^2)) ' dB']);

%% -------------------------------------------------------------------------
% 7. Plot the transmit beams steering in the azimuth plane.
% -------------------------------------------------------------------------
idx_azimuth = find(tx_azel(:,2) == 0);
atx.show_beam_codebook_azimuth(F(:,idx_azimuth));
rlim([0,Nt]);

%% -------------------------------------------------------------------------
% 8. Compare the coupling matrix with and without LoneSTAR.
% -------------------------------------------------------------------------
% for plotting
min_val = min([10*log10(abs(E_init(:)).^2); 10*log10(abs(E_post(:)).^2)]); 
max_val = max([10*log10(abs(E_init(:)).^2); 10*log10(abs(E_post(:)).^2)]); 

% compare coupling matrix with and without LoneSTAR
figure(3);
subplot(1,2,1);
imagesc(10*log10(abs(E_init).^2));
title('Without LoneSTAR');
xlabel('Transmit Beam Index');
ylabel('Receive Beam Index');
c = colorbar('EastOutside');
c.Label.Interpreter = 'latex';
c.Label.String = 'Self-Inteference Coupling Factor (dB)';
caxis([min_val,max_val]);
axis equal tight;
subplot(1,2,2);
imagesc(10*log10(abs(E_post).^2));
title('With LoneSTAR');
xlabel('Transmit Beam Index');
ylabel('Receive Beam Index');
c = colorbar('EastOutside');
c.Label.Interpreter = 'latex';
c.Label.String = 'Self-Inteference Coupling Factor (dB)';
caxis([min_val,max_val]);
axis equal tight;

%% -------------------------------------------------------------------------
% 9. Compare the CDFs of coupling with and without LoneSTAR.
% -------------------------------------------------------------------------
figure(4);
[f,x] = ecdf(10*log10(abs(E_init(:)).^2));
plot(x,f,'k--','DisplayName','Without LoneSTAR');
hold on;
[f,x] = ecdf(10*log10(abs(E_post(:)).^2));
plot(x,f,'k-','DisplayName','With LoneSTAR');
hold off;
xlabel('Self-Interference Coupling Factor (dB)');
ylabel('Cumulative Probability');
grid on;
grid minor;
legend('Location','Southeast','FontSize',10);
