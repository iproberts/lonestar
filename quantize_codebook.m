function G = quantize_codebook(F,bits_phase,bits_amplitude,dB_step)

if bits_phase == Inf && bits_amplitude == Inf
    G = F;
else
    % normalize to beamforming weight with maximum magnitude
    % (assumes high-resolution attenuator/VGA before beamforming network)
    mag = max(abs(F(:)));
    F = F ./ mag;

    % all possible phase settings
    Mp = 2^bits_phase;
    phases = (0:Mp-1) .* 2*pi/Mp - pi;

    % all possible amplitude settings
    Ma = 2^bits_amplitude;
    if Ma > 1
        amplitudes = (0:(Ma-1));
    else
        amplitudes = 0;
    end
    amplitudes = amplitudes * -1 .* dB_step;
    amplitudes = 10.^(amplitudes./20); % dB_step is a power attenuation

    % quantize each beamforming weight in codebook (could probably be more efficient)
    [M,N] = size(F);
    G = F;
    for n = 1:N % for each beam in codebook
        for m = 1:M % for each weight in each beam
            % get phase and magnitude
            f = F(m,n);
            ph = angle(f);
            am = abs(f);

            % project to closest phase
            idx = argmin(abs(exp(1j.*ph)-exp(1j.*phases)));
            ph = phases(idx);

            % project to closest amplitude
            idx = argmin(abs(am-amplitudes));
            am = amplitudes(idx);

            % quantized beamforming weight
            G(m,n) = am * exp(1j*ph);
        end
    end

    % renormalize
    G = G .* mag;
end

end