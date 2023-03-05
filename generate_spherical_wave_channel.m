function H = generate_spherical_wave_channel(atx,arx)

% carrier frequency of arrays
carrier_frequency_Hz = 30e9;

% number of antennas
Nt = atx.num_antennas;
Nr = arx.num_antennas;

% create copies of the transmit and receive arrays
atx_pos = copy_object(atx);
arx_pos = copy_object(arx);

% vertically stack transmit and receive arrays
distance_wavelengths = 10;
atx_pos.translate_array(0,0,distance_wavelengths); % transmit array above receive array

% show arrays
if true
    atx_pos.set_marker('bx');
    arx_pos.set_marker('rx');
    show_arrays_3d([atx_pos,arx_pos]);
    axis equal;
end

% create channel object
channel_object = channel.create('spherical-wave');
channel_object.set_arrays(atx_pos,arx_pos);
channel_object.set_carrier_frequency(carrier_frequency_Hz);

% ensure channel energy (squared Frobenius norm) is normalized to Nt*Nr
channel_object.set_force_channel_energy_normalization(true);
channel_object.set_normalized_channel_energy(Nt*Nr);

% channel realization (will be deterministic for this channel)
H = channel_object.channel_realization();

end