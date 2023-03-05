# About

This repository contains code related to the following paper on LoneSTAR, an analog beamforming codebook design methodology for full-duplex millimeter wave (mmWave) wireless communication systems.

[1] I. P. Roberts, S. Vishwanath, and J. G. Andrews, "LoneSTAR: Analog Beamforming Codebooks for Full-Duplex Millimeter Wave Systems," _IEEE Trans. Wireless Commun._, 2023, [PDF](https://ianproberts.com/pdf/pub/lonestar.pdf).

Using the code in this repo, which is based on the work presented in [1], users can create LoneSTAR beamforming codebooks in order to conduct research on full-duplex mmWave communication systems.

If you use this code or our paper in your work, please cite [1] with the following BibTeX.

```
@ARTICLE{roberts_lonestar_2023,
    author={I. P. Roberts and S. Vishwanath and J. G. Andrews},
    title={\textsc{LoneSTAR}: Analog Beamforming Codebooks for Full-Duplex Millimeter Wave Systems},
    journal={IEEE Trans. Wireless Commun.},
    year=2023,
    note={(early access)},
}
```



# Codebook-Based Beam Alignment

Today's mmWave systems---which have operated in a half-duplex fashion thus far---rely on beam alignment to identify beams that deliver high beamforming gain when serving a particular user. 
This is typically done by sweeping a set of candidate beams (called a codebook) and measuring the reference signal received power (RSRP) for each beam candidate in order to select which is used to close a link. 
This procedure allows a mmWave system to reliably deliver high beamforming gain without the need for downlink/uplink MIMO channel knowledge, which is expensive to obtain in practice.

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/222982499-dd6d1b01-1d13-4f55-8e89-1c1d631ab9b1.svg"/>
</p>

# What is Self-Interference? 

When a transceiver (a wireless device) attempts to transmit and receive at the same time using the same frequency spectrum, some of its transmitted signal will leak into its receiver, corrupting reception of a desired signal.
This undesired leakage (or coupling) is called "self-interference". 
Transmit and receiving at the same time using the same frequency spectrum is called "full-duplex" operation.

This work is particularly focused on self-interference in full-duplex systems operating at mmWave frequencies (roughly 30 GHz to 100 GHz).
In mmWave systems, dense antenna arrays containing dozens or even hundreds of individual antennas are used to overcome high path loss at these high carrier frequencies. 
They do this by forming very directional beams, focusing energy in a particular direction to increase received signal power. 

Suppose a mmWave system employs separate transmit and receive arrays and attempts to operate in a full-duplex fashion.
Codebook-based beam alignment will presumably be conducted on its downlink and its uplink.
In other words, it will select a transmit beam from some transmit codebook and a receive beam from some receive codebook.
The degree of self-interference coupled depends on the selected beams and the underlying self-interference channel manifesting between the transmit and receive arrays.

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/222982526-5be5c14f-915b-48c6-a14f-156515d7f816.svg"/>
</p>

# What is LoneSTAR?

In this work [2], we construct the first beam codebooks for full-duplex mmWave systems, called LoneSTAR codebooks.
LoneSTAR codebooks deliver high beamforming gain and broad coverage while significantly reducing the self-interference coupled between the transmit and receive beams of a full-duplex mmWave system.
An example of a transmit codebook designed with LoneSTAR is shown below.

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/222982332-bd648e5d-db44-4ca1-b207-47f34ae550ef.svg"/>
</p>

To design codebooks for full-duplex mmWave systems, we formulated an optimization problem that aims to minimize self-interference between all possible beams within the codebooks while constraining the coverage quality the codebooks provide over some pre-defined coverage regions.

A full-duplex mmWave transceiver can independently select any transmit beam and any receive beam from LoneSTAR codebooks and, with reasonable confidence, can expect low self-interference. 
LoneSTAR codebooks are designed at the full-duplex transceiver following estimation of the self-interference MIMO channel (which need not be perfect), do not require downlink/uplink MIMO channel knowledge, and do not demand any over-the-air feedback to/from users.

After their construction, LoneSTAR codebooks can be used to conduct beam alignment and serve any downlink-uplink user pair in a full-duplex fashion thereafter with reduced self-interference. 
Importantly, LoneSTAR accounts for limited phase and amplitude control present in practical analog beamforming networks.

# Contents

This repo contains the following MATLAB code:
 - a main script `main.m` illustrating example usage
 - an `array.m` object that can be used to construct and interface with antenna arrays
 - a function `quantize_codebook.m` that quantizes a codebook to some phase and amplitude resolution
 - a function `generate_spherical_wave_channel.m` that produces a realization of the self-interference channel matrix based on some relative transmit-receive array geometry
 - `channel.m` and `channel_spherical_wave.m` objects that are used by `generate_spherical_wave_channel.m`

The `array.m`, `channel.m`, and `channel_spherical_wave.m` objects are all from [MIMO for MATLAB](https://mimoformatlab.com). 

# Example Usage

We will now walk through `main.m` to summarize its usage.

Suppose a full-duplex mmWave base station employs a codebook of transmit beams that it uses to conduct downlink beam alignment and a codebook of receive beams that it uses to conduct uplink beam alignment.
This script will design these transmit and receive codebooks based on LoneSTAR.

When running LoneSTAR, there are five design parameters:
- transmit coverage region
- receive coverage region
- transmit coverage variance
- receive coverage variance
- self-interference channel estimation error variance

### Construct Transmit and Receive Arrays

The transmit and receive arrays are assumed to be identical in this example, but this can be generalized straightforwardly. 
As shown below, the arrays can be created by simply specifying the number of rows and columns in the planar arrays.
Below, we have assumed 8x8 uniform planar arrays for the transmitter and receiver.

```
% planar array size
M = 8; % number of rows
N = 8; % number of columns

% create transmit and receive arrays (custom array object)
atx = array.create(M,N); % transmit array
arx = array.create(M,N); % receive array

% number of transmit and receive antennas
Nt = M * N;
Nr = M * N;
```

### Phase and Amplitude Control

Practically, analog beamforming networks use digitally-controlled phase shifters and attenuators to form beams.
Given this, each phase shifter and attenuator has limited resolution, which must be accounted for. 
The phase shifter resolution and attenuator resolution (in bits) can be specified using the following.

```
bits_phase = 6; % phase shifter resolution (Inf for infinite resolution)
bits_amp = 6; % attenuator resolution (Inf for infinite resolution)
dB_step = 0.5; % attenuation (in dB) per LSB (log-stepped attenuators)
```

Along with its resolution, the attenuator's step size is specified using `dB_step`, which is the attenuation (in dB) per least-significant bit.
The phase shifters have uniform phase quantization.

### Define Transmit and Receive Coverage Regions

Before running LoneSTAR, the transmit and receive codebooks used for conventional beam alignment at the full-duplex base station must be defined.
This can be done by defining the steering directions of the codebooks' beams. 
The steering direction of each beam contains a component in azimuth and elevation, as illustrated below.

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/222922576-377d2020-43aa-4f13-9db8-b9ada9246846.svg"/>
</p>

In the example below, the transmit and receive coverage regions are identical. Each coverage region spans in azimuth from -56 deg. to 56 deg. in 8 deg. steps and in elevation from -8 deg. to 8 deg. in 8 deg. steps. This amounts to a total of 45 directions in each coverage region. The codebooks we will design will steer in these specified directions, meaning each codebook will have 45 beams.

```
% coverage regions span azimuth and elevation
az_deg = [-56:8:56].';
el_deg = [-8:8:8].';

% assume TX and RX coverage regions to be the same
tx_az = az_deg;
tx_el = el_deg;
rx_az = az_deg;
rx_el = el_deg;
```

The azimuth-elevation of each transmit and receive coverage direction can then be populated as follows.

```
% full TX coverage region in az-el
num_tx_az = length(tx_az);
num_tx_el = length(tx_el);
tx_azel = [repmat(tx_az,num_tx_el,1) repelem(tx_el,num_tx_az,1)];

% full RX coverage region in az-el
num_rx_az = length(rx_az);
num_rx_el = length(rx_el);
rx_azel = [repmat(rx_az,num_rx_el,1) repelem(rx_el,num_rx_az,1)];
```

### Baseline Codebooks: Conjugate Beamforming Codebooks

A simple way to serve the coverage regions we have defined is using conjugate beamforming (CBF) codeboks. 
These can be formed as follows by simply fetching the array response vectors in each direction of our coverage regions.

```
% get array response vectors
Atx = atx.get_array_response(tx_azel(:,1)*pi/180,tx_azel(:,2)*pi/180);
Arx = arx.get_array_response(rx_azel(:,1)*pi/180,rx_azel(:,2)*pi/180);

% conjugate beamforming codebooks
F_cbf = Atx;
W_cbf = Arx;
```

We must quantize the beamforming weights in each codebook to ensure they are physically realizable based on our phase shifter and attenuator resolution.
This is done using the following.

```
% quantize phase and amplitude of beams
F_cbf = quantize_codebook(F_cbf,bits_phase,bits_amp,dB_step);
W_cbf = quantize_codebook(W_cbf,bits_phase,bits_amp,dB_step);
```

The azimuth cut of these beams can be visualized using the following, which produces the figure below.

```
idx_azimuth = find(tx_azel(:,2) == 0);
atx.show_beam_codebook_azimuth(F_cbf(:,idx_azimuth));
rlim([0,Nt]);
```

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/222983808-61270abe-05fb-4560-9a60-e2dc0565e449.svg"/>
</p>

### Setting the Size and Resolution of the Transmit and Receive Spatial Neighborhoods

Recall, STEER will search across some spatial neighborhoods to select transmit and receive beams that offer low self-interference. 
The size and resolution of the spatial neighborhoods can be defined respectively as follows.

```
% neighborhood size
Delta_az = 2;
Delta_el = 2;

% neighborhood resolution
delta_az = 1;
delta_el = 1;
```

The transmit and receive spatial neighborhoods have been assumed the same in this example but this can be generalized as desired.

### Constructing Spatial Neighborhoods

The transmit and receive spatial neighborhoods can be independently constructed with the following.

```
% TX and RX neighborhoods
nbr_tx = construct_neighborhood(Delta_az_tx,Delta_el_tx,delta_az_tx,delta_el_tx);
nbr_rx = construct_neighborhood(Delta_az_rx,Delta_el_rx,delta_az_rx,delta_el_rx);
```

The entire joint TX-RX neighborhood is then constructed as follows.

```
% full TX-RX neighborhood
nbr = [repmat(nbr_tx,length(nbr_rx(:,1)),1) repelem(nbr_rx,length(nbr_tx(:,1)),1)];
```

Finally, we sort the neighborhood based on some distance metric. Here we use the sum squared distance of azimuth and elevation of the transmit and receive neighborhoods. Note that this is slightly different than the distance metric used in [1]; other distance metrics can be used as desired.

```
% sort full neighborhood by distance (can modify distance metric as desired)
[~,idx] = sort(sum(nbr.^2,2));
nbr = nbr(idx,:);
```

### Executing STEER

Before running STEER, a self-interference target should be specified. For example, we used an INR target of -7 dB below (and in [1]).
Note that lower thresholds will require longer execution time and more measurement overhead in practice.

```
% set INR target (design parameter)
INR_tgt_dB = -7; % lower = stricter threshold
```

STEER is then executed by the following chunk of code. For each possible initial transmit-receive beam selection based on the defined codebooks, STEER searches for the beam pair within the surrounding spatial neighborhood that meets the desired self-interference threshold while minimally deviating from the initial selection. If this threshold cannot be found, STEER will default to the beam pair that minimized self-interference. 

```
% reset counter
idx_row = 0;

% reset lookup table
lut = [];

% for each TX beam in codebook
for idx_tx = 1:length(txcb_azel)
    % initial TX steering direction
    tx_azel = txcb_azel(idx_tx,:);
    
    % for each RX beam in codebook
    for idx_rx = 1:length(rxcb_azel)
        % initial RX steering direction
        rx_azel = rxcb_azel(idx_rx,:);
        
        % reset min INR
        INR_min_dB = Inf;
        
        % reset counter
        num_meas = 0;
                
        % search over neighborhood
        for idx_nbr = 1:num_nbr
            % shift TX direction
            tx_azel_shift = nbr(idx_nbr,1:2);
            tx_azel_nbr = tx_azel + tx_azel_shift;
            
            % shift RX direction
            rx_azel_shift = nbr(idx_nbr,3:4);
            rx_azel_nbr = rx_azel + rx_azel_shift;
            
            % measure INR (this depends on either a model or measurements)
            INR_meas_dB = 3*randn(1);
            
            % record number of measurements required
            num_meas = num_meas + 1;
            
            % record nominal INR (with initial selection)
            if idx_nbr == 1
                INR_nom_dB = INR_meas_dB;
            end
            
            % check if new minimum INR found
            if INR_meas_dB < INR_min_dB
                % record new best INR and steering directions
                INR_min_dB = INR_meas_dB;
                tx_azel_opt = tx_azel_nbr;
                rx_azel_opt = rx_azel_nbr;
                
                % check if target met
                if INR_meas_dB <= INR_tgt_dB
                    break;
                end
            end
        end
        
        % add to lookup table
        idx_row = idx_row + 1;
        lut(idx_row,:) = [tx_azel, tx_azel_opt, rx_azel, rx_azel_opt, INR_nom_dB, INR_min_dB, num_meas];
    end
end
```

Note that the following lines should be replaced with measurements or with a model that you have chosen.

```
% measure INR (this depends on either a model or measurements)
INR_meas_dB = 3*randn(1);
```

Our previous work [2] or [3] could be used to statistically model `INR_meas_dB`. Or you may want to compute `INR_meas_dB` based on some self-interference channel model. If you have collected measurements, you would need to fetch the measured INR corresponding to the particular transmit and receive beams of interest. If you have any questions about this, please reach out to the corresponding author of [1].

### Output: Lookup Table

The results of STEER are saved in the lookup table `lut`, where each row is a different initial transmit-receive beam selection. Its columns contain the following:
- cols 1-2: initial TX direction in az-el (degrees)
- cols 3-4: shifted TX direction in az-el (degrees)
- cols 5-6: initial RX direction in az-el (degrees)
- cols 7-8: shifted RX direction in az-el (degrees)
- col 9: INR (in dB) with initial TX-RX directions
- col 10: INR (in dB) with shifted TX-RX directions
- col 11: number of measurements taken

An example row from `lut` is shown below. 

| Initial TX Az. | Initial TX El. | Shifted TX Az. | Shifted TX El. | Initial RX Az. | Initial RX El. | Shifted RX Az. | Shifted RX El. | Initial INR | Resulting INR | Number of Measurements |
| :----: | :----: |  :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: |
| 16 | 8 | 17 | 6 | 32 | 0 | 31 | -2 | 7.1 | -10.9 | 435 |

The initial transmit direction was (16,8). The shifted transmit direction is (17,6). The initial receive direction was (32,0). The shifted receive direction is (31,-2). A shift of (+1,-2) was applied to the transmit beam. A shift of (-1,-2) was applied to the receive beam. These shifts resulted in an INR reduction from 7.1 dB to -10.9 dB. This was the first beam pair found that offered an INR below `INR_tgt_dB` and required 435 measurements to find.

### Output: Self-Interference Distribution

If we plot the CDFs of columns 9 and 10, we can observe the reduction in self-interference offered by STEER across all possible initial transmit-receive beam selections.

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/222924224-ab2d5564-72c6-4db5-a523-d1a41141511f.svg"/>
</p>

From this plot, we can see that a median user without STEER has an INR of about 0 dB, whereas with STEER, a median user has an INR of around -7 dB.


### Output: Number of Measurements

If we plot the CDF of column 11, we can see how many measurements are required when executing STEER across all possible initial transmit-receive beam selections.

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/222924306-d203a618-8aed-4e82-bf5f-3ef3776ab8fd.svg"/>
</p>

From this plot, we can see that 50% of the time, STEER requires 65 measurements. Exhaustively searching the entire neighborhood in this example requires 625 measurements (25 candiates in transmit neighborhood times 25 candidates in receive neighborhood).

# Questions and Feedback

Feel free to reach out to the corresponding author of [1] with any questions or feedback.

# Our Related Work

If you're interested in LoneSTAR, you may also be interested in our other related work, listed below.

[2] I. P. Roberts, A. Chopra, T. Novlan, S. Vishwanath, and J. G. Andrews, "Beamformed Self-Interference Measurements at 28 GHz: Spatial Insights and Angular Spread," _IEEE Trans. Wireless Commun._, Nov. 2022, [PDF](https://ianproberts.com/pdf/pub/bfsi.pdf), [GitHub](https://ianproberts.com/bfsi).

[3] I. P. Roberts, A. Chopra, T. Novlan, S. Vishwanath, and J. G. Andrews, "Spatial and Statistical Modeling of Multi-Panel Millimeter Wave Self-Interference," Submitted to _IEEE J. Sel. Areas Commun._, 2023, [PDF](https://ianproberts.com/pdf/pub/simodel.pdf), [GitHub](https://ianproberts.com/simodel).

[4] I. P. Roberts, A. Chopra, T. Novlan, S. Vishwanath, and J. G. Andrews, "STEER: Beam Selection for Full-Duplex Millimeter Wave Communication Systems," _IEEE Trans. Commun._, Oct. 2022, [PDF](https://ianproberts.com/pdf/pub/steer.pdf), [GitHub](https://ianproberts.com/steer).

[2] and [3] invovles the measurement and modeling of mmWave self-interference with 28 GHz phased arrays. [4] leverages a phenomenon observed in [2] to construct a beamforming-based full-duplex solution called STEER.

This related work can be found at https://ianproberts.com.

# Acknowledgments

This work has been supported by the National Science Foundation Graduate Research Fellowship Program (Grant No. DGE-1610403). Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
