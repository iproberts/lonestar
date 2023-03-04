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

# Our Related Work

The work of [1] was inpsired by phenomena observed when analyzing our nearly 6.5 million measurements of self-interference in the following papers.

[2] I. P. Roberts, A. Chopra, T. Novlan, S. Vishwanath, and J. G. Andrews, "Beamformed Self-Interference Measurements at 28 GHz: Spatial Insights and Angular Spread," _IEEE Trans. Wireless Commun._, Nov. 2022, [PDF](https://ianproberts.com/pdf/pub/bfsi.pdf), [GitHub](https://ianproberts.com/bfsi).

[3] I. P. Roberts, A. Chopra, T. Novlan, S. Vishwanath, and J. G. Andrews, "Spatial and Statistical Modeling of Multi-Panel Millimeter Wave Self-Interference," Submitted to _IEEE J. Sel. Areas Commun._, 2023, [PDF](https://ianproberts.com/pdf/pub/simodel.pdf), [GitHub](https://ianproberts.com/simodel).

[4] I. P. Roberts, A. Chopra, T. Novlan, S. Vishwanath, and J. G. Andrews, "STEER: Beam Selection for Full-Duplex Millimeter Wave Communication Systems," _IEEE Trans. Commun._, Oct. 2022, [PDF](https://ianproberts.com/pdf/pub/steer.pdf), [GitHub](https://ianproberts.com/steer).

These measurements of self-interference were taken at 28 GHz in an anechoic chamber using two colocated 256-element phased arrays mounted on separate sides of an equilateral triangular platform. Please see [2] and [3] for details for a summary of these measurements.

This related work can be found at https://ianproberts.com.

# What is Self-Interference? 

When a transceiver (a wireless device) attempts to transmit and receive at the same time using the same frequency spectrum, some of its transmitted signal will leak into its receiver, corrupting reception of a desired signal.
This undesired leakage (or coupling) is called "self-interference". 
Transmit and receiving at the same time using the same frequency spectrum is called "full-duplex" operation.

This work is particularly focused on self-interference in full-duplex systems operating at mmWave frequencies (roughly 30 GHz to 100 GHz).
In mmWave systems, dense antenna arrays containing dozens or even hundreds of individual antennas are used to overcome high path loss at these high carrier frequencies. 
They do this by forming very directional beams, focusing energy in a particular direction to increase received signal power. 

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/222208767-c359fd9d-0fe0-4814-a56b-46d1f3fd306d.svg"/>
</p>

This work is interested in how much self-interference is coupled in full-duplex mmWave systems when using particular transmit and receive beams. Some transmit and receive beams will couple higher self-interference than others; this depends on the steering direction of the beams and on the (unknown) underlying self-interference channel. 

# What is STEER?

Modern mmWave communication systems rely on beam alignment to deliver sufficient beamforming gain to close the link between devices. 
In [1], we present the first beam selection methodology for full-duplex mmWave systems, called STEER, that delivers high beamforming gain while significantly reducing the self-interference coupled between the transmit and receive beams. 

STEER does not necessitate changes to conventional beam alignment methodologies nor additional over-the-air feedback, making it compatible with existing cellular standards. 
Instead, STEER uses conventional beam alignment to identify the general directions beams should be steered, and then it makes use of a minimal number of self-interference measurements to jointly select transmit and receive beams that deliver high gain in these directions while coupling low self-interference. 

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/222922193-67646503-455b-45bf-a573-348fc06e3703.svg"/>
</p>

In [1], STEER was implemented on an actual industry-grade 28 GHz phased array platform, which showed that full-duplex operation with beams selected by STEER can notably outperform both half-duplex and full-duplex operation with beams chosen via conventional beam selection. 
For instance, STEER can reliably reduce self-interference by more than 20 dB and improve SINR by more than 10 dB, compared to conventional beam selection. 
Our experimental results highlight that beam alignment can be used not only to deliver high beamforming gain in full-duplex mmWave systems but also to mitigate self-interference to levels near or below the noise floor, rendering additional self-interference cancellation unnecessary with STEER. 

# How STEER Works

STEER leverages our observations from a recent measurement campaign of mmWave self-interference [2], [3], which showed that small shifts in
the steering directions of the transmit and receive beams (on the order of one degree) can lead to noteworthy reductions in self-interference. 
STEER makes use of self-interference measurements across small spatial neighborhoods to jointly select transmit and receive beams at the full-duplex device that offer reduced self-interference while delivering high beamforming gain on the uplink and downlink.

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/222922156-2ad6a1c4-55b7-4ee6-9b21-78f32d1b044a.svg"/>
</p>

Put simply, STEER slightly shifts the transmit and receive beams at a full-duplex mmWave transceiver to reduce self-interference---with the goal being that these slight shifts do not prohibitively degrade downlink or uplink quality (i.e., SNR).
By maintaining high SNR and reducing self-interference, high SINRs can be achieved with STEER.

# Contents

This repo contains the following MATLAB code:
 - a main script `main.m` illustrating example usage
 - a function `construct_neighborhood.m` that can be used to construct a spatial neighborhood

# Example Usage

We will now walk through `main.m` to summarize its usage.

Suppose a full-duplex mmWave base station employs a codebook of transmit beams that it uses to conduct downlink beam alignment and a codebook of receive beams that it uses to conduct uplink beam alignment.
Whichever beams get selected from these codebooks via conventional beam alignment will be used to initialize STEER.

For each transmit-receive beam pair, we can run STEER to jointly select slightly shifted versions of these beams.
When running STEER, there are five design parameters:
- size of the transmit neighborhood
- resolution of the transmit neighborhood
- size of the receive neighborhood
- resolution of the receive neighborhood
- a self-interference target

Note that in [1], we assumed the size and resolution of the transmit and receive neighborhoods to be equal but this can be generalized straightforwardly.

### Define Transmit and Receive Codebooks

Before running STEER, the transmit and receive codebooks used for conventional beam alignment at the full-duplex base station must be defined.
This can be done by defining the steering directions of the codebooks' beams. 
The steering direction of each beam contains a component in azimuth and elevation, as illustrated below.

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/222922576-377d2020-43aa-4f13-9db8-b9ada9246846.svg"/>
</p>

In the example below, the transmit and receive codebooks are identical. Each codebook has beams that span in azimuth from -56 deg. to 56 deg. in 8 deg. steps and in elevation from -8 deg. to 8 deg. in 8 deg. steps. This amounts to a total of 45 beams in each codebook, meaning there are 2025 transmit-receive beam pairs.

```
% codebooks span azimuth and elevation
az_deg = [-56:8:56].';
el_deg = [-8:8:8].';

% assume TX and RX codebooks to be the same
txcb_az = az_deg;
txcb_el = el_deg;
rxcb_az = az_deg;
rxcb_el = el_deg;
```

The azimuth-elevation of each transmit and receive steering direction can then be populated as follows.

```
% full TX codebook in az-el
num_tx_az = length(txcb_az);
num_tx_el = length(txcb_el);
txcb_azel = [repmat(txcb_az,num_tx_el,1) repelem(txcb_el,num_tx_az,1)];

% full RX codebook in az-el
num_rx_az = length(rxcb_az);
num_rx_el = length(rxcb_el);
rxcb_azel = [repmat(rxcb_az,num_rx_el,1) repelem(rxcb_el,num_rx_az,1)];
```

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

# Acknowledgments

This work has been supported by the National Science Foundation Graduate Research Fellowship Program (Grant No. DGE-1610403). Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
