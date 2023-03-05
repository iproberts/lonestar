classdef channel_spherical_wave < channel
    properties
        large_scale_gain;
    end
    methods
        function obj = channel_spherical_wave(name)
            % CHANNEL_SPHERICAL_WAVE Creates a spherical-wave MIMO 
            % channel object.
            % 
            % Usage:
            %  obj = CHANNEL_SPHERICAL_WAVE()
            %  obj = CHANNEL_SPHERICAL_WAVE(name)
            % 
            % Args:
            %  name: an optional name for the object
            % 
            % Returns:
            %  obj: an object representing a spherical-wave MIMO channel
            if nargin < 1 || isempty(name)
                name = 'channel-spherical-wave-mimo';
            end
            obj.name = name;
        end
        
        function r = get_range_between_elements(obj,idx_tx,idx_rx)
            % GET_RANGE_BETWEEN_ELEMENTS Returns the range between an
            % element in the transmit array and an element in the receive
            % array.
            %
            % Usage:
            %  r = GET_RANGE_BETWEEN_ELEMENTS(idx_tx,idx_rx)
            %
            % Args:
            %  idx_tx: transmit array element index
            %  idx_rx: receive array element index
            %
            % Returns:
            %  r: the distance between the two elements (meters)
            dx = obj.array_transmit.x(idx_tx) - obj.array_receive.x(idx_rx);
            dy = obj.array_transmit.y(idx_tx) - obj.array_receive.y(idx_rx);
            dz = obj.array_transmit.z(idx_tx) - obj.array_receive.z(idx_rx);
            r = sqrt(dx.^2 + dy.^2 + dz.^2) * obj.carrier_wavelength;
        end
        
        function H = realization(obj)
            H = obj.channel_realization();
            obj.set_channel_matrix(H);
            H = obj.get_channel_matrix();
        end
        
        function H = channel_realization(obj)
            % CHANNEL_REALIZATION Realizes a random instance of the channel
            % matrix based on the channel object's parameters.
            %
            % Usage:
            %  H = CHANNEL_REALIZATION()
            %
            % Returns:
            %  H: the resulting channel matrix
            Nr = obj.Nr;
            Nt = obj.Nt;
            lam = obj.carrier_wavelength;
            H = zeros(Nr,Nt);
            for i = 1:Nr
                for j = 1:Nt
                    r = obj.get_range_between_elements(j,i);
                    val = 1/r * exp(-1j*2*pi/lam*r);
                    H(i,j) = val;
                end
            end
            n = norm(H,'fro').^2;
            if isinf(n)
                warning('Frobenius norm squared of the channel is infinite. Normalization is going to return a zero matrix.');
            end
            H = H ./ sqrt(n) * sqrt(Nt * Nr); % normalize the channel's Frob norm
            % n = norm(H,'fro').^2;
            obj.set_channel_matrix(H);
            obj.large_scale_gain = sqrt(n ./ (Nt * Nr));
        end
    end
end