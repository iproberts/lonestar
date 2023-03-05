classdef channel < matlab.mixin.Copyable
    % CHANNEL A generic MIMO channel.
    properties
        name; % human-readable identifier
        channel_matrix; % channel matrix
        num_antennas_transmit; % number of transmit antennas
        num_antennas_receive; % number of receive antennas
        array_transmit; % transmit array at the channel input
        array_receive; % receive array at the channel output
        carrier_frequency; % carrier frequency (Hz)
        carrier_wavelength; % carrier wavelength (m)
        propagation_velocity; % propagation velocity (m/s)
        normalized_channel_energy; % total energy the channel matrix has after normalization
        force_channel_energy_normalization; % a boolean indicating if the channel matrix should be normalized
    end
    methods(Static)
        function c = create(type)
            % CREATE Creates a channel object of a specific type.
            %
            % Usage:
            %  c = channel.create()
            %  c = channel.create(type)
            %
            % Args:
            %  type: (optional) a string specifying what type of channel to
            %  create
            %
            % Returns:
            %  c: a channel object of the type specified
            if nargin < 1 || isempty(type)
                type = 'default';
            end
            if strcmpi(type,'default') || strcmpi(type,'')
                c = channel();
            elseif strcmpi(type,'Rayleigh')
                c = channel_rayleigh();
            elseif strcmpi(type,'LOS') || strcmpi(type,'line-of-sight')
                c = channel_line_of_sight();
            elseif strcmpi(type,'Rician') || strcmpi(type,'Ricean')
                c = channel_rician();
            elseif strcmpi(type,'spherical-wave') || strcmpi(type,'spherical_wave') || strcmpi(type,'spherical')
                c = channel_spherical_wave();
            elseif strcmpi(type,'ray-cluster') || strcmpi(type,'ray_cluster') || strcmpi(type,'ray/cluster')
                c = channel_ray_cluster();
            else
                error('Invalid channel type specifier.');
            end
        end
    end
    methods
        function obj = channel(name)
            % CHANNEL Creates a MIMO channel object.
            % 
            % Usage:
            %  obj = CHANNEL()
            %  obj = CHANNEL(name)
            % 
            % Args:
            %  name: an optional name for the object
            % 
            % Returns:
            %  obj: an object representing a MIMO channel
            if nargin < 1
                name = [];
            end
            obj.set_name(name);
            obj.initialize();
        end
        
        function initialize(obj)
            % INITIALIZE Initializes a channel.
            %
            % Usage:
            %  INITIALIZE()
            obj.set_arrays(array(1),array(1));
            obj.set_carrier_frequency(1);
            obj.set_propagation_velocity(1);
            obj.set_force_channel_energy_normalization(false);
        end
        
        % -----------------------------------------------------------------
        % Set functions.
        % -----------------------------------------------------------------
        
        function set_name(obj,name)
            % SET_NAME Sets the name of the channel.
            %
            % Usage:
            %  SET_NAME()
            %  SET_NAME(name)
            %
            % Args:
            %  name: (optional) a string; if not passed, 'channel' is
            %  the default name used
            if nargin < 2 || isempty(name)
                name = 'channel';
            end
            obj.name = name;
        end
        
        function set_propagation_velocity(obj,val)
            % SET_PROPAGATION_VELOCITY Sets the propagation velocity of the
            % channel.
            %
            % Usage:
            %  SET_PROPAGATION_VELOCITY(val)
            %
            % Args:
            %  val: propagation velocity (meters/sec)
            obj.propagation_velocity = val;
            obj.carrier_wavelength = val / obj.carrier_frequency;
        end
        
        function set_carrier_frequency(obj,fc)
            % SET_CARRIER_FREQUENCY Sets the carrier frequency of the
            % channel.
            %
            % Usage:
            %  SET_CARRIER_FREQUENCY(fc)
            %
            % Args:
            %  fc: carrier frequency (Hz)
            %
            % Notes:
            %  : Also updates carrier wavelength.
            obj.carrier_frequency = fc;
            obj.carrier_wavelength = obj.propagation_velocity / fc;
        end
        
        function set_arrays(obj,array_transmit,array_receive)
            % SET_ARRAYS Sets the transmit and receive arrays at the input
            % and output of the channel.
            %
            % Usage:
            %  SET_ARRAYS(array_transmit,array_receive)
            %
            % Args:
            %  array_transmit: an array object at the channel input
            %  array_receive: an array object at the channel output
            obj.set_array_transmit(array_transmit);
            obj.set_array_receive(array_receive);
        end
        
        function set_array_transmit(obj,array)
            % SET_ARRAY_TRANSMIT Sets the transmit array object. Also sets
            % the number of transmit antennas accordingly.
            % 
            % Usage:
            %  SET_ARRAY_TRANSMIT(array)
            % 
            % Args:
            %  array: an array object
            obj.array_transmit = array;
            obj.num_antennas_transmit = array.num_antennas;
            obj.set_normalized_channel_energy();
        end
        
        function set_array_receive(obj,array)
            % SET_ARRAY_RECEIVE Sets the receive array object. Also sets
            % the number of receive antennas accordingly.
            % 
            % Usage:
            %  SET_ARRAY_RECEIVE(array)
            % 
            % Args:
            %  array: an array object
            obj.array_receive = array;
            obj.num_antennas_receive = array.num_antennas;
            obj.set_normalized_channel_energy();
        end
                
        function set_channel_matrix(obj,H)
            % SET_CHANNEL_MATRIX Sets the channel matrix.
            %
            % Usage:
            %  SET_CHANNEL_MATRIX(H)
            % 
            % Args:
            %  H: channel matrix
            if obj.force_channel_energy_normalization
                H = obj.enforce_channel_energy_normalization(H);
            end
            obj.channel_matrix = H;
        end
        
        function set_force_channel_energy_normalization(obj,force)
            % SET_FORCE_CHANNEL_ENERGY_NORMALIZATION Sets the enforcement
            % of channel energy normalization. If true, the channel matrix
            % will always be normalized such that its energy is of the 
            % desired value.
            %
            % Usage:
            %  SET_FORCE_CHANNEL_ENERGY_NORMALIZATION(force)
            %
            % Args:
            %  force: a boolean indicating if the channel matrix should be
            %  normalized or not
            obj.force_channel_energy_normalization = logical(force);
        end
                
        function set_normalized_channel_energy(obj,E)
            % SET_NORMALIZED_CHANNEL_ENERGY Sets the normalized energy of
            % the channel.
            %
            % Usage:
            %  SET_NORMALIZED_CHANNEL_ENERGY()
            %  SET_NORMALIZED_CHANNEL_ENERGY(E)
            %
            % Args:
            %  E: (optional) the desired normalized channel energy; if not
            %  passed, the product of the number of transmit antennas and
            %  the number of receive antennas will be used
            if nargin < 2 || isempty(E)
                Nt = obj.num_antennas_transmit;
                Nr = obj.num_antennas_receive;
                E = Nt * Nr;
            end
            obj.normalized_channel_energy = E;
        end
        
        % -----------------------------------------------------------------
        % Get functions.
        % -----------------------------------------------------------------
        
        function H = get_channel_matrix(obj)
            % GET_CHANNEL_MATRIX Returns the channel matrix.
            %
            % Usage:
            %  H = GET_CHANNEL_MATRIX()
            % 
            % Returns:
            %  H: channel matrix
            H = obj.channel_matrix;
        end
        
        % -----------------------------------------------------------------
        % Enforce functions.
        % -----------------------------------------------------------------
        
        function G = enforce_channel_energy_normalization(obj,H)
            % ENFORCE_CHANNEL_ENERGY_NORMALIZATION Normalizes the channel
            % matrix so that its total energy (squared Frobenius norm) is
            % equal the current normalized channel energy property. The
            % default normalized channel energy is the product of the 
            % number of transmit antenans and the number of receive
            % antennas.
            %
            % Usage:
            %  G = ENFORCE_CHANNEL_ENERGY_NORMALIZATION()
            %  G = ENFORCE_CHANNEL_ENERGY_NORMALIZATION(H)
            %
            % Args:
            %  H: (optional) a channel matrix; if not passed, the current
            %  channel matrix will be used and overwritten with the 
            %  normalized version; if passed, the channel matrix property 
            %  will not be set
            %
            % Returns:
            %  G: the normalized channel matrix
            if nargin < 2 || isempty(H)
                H = obj.get_channel_matrix();
            end
            E = obj.normalized_channel_energy;
            val = norm(H,'fro');
            G = H ./ val * sqrt(E);
            if nargin < 2 || isempty(H)
                obj.set_channel_matrix(G); % only set if H was not passed
            end
        end
        
        % -----------------------------------------------------------------
        % Shorthand functions.
        % -----------------------------------------------------------------
        
        function val = Nt(obj)
            % Nt Returns the number of transmit antennas into the channel.
            %
            % Usage:
            %  val = Nt()
            %
            % Returns:
            %  val: the number of transmit antennas into the channel
            val = obj.num_antennas_transmit;
        end
        
        function val = Nr(obj)
            % Nr Returns the number of receive antennas out of the channel.
            %
            % Usage:
            %  val = Nr()
            %
            % Returns:
            %  val: the number of receive antennas out of the channel
            val = obj.num_antennas_receive;
        end
        
        function val = H(obj)
            % H Returns the channel matrix.
            %
            % Usage:
            %  val = H()
            %
            % Returns:
            %  val: the channel matrix
            val = obj.get_channel_matrix();
        end
        
        % -----------------------------------------------------------------
        % Legacy functions.
        % -----------------------------------------------------------------
        
        function set_transmit_array(obj,array)
            % SET_TRANSMIT_ARRAY Sets the transmit array object (LEGACY).
            % 
            % Usage:
            %  SET_TRANSMIT_ARRAY(array)
            % 
            % Args:
            %  array: an array object
            obj.set_array_transmit(array);
        end
        
        function set_receive_array(obj,array)
            % SET_RECEIVE_ARRAY Sets the receive array object (LEGACY).
            % 
            % Usage:
            %  SET_RECEIVE_ARRAY(array)
            % 
            % Args:
            %  array: an array object
            obj.set_array_receive(array);
        end
    end
end