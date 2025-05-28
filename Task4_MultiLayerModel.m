clearvars
clc
set(0,'defaulttextInterpreter','latex') 
mfile_name          = mfilename('fullpath');
if contains(mfile_name,'LiveEditorEvaluationHelper')
    mfile_name=matlab.desktop.editor.getActiveFilename;
end
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);
cd('..')
addpath(genpath(pwd))

%% Mercury's internal parameter

% Inner core (solid)
Interior_Model(1).R0 = 1e-3 * 1180; %dummy layer
Interior_Model(1).rho0 = 7380;

Interior_Model(2).R0    = 1180;         % meters
Interior_Model(2).rho0  = 7380;           % kg/m³
%Interior_Model(2).Ks0   = 1.91e11;        % Pa
Interior_Model(2).mu0   = 10.0e10;         % Pa
Interior_Model(2).eta0  = 1e20;           % Pa·s

% Outer core (liquid)
Interior_Model(3).R0    = 2010;         % meters
Interior_Model(3).rho0  = 6470;         % kg/m³
%Interior_Model(3).Ks0   = 1.15e11;     % Pa
%Interior_Model(3).mu0   = 3.3e-1;      % Pa
Interior_Model(3).ocean  = 1;           % very low (fluid)

% Mantle (viscoelastic layer)
Interior_Model(4).R0    = 2410;         % meters (planetary radius)
Interior_Model(4).rho0  = 3100;         % kg/m³ (can vary)
%Interior_Model(4).Ks0   = 1.3e11;      % Pa
Interior_Model(4).mu0   = 6e10;         % Pa shear modulus
Interior_Model(4).eta0  = 1e18;         % Pa·s

% Crust (elastic lid)
 Interior_Model(5).R0    = 2439.36;     % meters (surface)
 Interior_Model(5).rho0  = 2850;        % kg/m³
% Interior_Model(5).Ks0   = 6.0e10;     % Pa
 Interior_Model(5).mu0   = 5.5e10;      % Pa
 Interior_Model(5).eta0  = 1e23;        % Pa·s (effectively elastic)

%% Define Forcing
% Mercurys Frequency
T = 87.969216879*24*3600; %Mercury's orbital period
omega0 = 2*pi / T; %Mercury's orbital frequency

Forcing(1).Td=T;
Forcing(1).n=2; 
Forcing(1).m=0; 
Forcing(1).F=1;

%% Define Numerics
%radial discretization
Numerics.Nlayers = length(Interior_Model); % number of concentric layers. Including the core!
Numerics.method = 'variable'; % method of setting the radial points per layer, here fixed number of layers
Numerics.Nrbase = 200; % depending on the method this will determine the number of points per layer
%code parallelization
Numerics.parallel_sol = 0; % Use a parfor-loop to call get_Love, either 0 or 1
Numerics.parallel_gen = 0; % Calculate potential coupling files and the propagation inside get_solution using parfor-loops, either 0 or 1
% lateral variations
Numerics.perturbation_order = 2; %maximum order to which couplings are considered
Numerics.solution_cutoff = 12; % maximum degree of solution, not used if perturbation order is specified
Numerics.load_couplings = 1; % 0=no loading, 1=loading of specific file, 2=searches for big enough file
Numerics.Nenergy = 12; % maximum degree to which energy dissipation is expanded 
Numerics.rheology_cutoff = 2; % maximum order of difference (so in log) up to which rheology is still used 
[Numerics, Interior_Model] = set_boundary_indices(Numerics, Interior_Model,'verbose');


%% Compute the tidal response
% write interior model in the right format for the code
Interior_Model= get_rheology(Interior_Model,Numerics,Forcing);
[Love_Spectra,y]=get_Love(Interior_Model,Forcing,Numerics,'verbose');
disp(['k_2 Model: ' num2str(Love_Spectra.k)])
disp(['h Love number LOV3D: ' num2str(Love_Spectra.h)])

