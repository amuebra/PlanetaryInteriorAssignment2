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

%--- Work in dimensional units
R0 = 2439.36;% Radius in km
rho0 = 5429.3; % density in kg/m^3
mu0 = 55e9; % shear modulus in Pa
eta0 = 3.16e18; % viscosity in Pa.s
G0=6.67430E-11;

% Define interior model: one homogeneous viscoelastic layer
Interior_Model(1).R0 = 1e-3 * R0; %dummy layer
Interior_Model(1).rho0 = rho0;
Interior_Model(2).R0 = R0;      % surface radius in km
Interior_Model(2).rho0 = rho0;  % density in kg/m^3
Interior_Model(2).mu0 = mu0;    % shear modulus in Pa
Interior_Model(2).eta0 = eta0;  % viscosity in Pa.s


T = 87.969216879*24*3600; %Mercury's orbital period
omega0 = 2*pi / T; %Mercury's orbital frequency

% For fluid-like behaviour use these values
% omega0=1e-12;
% T = 2*pi/omega0;

% Spherically-symmetric model assumed
% compute effective shear modulus
mu_eff=mu0/(rho0^2*(R0*1e3)^2*4/3*pi*G0);

%--- Forcing potential is defined, degree of forcing can be changed
Forcing(1).Td=T;
Forcing(1).n=2; 
Forcing(1).m=0; 
Forcing(1).F=1;

% Define numerics used to compute the tidal response
Numerics.Nlayers = length(Interior_Model); %for simple model
Numerics.method = 'variable'; % method of setting the radial points per layer
Numerics.Nrbase = 200; % depending on the method this will determine the number of points per layer
%code parallelization
Numerics.parallel_sol = 0; % Use a parfor-loop to call get_Love, either 0 or 1
Numerics.parallel_gen = 0; % Calculate potential coupling files and the propagation inside get_solution using parfor-loops, either 0 or 1
% lateral variations
Numerics.perturbation_order = 2; %set to 0, due to sphericallyharmonic body, maximum order to which couplings are considered
Numerics.solution_cutoff = 12; % maximum degree of solution, not used if perturbation order is specified
Numerics.load_couplings = 1; % 0=no loading, 1=loading of specific file, 2=searches for big enough file
Numerics.Nenergy = 12; % maximum degree to which energy dissipation is expanded 
Numerics.rheology_cutoff = 2; % maximum order of difference (so in log) up to which rheology is still used 
[Numerics, Interior_Model] = set_boundary_indices(Numerics, Interior_Model,'verbose');

%---Compute the tidal response
% write interior model in the right format for the code
Interior_Model = get_rheology(Interior_Model,Numerics,Forcing);
[Love_Spectra,y]=get_Love(Interior_Model,Forcing,Numerics,'verbose');

%---Compare against analytical expression
% obtain the Fourier-transformed effective shear modulus
mu_eff_hat=mu_eff*Interior_Model(2).muC;
% compute Love numbers using analytical expression
n=Forcing.n; 
mu_n=(2*n^2+4*n+3)/n*mu_eff_hat;
k2_analytic=1/(1+mu_n)*3/2/(n-1);
h2_analytic=1/(1+mu_n)*(2*n+1)/2/(n-1);
% display results
disp(['k Love number analytical expression: ' num2str(k2_analytic)])
disp(['k Love number LOV3D: ' num2str(Love_Spectra.k)])
disp(['Normalized difference: ' num2str((Love_Spectra.k-k2_analytic)/k2_analytic*100)  '%'])
disp(['h Love number analytical expression: ' num2str(h2_analytic)])
disp(['h Love number LOV3D: ' num2str(Love_Spectra.h)])
disp(['Normalized difference: ' num2str((Love_Spectra.h-h2_analytic)/h2_analytic*100)  '%'])
