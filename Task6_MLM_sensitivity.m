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

% --- Define ranges
r_range = logspace(log10(0.9), log10(1.1), 30) * 2010;            % Shear modulus range
eta_range = logspace(-5, 5, 30) * 1e18;            % Viscosity range

% --- Initialize arrays for Love numbers (real and imaginary)
k2_real = zeros(length(r_range), length(eta_range));
k2_imag = zeros(length(r_range), length(eta_range));
h2_real = zeros(length(r_range), length(eta_range));
h2_imag = zeros(length(r_range), length(eta_range));


% Orbital parameters
T = 87.969216879 * 24 * 3600; % Period in seconds
omega0 = 2 * pi / T;

% Loop over mu and eta
for i = 1:length(r_range)
    for j = 1:length(eta_range)
        % Inner core (solid)
        Interior_Model(1).R0 = 1e-3 * 1180; %dummy layer
        Interior_Model(1).rho0 = 7380;
        
        Interior_Model(2).R0    = 1180;         % meters
        Interior_Model(2).rho0  = 7380;           % kg/m³
        Interior_Model(2).mu0   = 10.0e10;         % Pa
        Interior_Model(2).eta0  = 1e20;           % Pa·s
        
        % Outer core (liquid)
        Interior_Model(3).R0    = r_range(i);         % meters
        Interior_Model(3).rho0  = 6470;           % kg/m³
        Interior_Model(3).ocean  = 1;           % very low (fluid)
        
        % Mantle (viscoelastic layer)
        Interior_Model(4).R0    = 2410;         % meters (planetary radius)
        Interior_Model(4).rho0  = 3100;           % kg/m³ (can vary)
        Interior_Model(4).mu0   = 6e10;         % Pa shear modulus
        Interior_Model(4).eta0  = eta_range(j);           % Pa·s
        
        % Crust (elastic lid)
        Interior_Model(5).R0    = 2439.36;     % meters (surface)
        Interior_Model(5).rho0  = 2850;          % kg/m³
        Interior_Model(5).mu0   = 5.5e10;        % Pa
        Interior_Model(5).eta0  = 1e23;  

        % Forcing
        Forcing(1).Td = T;
        Forcing(1).n = 2;
        Forcing(1).m = 0;
        Forcing(1).F = 1;

        % Numerics setup
        Numerics.Nlayers = length(Interior_Model);
        Numerics.method = 'variable';
        Numerics.Nrbase = 200;
        Numerics.parallel_sol = 0;
        Numerics.parallel_gen = 0;
        Numerics.perturbation_order = 2;
        Numerics.solution_cutoff = 12;
        Numerics.load_couplings = 1;
        Numerics.Nenergy = 12;
        Numerics.rheology_cutoff = 2;

        % Setup boundary and rheology
        [Numerics, Interior_Model] = set_boundary_indices(Numerics, Interior_Model, 'verbose');
        Interior_Model_aux = get_rheology(Interior_Model, Numerics, Forcing);

        % Compute tidal response
        [Love_Spectra,~] = get_Love(Interior_Model_aux, Forcing, Numerics, 'verbose');

        % Store results
        k2_real(i,j) = real(Love_Spectra.k);
        k2_imag(i,j) = imag(Love_Spectra.k);
        h2_real(i,j) = real(Love_Spectra.h);
        h2_imag(i,j) = imag(Love_Spectra.h);

    end
end
%%
% --- Plotting ---
figure;

% --- Real part of k2 ---
subplot(2,2,1)
levels = linspace(0, 1.5, 10);
[~, h] = contourf(eta_range / 1e18, r_range / 2010, k2_real, levels);
shading interp
hold on
[C, h_lines] = contour(eta_range / 1e18, r_range / 2010, k2_real, levels, 'LineColor', 'k');
clabel(C, h_lines, 'FontSize', 10, 'Color', 'k', 'LabelSpacing', 400);
set(gca, 'XScale', 'log', 'FontSize', 12)
xlabel('Viscosity $\eta / \eta_0$', 'Interpreter', 'latex')
xlim([10^-5, 10^5])
ylabel('Core Radius $R / R_0$', 'Interpreter', 'latex')
ylim([0.9, 1.1])
title('Real part of $k_2$', 'Interpreter', 'latex')
colorbar

% --- Imaginary part of k2 ---
subplot(2,2,2)
levels = linspace(min(k2_imag(:)), 0, 5);
[~, h] = contourf(eta_range / 1e18, r_range / 2010, k2_imag, levels);
shading interp
hold on
[C, h_lines] = contour(eta_range / 1e18, r_range / 2010, k2_imag, levels, 'LineColor', 'k');
clabel(C, h_lines, 'FontSize', 10, 'Color', 'k', 'LabelSpacing', 400);
clim([min(k2_imag(:)), 0]);
set(gca, 'XScale', 'log', 'FontSize', 12)
xlabel('Viscosity $\eta / \eta_0$', 'Interpreter', 'latex')
xlim([10^-5, 10^5])
ylim([0.9, 1.1])
ylabel('Core Radius $R / R_0$', 'Interpreter', 'latex')
title('Imaginary part of $k_2$', 'Interpreter', 'latex')
colorbar

% --- Real part of h2 ---
subplot(2,2,3)
levels = linspace(0, 2.5, 10);
[~, h] = contourf(eta_range / 1e18, r_range / 2010, h2_real, levels);
shading interp
hold on
[C, h_lines] = contour(eta_range / 1e18, r_range / 2010, h2_real, levels, 'LineColor', 'k');
clabel(C, h_lines, 'FontSize', 10, 'Color', 'k', 'LabelSpacing', 400);
set(gca, 'XScale', 'log', 'FontSize', 12)
xlabel('Viscosity $\eta / \eta_0$', 'Interpreter', 'latex')
ylabel('Core Radius $R / R_0$', 'Interpreter', 'latex')
xlim([10^-5, 10^5])
%ylim([0.9, 1.1])
title('Real part of $h_2$', 'Interpreter', 'latex')
colorbar

% --- Imaginary part of h2 ---
subplot(2,2,4)
levels = linspace(min(h2_imag(:)), 0, 5);
[~, h] = contourf(eta_range / 1e18, r_range / 2010, h2_imag, levels);
shading interp
hold on
[C, h_lines] = contour(eta_range / 1e18, r_range / 2010, h2_imag, levels, 'LineColor', 'k');
clabel(C, h_lines, 'FontSize', 10, 'Color', 'k', 'LabelSpacing', 400);
clim([min(h2_imag(:)), 0]);
set(gca, 'XScale', 'log', 'FontSize', 12)
xlabel('Viscosity $\eta / \eta_0$', 'Interpreter', 'latex')
ylabel('Core Radius $R / R_0$', 'Interpreter', 'latex')
xlim([10^-5, 10^5])
%ylim([0.9, 1.1])
title('Imaginary part of $h_2$', 'Interpreter', 'latex')
colorbar