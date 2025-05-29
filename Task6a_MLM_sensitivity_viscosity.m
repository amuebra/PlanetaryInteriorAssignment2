clearvars
clc
set(0,'defaulttextInterpreter','latex') 

% Define range
N = 30;
% Define constant viscosity
mantle_const = 1e18;
crust_const = 1e20;
inner_const = 1e23;

inner_range = logspace(-5, 5, N) * inner_const;
mantle_range = logspace(-5, 5, N) * mantle_const;
crust_range = logspace(-5, 5, N) * crust_const;

% Initialize result arrays
k2_inner = zeros(N,1);
k2_crust = zeros(N,1);
k2_mantle = zeros(N,1);

%% LOOP 1: Vary inner core viscosity
for i = 1:N
    Interior_Model = setup_model(inner_range(i), mantle_const, crust_const);
    [k2_inner(i), h2_inner(i)] = compute_k2(Interior_Model);
end

%% LOOP 2: Vary outer core viscosity
for i = 1:N
    Interior_Model = setup_model(inner_const, mantle_range(i), crust_const);
    [k2_mantle(i), h2_mantle(i)] = compute_k2(Interior_Model);
end

%% LOOP 3: Vary mantle viscosity
for i = 1:N
    Interior_Model = setup_model(inner_const, mantle_const, crust_range(i));
    [k2_crust(i), h2_crust(i)] = compute_k2(Interior_Model);
end

%% Plotting k_2
close all;
figure('units', 'points', 'Position', [0, 0, 455.2441, 0.8*455.2441]);
hold on
% Normalized x-axis: divide each by its central value
plot(inner_range / inner_const, real(k2_inner), 'r-', 'LineWidth', 2)
plot(mantle_range / mantle_const, real(k2_mantle), 'b--', 'LineWidth', 2)
plot(crust_range / crust_const, real(k2_crust), 'g-.', 'LineWidth', 2)
xlabel('Normalized Viscosity ($\eta / \eta_0$)', 'Interpreter', 'latex', 'FontSize', 14)
ylim([0,1]);
set(gca, 'XScale', 'log')
ylabel('Real part of $k_2$', 'Interpreter', 'latex', 'FontSize', 14)
legend({'Inner Core', 'Mantle', 'Crust'}, 'Location', 'best')
% title('$k_2$ sensitivity to layer viscosity', 'Interpreter', 'latex')
grid on
savefig('Assignment2/2.1_viscosity_norm.fig')
saveas(gcf, 'Assignment2/2.1_viscosity_norm.svg')

figure('units', 'points', 'Position', [0, 0, 455.2441, 0.8*455.2441]);
hold on
% Normalized x-axis: divide each by its central value
plot(inner_range / inner_const, imag(k2_inner), 'r-', 'LineWidth', 2)
plot(mantle_range / mantle_const, imag(k2_mantle), 'b--', 'LineWidth', 2)
plot(crust_range / crust_const, imag(k2_crust), 'g-.', 'LineWidth', 2)
xlabel('Normalized Viscosity ($\eta / \eta_0$)', 'Interpreter', 'latex', 'FontSize', 14)
% ylim([0,1]);
set(gca, 'XScale', 'log')
ylabel('Imaginary part of $k_2$', 'Interpreter', 'latex', 'FontSize', 14)
legend({'Inner Core', 'Mantle', 'Crust'}, 'Location', 'best')
% title('$k_2$ sensitivity to layer viscosity', 'Interpreter', 'latex')
grid on
savefig('Assignment2/2.1_viscosity_norm_imag.fig')
saveas(gcf, 'Assignment2/2.1_viscosity_norm_imag.svg')

%% Plotting h_2
figure('units', 'points', 'Position', [0, 0, 455.2441, 0.8*455.2441]);
hold on
% Normalized x-axis: divide each by its central value
plot(inner_range / inner_const, real(h2_inner), 'r-', 'LineWidth', 2)
plot(mantle_range / mantle_const, real(h2_mantle), 'b--', 'LineWidth', 2)
plot(crust_range / crust_const, real(h2_crust), 'g-.', 'LineWidth', 2)
xlabel('Normalized Viscosity ($\eta / \eta_0$)', 'Interpreter', 'latex', 'FontSize', 14)
% ylim([0,1]);
set(gca, 'XScale', 'log')
ylabel('Real part of $h_2$', 'Interpreter', 'latex', 'FontSize', 14)
legend({'Inner Core', 'Mantle', 'Crust'}, 'Location', 'best')
% title('$k_2$ sensitivity to layer viscosity', 'Interpreter', 'latex')
grid on
savefig('Assignment2/2.1_viscosity_norm_h.fig')
saveas(gcf, 'Assignment2/2.1_viscosity_norm_h.svg')

figure('units', 'points', 'Position', [0, 0, 455.2441, 0.8*455.2441]);
hold on
% Normalized x-axis: divide each by its central value
plot(inner_range / inner_const, imag(h2_inner), 'r-', 'LineWidth', 2)
plot(mantle_range / mantle_const, imag(h2_mantle), 'b--', 'LineWidth', 2)
plot(crust_range / crust_const, imag(h2_crust), 'g-.', 'LineWidth', 2)
xlabel('Normalized Viscosity ($\eta / \eta_0$)', 'Interpreter', 'latex', 'FontSize', 14)
% ylim([0,1]);
set(gca, 'XScale', 'log')
ylabel('Imaginary part of $h_2$', 'Interpreter', 'latex', 'FontSize', 14)
legend({'Inner Core', 'Mantle', 'Crust'}, 'Location', 'best')
% title('$k_2$ sensitivity to layer viscosity', 'Interpreter', 'latex')
grid on
savefig('Assignment2/2.1_viscosity_norm_h_imag.fig')
saveas(gcf, 'Assignment2/2.1_viscosity_norm_h_imag.svg')

%% Functions
function Interior_Model = setup_model(eta_inner, eta_mantle, eta_crust)
    % Inner core (solid)
    Interior_Model(1).R0 = 1e-3 * 1180; % dummy layer
    Interior_Model(1).rho0 = 7380;

    Interior_Model(2).R0 = 1180;
    Interior_Model(2).rho0 = 7380;
    Interior_Model(2).mu0 = 10e10;
    Interior_Model(2).eta0 = eta_inner;

    % Outer core
    Interior_Model(3).R0 = 2010;
    Interior_Model(3).rho0 = 6470;
    Interior_Model(3).ocean = 1;

    % Mantle
    Interior_Model(4).R0 = 2410;
    Interior_Model(4).rho0 = 3100;
    Interior_Model(4).mu0 = 6e10;
    Interior_Model(4).eta0 = eta_mantle;

    % Crust
    Interior_Model(5).R0 = 2439.36;
    Interior_Model(5).rho0 = 2850;
    Interior_Model(5).mu0 = 5.5e10;
    Interior_Model(5).eta0 = eta_crust;
end

function [k2_value, h2_value] = compute_k2(Interior_Model)
    T = 87.969216879*24*3600;
    Forcing.Td = T;
    Forcing.n = 2; 
    Forcing.m = 0; 
    Forcing.F = 1;

    Numerics.Nlayers = length(Interior_Model);
    Numerics.method = 'variable';
    Numerics.Nrbase = 200;
    Numerics.parallel_sol = 0;
    Numerics.parallel_gen = 0;
    Numerics.perturbation_order = 2;

    [Numerics, Interior_Model] = set_boundary_indices(Numerics, Interior_Model, 'verbose');
    Interior_Model_aux = get_rheology(Interior_Model, Numerics, Forcing);
    [Love_Spectra, ~] = get_Love(Interior_Model_aux, Forcing, Numerics, 'verbose');
    k2_value = Love_Spectra.k;
    h2_value = Love_Spectra.h;
end
