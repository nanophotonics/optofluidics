clc
clear 
close all

%% Initial parameters
% *************************************************************************
eta = 377; % Ohm. Characteristic impedance of free space
eps0 = 8.85e-12; % F/m. Permittivity of free space
kb = 1.38e-23; % m^2*kg*s^-2*K^-1. Boltzmann constant

R = 9; % um. Radius of the HCPCF core
nmed = 1.33; % Refractive index of water
T = 273 + 70; % K. Ambient temperature

d_limits = [40, 20, 100]; % nm
lambda_limits = [500, 10, 1000]; % nm
power_limits = [100, 100, 1200]; % mW

limits = [d_limits; ...
          lambda_limits; ...
          power_limits;...
          ];

column_titles = {'Start:', 'Step:', 'End:'};
row_titles = {'Diameter (nm):',...
              'Wavelength (nm):',...
              'Power (mW):',...
              };
          
limits = dialog_table(row_titles,column_titles,limits);

d_limits = limits(1,:);
lambda_limits = limits(2,:);
power_limits = limits(3,:);

d = d_limits(1): d_limits(2) : d_limits(3); % nm. Diameter of the Au NP
lambda = lambda_limits(1) : lambda_limits(2) : lambda_limits(3); % nm. Wavelength of the laser in vacuum
power = power_limits(1) : power_limits(2) : power_limits(3); % mW of total laser power


%% Read Gold refractive index file and interpolate
% *************************************************************************
foldername = '';
filename = 'Gold Refractive-Rakic.txt';
n_data = dlmread([foldername filename], '\t', 1, 0);
% n_data(:,1) = wave (nm)
% n_data(:,2) = n
% n_data(:,3) = k
n_real = spline(n_data(:,1), n_data(:,2), lambda);
n_imaginary = spline(n_data(:,1), n_data(:,3), lambda);
n_particle = n_real + 1i*n_imaginary;

% figure
% plot(n_data(:,1), n_data(:,2), 'LineWidth', 2), hold all
% plot(n_data(:,1), n_data(:,3), 'LineWidth', 2), hold all
% xlabel('Wavelength (nm)')
% ylabel('Refractive index = n + i*k')
% legend('n','k', 'Location', 'NW')
% title(filename)
% xlim([200,1600])
% grid on

%% Calculate the polarisability
% *************************************************************************

% standard polarisability (SI units)
alpha_0 = zeros(size(d,2),size(lambda,2));
for i = 1:1:size(d,2)
    alpha_0(i,:) = 4*pi*(d(i)*1e-9/2)^3 .* ...
        (n_particle.^2 - nmed^2) ./ (n_particle.^2 + 2*nmed^2); 
end

% polarisability (SI units)
alpha = alpha_0;

%% Calculate the potential
% *************************************************************************

KE = kb*T;  % J, kinetic energy

% menu_profile = 2;
menu_profile = menu('Beam Profile?', 'Bessel', 'Gaussian');

if menu_profile == 1 % Bessel
    r = linspace(0, R, 100); % um. Radial coordinate inside fibre 
%     r = 0; % um. Radial coordinate inside fibre
    j12 = 0.269514; % Square of the Bessel function 1st kind 1st order at the first zero of the 0th order
    j01 = 2.40483; % First zero of the Bessel function 1st kind 0th order
    J02 = besselj(0,j01*r/R).^2; % Square of the Bessel function 1st kind 0th order
    z = 0;
elseif menu_profile == 2 % Gaussian
    r = linspace(0, 15, 100); % um. Radial coordinate
    z = 0:20:200; % um. Distance from the end of the fibre
    w0 = 5; % um. Beam diameter at the end of the fibre
    wz = zeros(size(lambda,2),size(z,2));
end

intensity = zeros(size(lambda,2),size(power,2),size(r,2),size(z,2),size(d,2));
potential = zeros(size(lambda,2),size(power,2),size(r,2),size(z,2),size(d,2));

for i = 1:1:size(lambda,2)
    disp([num2str(i/size(lambda,2)*100, '%.0f') '%'])
    for j = 1:1:size(power,2)
        for k = 1:1:size(r,2)
            for l = 1:1:size(z,2)
                if menu_profile == 1 % Bessel
                    intensity(i,j,k,l) = (power(j)*1e-3 / (pi * (R*1e-6)^2)) / j12 * J02(k);
                elseif menu_profile == 2 % Gaussian
                    wz(i,l) = w0 * sqrt(1+ (lambda(i)/1e3*z(l)/pi/w0^2)^2); % um. Beam diameter at distance z from the end of the fibre
                    intensity(i,j,k,l) = 2*power(j)*1e-3 / (pi*(wz(i,l)*1e-6)^2) * exp(-2*r(k)^2/wz(i,l)^2);
                end
                for m = 1:1:size(d,2)
                    potential(i,j,k,l,m) = - 1/2 * eps0 * eta * nmed * real(alpha(m,i)) * intensity(i,j,k,l) / KE;
                end
            end
        end
    end
end
% potential(lambda,power,r,z,d)

%% Options
% *************************************************************************

options = {['Wavelength: ' num2str(min(lambda)) ' to ' num2str(max(lambda)) ' nm'], ...
           ['Power: ' num2str(min(power)) ' to ' num2str(max(power)) ' mW'], ...
           ['Radial coordinate r : ' num2str(min(r)) ' to ' num2str(max(r)) ' um'], ...
           ['Axial coordinate z: ' num2str(min(z)) ' to ' num2str(max(z)) ' um'], ...
           ['Particle diameter d: ' num2str(min(d)) ' to ' num2str(max(d)) ' nm'], ...
           };
       
parameters = {lambda, ...
              power, ...
              r, ...
              z, ...
              d, ...
              };

parameter_names = {'\lambda', ...
                   'P', ...
                   'r', ...
                   'z', ...
                   'd', ...
                   };
               
parameter_units = {'nm', ...
                   'mW', ...
                   '\mum', ...
                   '\mum', ...
                   'nm', ...
                   };

fixed_parameters = {num2str(min(lambda)), ...
                    num2str(min(power)), ...
                    num2str(min(r)), ...
                    num2str(min(z)), ...
                    num2str(min(d)), ...
                    };
               

%% Plotting
% *************************************************************************

x_index = 1;
y_index = 1;
while x_index == y_index
    [x_index, y_index] = dialog_two_lists('Select axis and legend for line plot:', ...
                                          'X Axis:', options, ...
                                          'Legend:', options);
end
fixed_parameter_indices = 1:1:size(options,2);
fixed_parameter_indices([x_index, y_index]) = [];

input_title = 'Select fixed parameter values';
input_data = options(fixed_parameter_indices);
limits = fixed_parameters(fixed_parameter_indices);
dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
answer = inputdlg(input_data,input_title,dim,limits,dlg_options);
fixed_parameters(fixed_parameter_indices) = answer;

indices = cell(size(options));
title_text = 'Au NP';
for i = 1:1:size(options,2)
    if i == x_index || i == y_index
        indices{i} = 1:1:size(potential,i);
    else
        [~,indices{i}] = min(abs(parameters{i}-str2double(fixed_parameters{i})));
        title_text = strcat(title_text, ...
            [', ' parameter_names{i} ' = ' num2str(parameters{i}(indices{i})) ' ' parameter_units{i}]);
    end
end

figure_line = figure('Units','normalized','Position',[0.1 0.1 0.8 0.7], 'tag', 'figure_line');

h = plot(parameters{x_index}, ...
         squeeze(potential(indices{1}, indices{2}, indices{3}, indices{4}, indices{5})), ...
         'LineWidth', 2); hold all
plot_legend = cell(size(h,1),1);
for i = indices{y_index}
    plot_legend{i} = [parameter_names{y_index} ' = ' num2str(parameters{y_index}(i)) ' ' parameter_units{y_index}];
end
xlabel([parameter_names{x_index} ' (' parameter_units{x_index} ')'])
ylabel('Energy (kT@294K)')
set(gca,'FontSize', 14)
legend(plot_legend, 'Location', 'EO')
title(title_text)
grid on

%% colour scheme
figure(figure_line)

colour_type = {'DEFAULT', ...
               'parula', 'jet', 'hsv', 'cool', ...
               'spring', 'summer', 'autumn', 'winter', ...
               'gray', 'copper',...
               'red', 'green', 'aqua', 'blue', 'purple',...
               };

selected_colour = 3;
[selected_colour, ~] = listdlg('PromptString', 'Colour scheme:',...
                           'SelectionMode', 'single', ...
                           'ListString', colour_type,...
                           'InitialValue', selected_colour);
for i = 1:1:size(h,1)
    if selected_colour > 1 
        colour_RGB = colour_gradient(i, size(h,1), colour_type(selected_colour));
        h(i).Color = colour_RGB;  
    end
    h(i).MarkerSize = 1;
    h(i).LineStyle = '-';
    h(i).LineWidth = 2;
end


%% Contour Plot
% *************************************************************************

x_index = 1;
y_index = 1;
while x_index == y_index
    [x_index, y_index] = dialog_two_lists('Select axis for 2D contour plot:', ...
                                          'X Axis:', options, ...
                                          'Y Axis:', options);
end
fixed_parameter_indices = 1:1:size(options,2);
fixed_parameter_indices([x_index, y_index]) = [];

input_title = 'Select fixed parameter values';
input_data = options(fixed_parameter_indices);
limits = fixed_parameters(fixed_parameter_indices);
dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
answer = inputdlg(input_data,input_title,dim,limits,dlg_options);
fixed_parameters(fixed_parameter_indices) = answer;

indices = cell(size(options));
title_text = 'Au NP';
for i = 1:1:size(options,2)
    if i == x_index || i == y_index
        indices{i} = 1:1:size(potential,i);
    else
        [~,indices{i}] = min(abs(parameters{i}-str2double(fixed_parameters{i})));
        title_text = strcat(title_text, ...
            [', ' parameter_names{i} ' = ' num2str(parameters{i}(indices{i})) ' ' parameter_units{i}]);
    end
end

contour_values = squeeze(potential(indices{1}, indices{2}, indices{3}, indices{4}, indices{5}));
if size(contour_values,1) == size(parameters{x_index},2)
    contour_values = contour_values';
end

figure('Units','normalized','Position',[0.1 0.08 0.8 0.8]);
contour_levels = linspace(min(min(contour_values)), max(max(contour_values)), 100);
contourf(parameters{x_index},...
         parameters{y_index},...
         contour_values,...
         'LineStyle', 'none',...
         'LevelListMode', 'manual', ...
         'LevelList', contour_levels);
colormap(flipud(jet))
colorbar
xlabel([parameter_names{x_index} ' (' parameter_units{x_index} ')'])
ylabel([parameter_names{y_index} ' (' parameter_units{y_index} ')'])
title({'Energy (kT@294K)', title_text})
set(gca, 'FontSize', 16);
