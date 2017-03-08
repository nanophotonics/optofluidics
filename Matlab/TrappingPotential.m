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

d.name = 'Diameter';
d.units = 'nm';
d.min = 8; 
d.max = 10; 
d.step = 2; 
d.values = d.min : d.step : d.max;
d.size = size(d.values,2);

lambda.name = 'Wavelength';
lambda.units = 'nm';
lambda.min = 725;
lambda.max = 795;
lambda.step = 5;
lambda.values = lambda.min : lambda.step : lambda.max;
lambda.size = size(lambda.values,2);

power.name = 'Power';
power.units = 'mW';
power.min = 100;
power.max = 500;
power.step = 100;
power.values = power.min : power.step : power.max;
power.size = size(power.values,2);

parameters = [d,lambda,power];

column_titles = {'Min:', 'Max:', 'Step:'};
variable_limits = zeros(size(parameters,2),size(column_titles,2));
row_titles = cell(size(parameters,2),1);
count = 0;
for i = parameters
    count = count + 1;
    row_titles{count} = [i.name ' (' i.units ')'];
    variable_limits(count,:) = [i.min, i.max, i.step];
end
          
% variable_limits = dialog_table(row_titles,column_titles,variable_limits);

count = 0;
for i = parameters
    count = count + 1;
    i.min = variable_limits(count,1);
    i.max = variable_limits(count,2);
    i.step = variable_limits(count,3);
    i.values = i.min : i.step : i.max;
    i.size = size(i.values,2);
    parameters(count) = i;
    disp(i)
end

d = parameters(1);
lambda = parameters(2);
power = parameters(3);


%% Read Gold refractive index file and interpolate
% *************************************************************************
foldername = '';
filename = 'Gold Refractive-Rakic.txt';
n_data = dlmread([foldername filename], '\t', 1, 0);
% n_data(:,1) = wave (nm)
% n_data(:,2) = n
% n_data(:,3) = k
n_real = spline(n_data(:,1), n_data(:,2), lambda.values);
n_imaginary = spline(n_data(:,1), n_data(:,3), lambda.values);
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
alpha_0 = zeros(d.size,lambda.size);
for i = 1:1:d.size
    alpha_0(i,:) = 4*pi*(d.values(i)*1e-9/2)^3 .* ...
        (n_particle.^2 - nmed^2) ./ (n_particle.^2 + 2*nmed^2); 
end

% polarisability (SI units)
alpha = alpha_0;

%% Calculate the intensity and potential
% *************************************************************************

KE = kb*T;  % J, kinetic energy

menu_profile = 1;
% menu_profile = menu('Beam Profile?', 'Bessel', 'Gaussian');
    
if menu_profile == 1 % Bessel
    r.name = 'Radial coordinate';
    r.units = 'um';
    r.min = -R; 
    r.max = R; 
    r.size = 200;
    r.values = linspace(r.min,r.max,r.size);
    r.step = r.values(2)-r.values(1);     

%     r = linspace(-R, R, 200); % um. Radial coordinate inside fibre 
%     r = 0; % um. Radial coordinate inside fibre

    j12 = 0.269514; % Square of the Bessel function 1st kind 1st order at the first zero of the 0th order
    j01 = 2.40483; % First zero of the Bessel function 1st kind 0th order
    J02 = besselj(0,j01*r.values/R).^2; % Square of the Bessel function 1st kind 0th order

%     z = 0;
%     z = 0:1:5;
    z.name = 'Axial coordinate';
    z.units = 'um';
    z.min = 0; 
    z.max = 0; 
    z.size = 200;
    z.values = linspace(z.min,z.max,z.size);
    z.step = z.values(2)-z.values(1);
    
elseif menu_profile == 2 % Gaussian
    r = linspace(0, 15, 100); % um. Radial coordinate
    z = 0:20:200; % um. Distance from the end of the fibre
    w0 = 5; % um. Beam radius at the end of the fibre
    % MAKE POP UP WINDOW
    wz = zeros(lambda.size,z.size);
end

intensity = zeros(lambda.size,power.size,r.size,z.size,d.size);
potential = zeros(lambda.size,power.size,r.size,z.size,d.size);

for i = 1:1:lambda.size
    clc, disp([num2str(i/lambda.size*100, '%.0f') '%'])
    for j = 1:1:power.size
        for k = 1:1:r.size
            for l = 1:1:z.size
                if menu_profile == 1 % Bessel
                    intensity(i,j,k,l) = (power.values(j)*1e-3 / (pi * (R*1e-6)^2)) / j12 * J02(k); % mW / nm^2
                elseif menu_profile == 2 % Gaussian
                    wz(i,l) = w0 * sqrt(1+ (lambda.values(i)/1e3*z.values(l)/pi/w0^2)^2); % um. Beam diameter at distance z from the end of the fibre
                    intensity(i,j,k,l) = 2*power.values(j)*1e-3 / (pi*(wz(i,l)*1e-6)^2) * exp(-2*r.values(k)^2/wz(i,l)^2);
                end
                for m = 1:1:d.size
                    potential(i,j,k,l,m) = - 1/2 * eps0 * eta * nmed * real(alpha(m,i)) * intensity(i,j,k,l) / KE;
                end
            end
        end
    end
end
% intensity(lambda,power,r,z,d)
% potential(lambda,power,r,z,d)

%% Calculate the gradient force, velocity, and time
% *************************************************************************
F_grad = zeros(size(intensity));
for i = 1:1:lambda.size
    for j = 1:1:power.size
        for l = 1:1:z.size
            for m = 1:1:d.size
                F_grad(i,j,:,l,m) = 1/2 * eps0 * eta * nmed * real(alpha(m,i)) * ...
                    gradient(squeeze(intensity(i,j,:,l,m)),r.values);
            end
        end
    end
end

velocity = zeros(size(F_grad));
for m = 1:1:d.size
    velocity(:,:,:,:,m) = F_grad(:,:,:,:,m) ./ (6*pi*eta*d(m)/2);
end

% time = zeros(size(velocity));
% for i = 1:1:lambda.size
%     for j = 1:1:power.size
%         for k = 1:1:r.size
%             time(i,j,:) = trapz(r(1:1:k),velocity(i,j,1:1:k));
%         end
%     end
% end


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
               

%% Line Plot
% *************************************************************************

menu_plot = 1;
menu_plot = menu('Which variable to plot?', ...
    'Potential', ...
    'Intensity', ...
    'Gradient Force',...
    'Velocity');

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
variable_limits = fixed_parameters(fixed_parameter_indices);
dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
answer = inputdlg(input_data,input_title,dim,variable_limits,dlg_options);
fixed_parameters(fixed_parameter_indices) = answer;

indices = cell(size(options));
if menu_profile == 1 % bessel
    title_text = 'Au NP';
elseif menu_profile == 2 % gaussian
    title_text = ['Au NP, w0 = ' num2str(w0) ' \mum'];
end
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

if menu_plot == 1 % potential
    h = plot(parameters{x_index}, ...
             squeeze(potential(indices{1}, indices{2}, indices{3}, indices{4}, indices{5})), ...
             'LineWidth', 2); hold all
elseif menu_plot == 2 % intensity
    h = plot(parameters{x_index}, ...
             squeeze(intensity(indices{1}, indices{2}, indices{3}, indices{4})), ...
             'LineWidth', 2); hold all
elseif menu_plot == 3 % gradient force
    h = plot(parameters{x_index}, ...
             squeeze(F_grad(indices{1}, indices{2}, indices{3}, indices{4}, indices{5})), ...
             'LineWidth', 2); hold all
elseif menu_plot == 4 % velocity
    h = plot(parameters{x_index}, ...
             squeeze(velocity(indices{1}, indices{2}, indices{3}, indices{4}, indices{5})), ...
             'LineWidth', 2); hold all
end
plot_legend = cell(size(h,1),1);
for i = indices{y_index}
    plot_legend{i} = [parameter_names{y_index} ' = ' num2str(parameters{y_index}(i)) ' ' parameter_units{y_index}];
end
xlabel([parameter_names{x_index} ' (' parameter_units{x_index} ')'])
if menu_plot == 1 % potential
    ylabel(['Energy (kT@' num2str(T) 'K)'])
elseif menu_plot == 2 % intensity
    ylabel('Intensity (mW/nm^2)')
elseif menu_plot == 3 % gradient force
    ylabel('F_{grad} (???)')
elseif menu_plot == 4 % velocity
    ylabel('Velocity (???)')
end

set(gca,'FontSize', 14)
legend(plot_legend, 'Location', 'EO')
title(title_text)
grid on

% colour scheme
figure(figure_line)

colour_type = {'DEFAULT', ...
               'parula', 'jet', 'hsv', 'cool', ...
               'spring', 'summer', 'autumn', 'winter', ...
               'gray', 'copper',...
               'red', 'green', 'aqua', 'blue', 'purple',...
               };

selected_colour = 3;
% [selected_colour, ~] = listdlg('PromptString', 'Colour scheme:',...
%                            'SelectionMode', 'single', ...
%                            'ListString', colour_type,...
%                            'InitialValue', selected_colour);
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

menu_plot = 1;
menu_plot = menu('Which variable to plot?', ...
    'Potential', ...
    'Intensity', ...
    'Gradient Force',...
    'Velocity');

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
variable_limits = fixed_parameters(fixed_parameter_indices);
dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
answer = inputdlg(input_data,input_title,dim,variable_limits,dlg_options);
fixed_parameters(fixed_parameter_indices) = answer;

indices = cell(size(options));
if menu_profile == 1 % bessel
    title_text = 'Au NP';
elseif menu_profile == 2 % gaussian
    title_text = ['Au NP, w0 = ' num2str(w0) ' \mum'];
end
for i = 1:1:size(options,2)
    if i == x_index || i == y_index
        indices{i} = 1:1:size(potential,i);
    else
        [~,indices{i}] = min(abs(parameters{i}-str2double(fixed_parameters{i})));
        title_text = strcat(title_text, ...
            [', ' parameter_names{i} ' = ' num2str(parameters{i}(indices{i})) ' ' parameter_units{i}]);
    end
end

if menu_plot == 1 % potential
    contour_values = squeeze(potential(indices{1}, indices{2}, indices{3}, indices{4}, indices{5}));
elseif menu_plot == 2 % intensity
    contour_values = squeeze(intensity(indices{1}, indices{2}, indices{3}, indices{4}));
elseif menu_plot == 3 % gradient force
    contour_values = squeeze(F_grad(indices{1}, indices{2}, indices{3}, indices{4}));
end

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
if menu_plot == 1 % potential
    title({['Energy (kT@' num2str(T) 'K)'], title_text})
elseif menu_plot == 2 % intensity
    title({'Intensity (mW/nm^2)', title_text}) 
elseif menu_plot == 3 % gradient force
    title({'F_{grad} (???)', title_text}) 
end
set(gca, 'FontSize', 16);
