clc
clear 
close all
figures = {};

%% Initial parameters
% *************************************************************************
eta = 377; % Ohm. Characteristic impedance of free space
eps0 = 8.85e-12; % F/m. Permittivity of free space
kb = 1.38e-23; % m^2*kg*s^-2*K^-1. Boltzmann constant

R = 9; % um. Radius of the HCPCF core
nmed = 1.33; % Refractive index of water
T = 273 + 70; % K. Ambient temperature

d.name = 'Diameter';
d.symbol = 'd';
d.units = 'nm';
d.min = 8; 
d.max = 10; 
d.step = 2; 
d.values = d.min : d.step : d.max;
d.size = size(d.values,2);
d.selected = d.min;

lambda.name = 'Wavelength';
lambda.units = 'nm';
lambda.symbol = '\lambda';
lambda.min = 725;
lambda.max = 975;
lambda.step = 50;
lambda.values = lambda.min : lambda.step : lambda.max;
lambda.size = size(lambda.values,2);
lambda.selected = lambda.min;

power.name = 'Power';
power.symbol = 'P';
power.units = 'mW';
power.min = 100;
power.max = 500;
power.step = 100;
power.values = power.min : power.step : power.max;
power.size = size(power.values,2);
power.selected = power.min;

parameters = [d,lambda,power];

column_titles = {'Min:', 'Max:', 'Step:'};
default_values = zeros(size(parameters,2),size(column_titles,2));
row_titles = cell(size(parameters,2),1);
for i = 1:1:size(parameters,2)
    row_titles{i} = [parameters(i).name ' (' parameters(i).units ')'];
    default_values(i,:) = [parameters(i).min, parameters(i).max, parameters(i).step];
end
          
default_values = dialog_table(row_titles,column_titles,default_values);

i_counter = 0;
for i = parameters
    i_counter = i_counter + 1;
    i.min = default_values(i_counter,1);
    i.max = default_values(i_counter,2);
    i.step = default_values(i_counter,3);
    i.values = i.min : i.step : i.max;
    i.size = size(i.values,2);
    parameters(i_counter) = i;
%     disp(i)
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
for m = 1:1:d.size
    alpha_0(m,:) = 4*pi*(d.values(m)*1e-9/2)^3 .* ...
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
    r.symbol = 'r';
    r.units = '\mum';
    r.min = -R; 
    r.max = R; 
    r.size = 20;
    r.values = linspace(r.min,r.max,r.size);
    r.step = r.values(2)-r.values(1);     
    r.selected = r.min;
    
    j12 = 0.269514; % Square of the Bessel function 1st kind 1st order at the first zero of the 0th order
    j01 = 2.40483; % First zero of the Bessel function 1st kind 0th order
    J02 = besselj(0,j01*r.values/R).^2; % Square of the Bessel function 1st kind 0th order

    z.name = 'Axial coordinate';
    z.symbol = 'z';
    z.units = 'um';
    z.min = 0; 
    z.max = 0; 
    z.size = 1;
    z.values = linspace(z.min,z.max,z.size);
    z.selected = z.min;
    if z.size > 1
        z.step = z.values(2)-z.values(1);
    else
        z.step = 0;
    end
    
elseif menu_profile == 2 % Gaussian
    % NEEDS UPDATING TO THE NEW STRUCTURE FORMAT!!!
    r = linspace(0, 15, 100); % um. Radial coordinate
    z = 0:20:200; % um. Distance from the end of the fibre
    w0 = 5; % um. Beam radius at the end of the fibre
    % MAKE POP UP WINDOW
    wz = zeros(lambda.size,z.size);
end

intensity.name = 'Intensity';
intensity.symbol = 'I';
intensity.units = 'mW/nm^2';
intensity.values = zeros(lambda.size,power.size,r.size,z.size,d.size);
    
potential.name = 'Potential Energy';
potential.symbol = 'U';
potential.units = ['kT@' num2str(T) 'K'];
potential.values = zeros(lambda.size,power.size,r.size,z.size,d.size);

% intensity = zeros(lambda.size,power.size,r.size,z.size,d.size);
% potential = zeros(lambda.size,power.size,r.size,z.size,d.size);


for i = 1:1:lambda.size;
    clc, disp([num2str(i/lambda.size*100, '%.0f') '%'])    
    for j = 1:1:power.size        
        for k = 1:1:r.size            
            for l = 1:1:z.size                
                if menu_profile == 1 % Bessel
                    intensity.values(i,j,k,l) = ...
                        (power.values(j)*1e-3 / (pi * (R*1e-6)^2)) / j12 * J02(k); % mW / nm^2
                elseif menu_profile == 2 % Gaussian
                    % NEEDS UPDATING TO THE NEW STRUCTURE FORMAT!!!
                    wz(i,l) = w0 * sqrt(1+ (lambda.values(i)/1e3*z.values(l)/pi/w0^2)^2); % um. Beam diameter at distance z from the end of the fibre
                    intensity.values(i,j,k,l) = 2*power.values(j)*1e-3 / (pi*(wz(i,l)*1e-6)^2) * exp(-2*r.values(k)^2/wz(i,l)^2);
                end
                for m = 1:1:d.size
                    potential.values(i,j,k,l,m) = ...
                        - 1/2 * eps0 * eta * nmed * real(alpha(m,i)) * ...
                        intensity.values(i,j,k,l) / KE;
                end
            end
        end
    end
end
% disp([i,j,k,l,m])
% intensity(lambda,power,r,z,d)
% potential(lambda,power,r,z,d)

%% Calculate the gradient force, velocity, and time
% *************************************************************************
F_grad.name = 'Gradient Force';
F_grad.symbol = 'F_{grad}';
F_grad.units = '???';
F_grad.values = zeros(size(intensity.values));

for i = 1:1:lambda.size
    for j = 1:1:power.size
        for l = 1:1:z.size
            for m = 1:1:d.size
                F_grad.values(i,j,:,l,m) = 1/2 * eps0 * eta * nmed * real(alpha(m,i)) * ...
                    gradient(squeeze(intensity.values(i,j,:,l,m)),r.values);
                % IT's NOT WORKING WELL
            end
        end
    end
end

velocity.name = 'Particle velocity';
velocity.symbol = 'v';
velocity.units = '???';
velocity.values = zeros(size(F_grad.values));

for m = 1:1:d.size
    velocity.values(:,:,:,:,m) = F_grad.values(:,:,:,:,m) ./ (6*pi*eta*d.values(m)/2);
end

time.name = 'Response time';
time.symbol = 'v';
time.units = '???';
time.values = zeros(size(velocity.values));

for i = 1:1:lambda.size
    for j = 1:1:power.size
        for k = 1:1:r.size
            for l = 1:1:z.size
                for m = 1:1:d.size
%                     time.values(i,j,k,l,m) = trapz(r.values(1:1:k),velocity.values(i,j,1:1:k,l,m));
                    time.values(i,j,k,l,m) = sum(velocity.values(i,j,1:1:k,l,m));
                    % UPDATE THIS FORMULA TO INTEGRATE FROM r = 0 to r = r.values(k)
                end
            end
        end
    end
end


%% Options
% *************************************************************************

parameters = [lambda,power,r,z,d];
parameters_options = cell(size(parameters));
i_counter = 0;
for i = parameters
    i_counter = i_counter + 1;
    parameters_options{i_counter} = [i.name ': ' i.symbol ': ' ...
        num2str(i.min) ' to ' num2str(i.max) ' ' i.units];
end
x_index = 3;
y_index = 2;

variables = [potential, intensity, F_grad, velocity, time];
variables_options = cell(size(variables));
for i = 1:1:size(variables,2)
    variables_options{i} = variables(i).name;
end
menu_variable = 3;

plot_type = {'Line','Contour'};
menu_plot = 1;


%% Plot
% *************************************************************************

[menu_plot, menu_variable] = dialog_two_lists('Select plot options:', ...
                                      'Type', plot_type, menu_plot,...
                                      'Variable', variables_options, menu_variable);

if menu_plot == 1 % line plot
    plot_options{1} = 'X Axis';
    plot_options{2} = 'Legend';
elseif menu_plot == 2 % contour plot
    plot_options{1} = 'X Axis';
    plot_options{2} = 'Y Axis';
end

okay = 0;
while okay == 0
    [x_index, y_index] = dialog_two_lists('Select plot options:', ...
                                          plot_options{1}, parameters_options, x_index,...
                                          plot_options{2}, parameters_options, y_index);
    if x_index == y_index
        waitfor(errordlg([plot_options{1} ' and ' plot_options{2} ' must be different'], 'Error'));
    else
        okay = 1;
    end
end

fixed_parameter_indices = 1:1:size(parameters_options,2);
fixed_parameter_indices([x_index, y_index]) = [];

default_values = cell(size(fixed_parameter_indices));
i_counter = 0;
for i = fixed_parameter_indices
    i_counter = i_counter + 1;
    default_values{i_counter} = num2str(parameters(i).selected);
end

input_title = 'Select fixed parameter values';
input_data = parameters_options(fixed_parameter_indices);
dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
answer = inputdlg(input_data,input_title,dim,default_values,dlg_options);

i_counter = 0;
for i = fixed_parameter_indices
    i_counter = i_counter + 1;
    parameters(i).selected = str2double(answer{i_counter});
end

indices = cell(size(parameters_options));
if menu_profile == 1 % bessel
    title_text = 'Au NP';
elseif menu_profile == 2 % gaussian
    title_text = ['Au NP, w0 = ' num2str(w0) ' \mum'];
end

for i = 1:1:size(parameters,2)  
    if i == x_index || i == y_index
        indices{i} = 1:1:size(variables(1).values,i);
    else
        [~,indices{i}] = min(abs(parameters(i).values - parameters(i).selected));
        parameters(i).selected = parameters(i).values(indices{i});
        title_text = strcat(title_text, [', ' parameters(i).symbol ' = ' ...
            num2str(parameters(i).selected) ' ' parameters(i).units]);
    end
end


figures{end+1} = figure('Units','normalized','Position',[0.1 0.1 0.8 0.7], 'tag', 'figure_line');

if menu_plot == 1 % line plot
    h = plot(parameters(x_index).values, ...
             squeeze(variables(menu_variable).values(indices{1}, indices{2}, indices{3}, indices{4}, indices{5})), ...
             'LineWidth', 2); hold all

    plot_legend = cell(size(h,1),1);
    for i = indices{y_index}
        plot_legend{i} = [parameters(y_index).symbol ' = ' num2str(parameters(y_index).values(i)) ' ' parameters(y_index).units];
    end
    xlabel([parameters(x_index).name ': ' parameters(x_index).symbol ' (' parameters(x_index).units ')'])
    ylabel([variables(menu_variable).name ': ' variables(menu_variable).symbol ' (' variables(menu_variable).units ')'])

    set(gca,'FontSize', 14)
    legend(plot_legend, 'Location', 'EO')
    title(title_text)
    grid on
    
    % colour scheme
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

elseif menu_plot == 2 % contour plot
    contour_values = squeeze(variables(menu_variable).values(indices{1}, indices{2}, indices{3}, indices{4}, indices{5}));
    
    if size(contour_values,1) == parameters(x_index).size
        contour_values = contour_values';
    end
    
    contour_levels = linspace(min(min(contour_values)), max(max(contour_values)), 100);
    
    contourf(parameters(x_index).values,...
             parameters(y_index).values,...
             contour_values,...
             'LineStyle', 'none',...
             'LevelListMode', 'manual', ...
             'LevelList', contour_levels);
    colormap(flipud(jet))
    colorbar
    
    xlabel([parameters(x_index).name ': ' parameters(x_index).symbol ' (' parameters(x_index).units ')'])
    ylabel([parameters(y_index).name ': ' parameters(y_index).symbol ' (' parameters(y_index).units ')'])
    title({title_text,...
        [variables(menu_variable).name ': ' variables(menu_variable).symbol ' (' variables(menu_variable).units ')']})
    set(gca, 'FontSize', 16);
end
