% Written by Ana Andres-Arroyo (aa938)
% Calculates the optical potential of a Bessel/Gaussian beam
% Calculates the gradient force and particle velocity
% Fits the gradient force to Hooke's law for small displacements

clc
clear 
close all
figures = {};

%% Initial parameters
% *************************************************************************
eta = 377; % Ohm. Characteristic impedance of free space
eps0 = 8.85e-12; % F/m. Permittivity of free space
kb = 1.38e-23; % m^2*kg*s^-2*K^-1. Boltzmann constant
c=3E8; %light speed
nmed = 1.33; % Refractive index of water
T = 273.15 + 20; % K. Ambient temperature
viscosity_log = (1.3272*(293.15-T)-0.001053*(T-293.15).^2) ./ (T-168.15) - 2.999; % for T > 293.15 K = 20 C
viscosity = 10.^viscosity_log; % Pa*s

d.name = 'Diameter';
d.symbol = 'd';
d.units = 'nm';
d.min = 60; 
d.max = 60; 
d.step = 10; 
d.values = d.min : d.step : d.max;
d.size = size(d.values,2);
d.selected = 8;

lambda.name = 'Wavelength';
lambda.units = 'nm';
lambda.symbol = '\lambda';
lambda.min = 800;
lambda.max = 800;
lambda.step = 50;
lambda.values = lambda.min : lambda.step : lambda.max;
lambda.size = size(lambda.values,2);
lambda.selected = 800;

power.name = 'Power';
power.symbol = 'P';
power.units = 'mW';
power.min = 1000;
power.max = 1000;
power.step = 50;
power.values = power.min : power.step : power.max;
power.size = size(power.values,2);
power.selected = power.max;

parameters = [d,lambda,power];

column_titles = {'Min:', 'Max:', 'Step:'};
limits = zeros(size(parameters,2),size(column_titles,2));
row_titles = cell(size(parameters,2),1);
for i = 1:1:size(parameters,2)
    row_titles{i} = [parameters(i).name ' (' parameters(i).units ')'];
    limits(i,:) = [parameters(i).min, parameters(i).max, parameters(i).step];
end
          
limits = dialog_table(row_titles,column_titles,limits);

i_counter = 0;
for i = parameters
    i_counter = i_counter + 1;
    i.min = limits(i_counter,1);
    i.max = limits(i_counter,2);
    i.step = limits(i_counter,3);
    i.values = i.min : i.step : i.max;
    i.size = size(i.values,2);
    parameters(i_counter) = i;
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

% standard polarisability (CGS units: cm^3) ?????
alpha_0 = zeros(d.size,lambda.size);
for m = 1:1:d.size
    alpha_0(m,:) = nmed^2*eps0*4*pi*((d.values(m)*1e-9)/2)^3 .* ...
        (n_particle.^2 - nmed^2) ./ (n_particle.^2 + 2*nmed^2); %from  https://www.photonics.ethz.ch/fileadmin/user_upload/FILES2SHARE/ignatovitch02.pdf
end

% polarisability (CGS units: cm^3) ?????
% a correction to the polarisability may sometimes be needed
alpha = alpha_0;

%% CALCULATE: intensity, potential
% *************************************************************************

KE = kb*T;  % J, kinetic energy

menu_profile = 2;
menu_profile = menu('Beam Profile?', 'Bessel', 'Gaussian');

r.name = 'Radial coordinate';
r.symbol = 'r';
r.units = '\mum';

z.name = 'Axial coordinate';
z.symbol = 'z';
z.units = '\mum';    

if menu_profile == 1 % Bessel
    R = 9; % um. Radius of the HCPCF core
    input_title = 'Bessel beam parameters.';
    input_data = 'Radius of HCPCF core: R (um):';
    default_values = {num2str(R)};
    dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
    answer = inputdlg(input_data,input_title,dim,default_values,dlg_options);
    R = str2double(answer{1});
    
    r.min = -R; 
    r.max = R; 
    r.size = 201; % select and odd number
    r.values = linspace(r.min,r.max,r.size);
    r.step = r.values(2)-r.values(1);     
    r.selected = r.max;
    
    j12 = 0.269514; % Square of the Bessel function 1st kind 1st order at the first zero of the 0th order
    j01 = 2.40483; % First zero of the Bessel function 1st kind 0th order
    J02 = besselj(0,j01*r.values/R).^2; % Square of the Bessel function 1st kind 0th order

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
    w0 = 1; % um. Beam radius at the end of the fibre
       
    input_title = 'Gaussian beam parameters.';
    input_data = 'Beam waist: w0 (um):';
    default_values = {num2str(w0)};
    dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
    answer = inputdlg(input_data,input_title,dim,default_values,dlg_options);
    w0 = str2double(answer{1});
    
    r.min = -w0*3; 
    r.max = w0*3;    
    r.size = 101; % select and odd number
    r.values = linspace(r.min,r.max,r.size);
    r.step = r.values(2)-r.values(1);     
    r.selected = r.min;

    z.min = 0; 
    z.max = 20; 
    z.size = 201; % select and odd number
    z.values = linspace(z.min,z.max,z.size);
    z.step = z.values(2)-z.values(1);
    z.selected = z.min;    
    
    wz = zeros(lambda.size,z.size);
end

intensity.name = 'Beam Intensity';
intensity.symbol = 'I';
intensity.units = 'W/m^2'; %changed to W/m^2
intensity.values = zeros(lambda.size,power.size,r.size,z.size,d.size);
    
potential.name = 'Potential Energy';
potential.symbol = 'U';
potential.units = ['kT@' num2str(T) 'K'];
potential.values = zeros(lambda.size,power.size,r.size,z.size,d.size);

for i = 1:1:lambda.size;
    clc, disp([num2str(i/lambda.size*100, '%.0f') '%: intensity and potential'])    
    for j = 1:1:power.size        
        for k = 1:1:r.size            
            for l = 1:1:z.size                
                if menu_profile == 1 % Bessel
                    intensity.values(i,j,k,l,:) = ...
                        (power.values(j)*1e-3 / (pi * (R*1e-6)^2)) / j12 * J02(k); % W / m^2
                elseif menu_profile == 2 % Gaussian
                    % NEEDS UPDATING TO THE NEW STRUCTURE FORMAT!!!
                    wz(i,l) = w0 * sqrt(1+ (lambda.values(i)/1e3*z.values(l)/pi/w0^2)^2); % um. Beam diameter at distance z from the end of the fibre
                    intensity.values(i,j,k,l,:) = 2*power.values(j)*1e-3 / (pi*(wz(i,l)*1e-6)^2) * exp(-2*r.values(k)^2/wz(i,l)^2); %W/m^2
                end
                for m = 1:1:d.size
                    potential.values(i,j,k,l,m) = ...
                        -1/2*(2/(c*nmed*eps0))*real(alpha(m,i)) * ...
                    intensity.values(i,j,k,l,m)/KE; 
                     
                                                
                end
            end
        end
    end
end
% disp(size(intensity.values))
% intensity(lambda,power,r,z,d)
% potential(lambda,power,r,z,d)

%% CALCULATE: gradient force, velocity, time
% *************************************************************************
Fgrad.name = 'Gradient Force';
Fgrad.symbol = 'F_{grad}';
Fgrad.units = 'N'; 
Fgrad.values = zeros(size(intensity.values));

for i = 1:1:lambda.size
    clc, disp([num2str(i/lambda.size*100, '%.0f') '%: gradient force'])    
    for j = 1:1:power.size
        for l = 1:1:z.size
            for m = 1:1:d.size
                % polarisability in cm^3
                Fgrad.values(i,j,:,l,m) = 1/2*(2/(c*nmed*eps0))*real(alpha(m,i)) * ...
                    gradient(squeeze(intensity.values(i,j,:,l,m)),r.values*1E-6);        
            end
        end
    end
end

velocity.name = 'Particle Velocity';
velocity.symbol = 'v';
velocity.units = 'm/s';
velocity.values = zeros(size(Fgrad.values));

for m = 1:1:d.size
    velocity.values(:,:,:,:,m) = Fgrad.values(:,:,:,:,m) ./ (6*pi*viscosity*d.values(m)/2);
end

time.name = 'Response time';
time.symbol = 't';
time.units = 's';
time.values = zeros(size(velocity.values));

% MAKE SURE THERE IS AN ODD NUMBER OF r.values
[~,zero_index] = min(abs(r.values));
for i = 1:1:lambda.size
    clc, disp([num2str(i/lambda.size*100, '%.0f') '%: time'])    
    for j = 1:1:power.size
        for k = 1:1:r.size
            if k ~= zero_index
                if k < zero_index
                    integration_indices = k:1:zero_index;
                elseif k > zero_index
                    integration_indices = zero_index:1:k;
                end
                for l = 1:1:z.size
                    for m = 1:1:d.size
                        time.values(i,j,k,l,m) = abs(trapz(r.values(integration_indices),...
                            squeeze(velocity.values(i,j,integration_indices,l,m))));
                    end
                end
            end
        end
    end
end

%% LINEAR FIT: trap stiffness
% *************************************************************************

r_linear = r.max * 0.2;
[~, linear_indices] = intersect(find(r.values < r_linear), find(r.values > - r_linear));

kr.name = 'Radial Trap Stiffness';
kr.symbol = 'kr';
kr.units = 'pN/\mum/mW'; 
kr.values = zeros(lambda.size,z.size,d.size);



% Fgrad.values = zeros(lambda.size,power.size,r.size,z.size,d.size);
j = power.size;
for i = 1:1:lambda.size
    clc, disp([num2str(i/lambda.size*100, '%.0f') '%: trap stiffness'])    
%     if i == 1
%         figures{end+1} = figure('Units','normalized','Position',[0.1 0.1 0.8 0.7]);
%         plot_legend = {};
%     end
    for l = 1:1:z.size
        for m = 1:1:d.size
            fitting_parameters = polyfit(r.values(linear_indices), ...
                squeeze(Fgrad.values(i,j,linear_indices,l,m))', 1);
            kr.values(i,l,m) = -fitting_parameters(1) / power.values(j) * 1e3;
            
%             if i == 1
%                 plot(r.values, ...
%                     squeeze(Fgrad.values(i,j,:,l,m)), ...
%                     '-', 'LineWidth', 2), hold all   
%                 plot_legend{end+1} = [d.name ' = ' num2str(d.values(m)) ' ' d.units];
%                 plot(r.values(linear_indices), ...
%                     r.values(linear_indices)* fitting_parameters(1) + fitting_parameters(2), ...
%                     '--', 'LineWidth', 2), hold all      
%                 plot_legend{end+1} = ['kr = ' num2str(kr.values(i,l,m)) ' pN/\mum/mW'];
%                 legend(plot_legend, 'Location', 'EO')
%                 xlabel([r.name ': ' r.symbol ' (' r.units ')'])
%                 ylabel([Fgrad.name ': ' Fgrad.symbol ' (' Fgrad.units ')'])
%                 grid on
%                 set(gca, 'FontSize', 16);     
%             end
        end
    end
end

%% Options
% *************************************************************************

variables = [potential, intensity, Fgrad, kr, velocity, time];
variables_options = cell(size(variables));
for i = 1:1:size(variables,2)
    variables_options{i} = variables(i).name;
end
selected_variable = 2;

plot_styles = {'Line','Contour'};
selected_style = 1;

[selected_style, selected_variable] = dialog_two_lists('Select plot options:', ...
                                      'Style', plot_styles, selected_style,...
                                      'Variable', variables_options, selected_variable);

if strcmp(variables(selected_variable).symbol, 'kr')                        
    parameters = [lambda,z,d];
    x_index = 3;
    y_index = 1;
else
    parameters = [lambda,power,r,z,d];
    x_index = 3;
    y_index = 2;
end

parameters_options = cell(size(parameters));
i_counter = 0;
for i = parameters
    i_counter = i_counter + 1;
    parameters_options{i_counter} = [i.name ': ' i.symbol ': ' ...
        num2str(i.min) ' to ' num2str(i.max) ' ' i.units];
end

if selected_style == 1 % line plot
    plot_options{1} = 'X Axis';
    plot_options{2} = 'Legend';
elseif selected_style == 2 % contour plot
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

fixed_parameter_indices = 1:1:size(parameters,2);
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
    title_text = ['Au NP, Bessel, R = ' num2str(R) ' \mum'];
elseif menu_profile == 2 % gaussian
    title_text = ['Au NP, Gaussian, w0 = ' num2str(w0) ' \mum'];
end

for i = 1:1:size(parameters,2)  
    if i == x_index || i == y_index
        indices{i} = 1:1:parameters(i).size;
    else
        [~,indices{i}] = min(abs(parameters(i).values - parameters(i).selected));
        parameters(i).selected = parameters(i).values(indices{i});
        title_text = strcat(title_text, [', ' parameters(i).symbol ' = ' ...
            num2str(parameters(i).selected) ' ' parameters(i).units]);
    end
end


%% Plot
% *************************************************************************

% selected_style = menu('Plot Style: ', plot_styles);
 [selected_style, selected_variable] = dialog_two_lists('Select plot options:', ...
                                       'Style', plot_styles, selected_style,...
                                       'Variable', variables_options, selected_variable);

figures{end+1} = figure('Units','normalized','Position',[0.1 0.1 0.8 0.7]);

if selected_style == 1 % line plot
    for i = 1:1:size(indices,2)
    end
    h = plot(parameters(x_index).values, ...
             squeeze(variables(selected_variable).values(indices{:})), ...
             'LineWidth', 2); hold all

    plot_legend = cell(size(h,1),1);
    for i = indices{y_index}
        plot_legend{i} = [parameters(y_index).symbol ' = ' num2str(parameters(y_index).values(i)) ' ' parameters(y_index).units];
    end
    xlabel([parameters(x_index).name ': ' parameters(x_index).symbol ' (' parameters(x_index).units ')'])
    ylabel([variables(selected_variable).name ': ' variables(selected_variable).symbol ' (' variables(selected_variable).units ')'])

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

elseif selected_style == 2 % contour plot
    contour_values = squeeze(variables(selected_variable).values(indices{:}));
    
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
        [variables(selected_variable).name ': ' variables(selected_variable).symbol ' (' variables(selected_variable).units ')']})
    set(gca, 'FontSize', 16);
end

