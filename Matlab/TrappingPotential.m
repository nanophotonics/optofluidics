% Written by Ana Andres-Arroyo (aa938)
% Calculates the optical potential of a Bessel/Gaussian beam
% Calculates the gradient force and particle vGrad
% Fits the gradient force to Hooke's law for small displacements

clc
clear 
close all
figures = {};

directory_save = 'R:\aa938\NanoPhotonics\Laboratory\2017.03.23 - gradient force calculations\';

%% Initial parameters
% *************************************************************************
eps0 = 8.85e-12; % F/m. Permittivity of free space
kb = 1.38e-23; % m^2*kg/s^2/K^1. Boltzmann constant
c = 3e8; % m/s. Speed of light
T = 273.15 + 20; % K. Ambient temperature
viscosity_log = (1.3272*(293.15-T)-0.001053*(T-293.15).^2) ./ (T-168.15) - 2.999; % for T > 293.15 K = 20 C
viscosity = 10.^viscosity_log; % Pa*s. Viscosity of water

units = 2;
units = menu('Units:', 'SI', 'CGS');
% units_SI = units_display * units_conversion

d.name = 'Diameter';
d.symbol = 'd';
d.units_SI = 'm';
d.units_display = 'nm';
d.units_conversion = 1e-9;
d.min = 80e-9; 
d.max = 80e-9; 
d.step = 1e-9; 
d.values = d.min : d.step : d.max;
d.size = size(d.values,2);
d.selected = 8e-9;

lambda.name = 'Wavelength';
lambda.symbol = '\lambda';
lambda.units_SI = 'm';
lambda.units_display = 'nm';
lambda.units_conversion = 1e-9;
lambda.min = 400e-9;
lambda.max = 720e-9;
lambda.step = 2e-9;
lambda.values = lambda.min : lambda.step : lambda.max;
lambda.size = size(lambda.values,2);
lambda.selected = 800e-9;

power.name = 'Power';
power.symbol = 'P';
power.units_SI = 'W';
power.units_display = 'mW';
power.units_conversion = 1e-3; 
power.min = 1000e-3;
power.max = 1000e-3;
power.step = 100e-3;
power.values = power.min : power.step : power.max;
power.size = size(power.values,2);
power.selected = power.max;

parameters = [d,lambda,power];

column_titles = {'Min:', 'Max:', 'Step:'};
limits = zeros(size(parameters,2),size(column_titles,2));
row_titles = cell(size(parameters,2),1);
for i = 1:1:size(parameters,2)
    row_titles{i} = [parameters(i).name ' (' parameters(i).units_display ')'];
    limits(i,:) = [parameters(i).min, parameters(i).max, parameters(i).step] ...
        / parameters(i).units_conversion;
end
          
limits = dialog_table(row_titles,column_titles,limits);

i_counter = 0;
for i = parameters
    i_counter = i_counter + 1;
    i.min = limits(i_counter,1) * i.units_conversion;
    i.max = limits(i_counter,2) * i.units_conversion;
    i.step = limits(i_counter,3) * i.units_conversion;
    i.values = i.min : i.step : i.max;
    i.size = size(i.values,2);
    parameters(i_counter) = i;
end

d = parameters(1);
lambda = parameters(2);
power = parameters(3);


%% READ: refractive index
% *************************************************************************
nmed = 1.33; % Refractive index of water

foldername = '';
filename = 'Gold Refractive-Rakic.txt';
n_data = dlmread([foldername filename], '\t', 1, 0);
% n_data(:,1) = wave (nm)
% n_data(:,2) = n
% n_data(:,3) = k
n_real = spline(n_data(:,1), n_data(:,2), lambda.values * 1e9);
n_imaginary = spline(n_data(:,1), n_data(:,3), lambda.values * 1e9);
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

%% CALCULATE: polarisability
% *************************************************************************

% standard polarisability
alpha_0 = zeros(d.size,lambda.size);
for m = 1:1:d.size
    if units == 1 % SI units: F.m^2
        alpha_0(m,:) = nmed^2*eps0*4*pi*((d.values(m))/2)^3 .* ...
            (n_particle.^2 - nmed^2) ./ (n_particle.^2 + 2*nmed^2); 
    elseif units == 2 % CGS units: m^3
        alpha_0(m,:) = 4*pi*(d.values(m)/2)^3 .* ...
            (n_particle.^2 - nmed^2) ./ (n_particle.^2 + 2*nmed^2);
    end
    
end

% polarisability
% a correction to the polarisability may sometimes be needed
alpha = alpha_0;

%% CALCULATE: intensity, potential
% *************************************************************************

menu_profile = 2;
menu_profile = menu('Beam Profile?', 'Bessel', 'Gaussian');

r.name = 'Radial coordinate';
r.symbol = 'r';
r.units_SI = 'm';
r.units_display = '\mum';
r.units_conversion = 1e-6;

z.name = 'Axial coordinate';
z.symbol = 'z';
z.units_SI = 'm';
z.units_display = '\mum';
z.units_conversion = 1e-6;

if menu_profile == 1 % Bessel
    R = 9e-6; % m. Radius of the HCPCF core
    input_title = 'Bessel beam parameters.';
    input_data = 'Radius of HCPCF core: R (um):';
    default_values = {num2str(R*1e6)};
    dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
    answer = inputdlg(input_data,input_title,dim,default_values,dlg_options);
    R = str2double(answer{1})*1e-6; % m
    
    r.min = -R; 
    r.max = R; 
    r.size = 200; % select and even number
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
%     w0 = 0.35e-6; % m. Beam radius
    w0 = 5.3e-6/2; % m. Beam radius
       
    input_title = 'Gaussian beam parameters.';
    input_data = 'Beam radius: w0 (um):';
    default_values = {num2str(w0*1e6)};
    dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
    answer = inputdlg(input_data,input_title,dim,default_values,dlg_options);
    w0 = str2double(answer{1})*1e-6; % m
    
    r.min = -w0*2; 
    r.max = w0*2;    
    r.size = 101; % select and odd number
    r.values = linspace(r.min,r.max,r.size);
    r.step = r.values(2)-r.values(1);     
    r.selected = r.min;

    z.min = 0; 
    z.max = 2*w0*10; % m
    z.size = 21; % select and odd number
    z.values = linspace(z.min,z.max,z.size);
    z.step = z.values(2)-z.values(1);
    z.selected = z.min;    
    
    wz = zeros(lambda.size,z.size);
end

intensity.name = 'Beam Intensity';
intensity.symbol = 'I';
intensity.units_SI = 'W/m^2';
intensity.units_display = 'W/m^2';
intensity.units_conversion = 1;
intensity.values = zeros(lambda.size,power.size,r.size,z.size,d.size);
    
potential.name = 'Potential Energy';
potential.symbol = 'U';
potential.units_SI = 'J';
potential.units_display = ['kT@' num2str(T) 'K'];
potential.units_conversion = kb*T; % kinetic energy
% potential.units_display = 'J';
% potential.units_conversion = 1; % kinetic energy
potential.values = zeros(lambda.size,power.size,r.size,z.size,d.size);

for i = 1:1:lambda.size;
    clc, disp([num2str(i/lambda.size*100, '%.0f') '%: intensity and potential'])    
    for j = 1:1:power.size        
        for k = 1:1:r.size            
            for l = 1:1:z.size                
                if menu_profile == 1 % Bessel
                    intensity.values(i,j,k,l,:) = ...
                        (power.values(j) / (pi * (R)^2)) / j12 * J02(k);
                elseif menu_profile == 2 % Gaussian
                    wz(i,l) = w0 * sqrt(1+ ...
                        (lambda.values(i)*z.values(l)/pi/w0^2)^2); % m. Beam diameter at distance z 
                    intensity.values(i,j,k,l,:) = 2*power.values(j) / ...
                        (pi*wz(i,l)^2) * exp(-2*r.values(k)^2/wz(i,l)^2);
                end
                for m = 1:1:d.size
                    if units == 1 % SI
                    potential.values(i,j,k,l,m) = ...
                        -1/2*(2/(c*nmed*eps0))*real(alpha(m,i)) * ...
                    intensity.values(i,j,k,l,m); 
%                     intensity.values(i,j,k,l,m)/KE; 
                    elseif units == 2 % CGS
                    end                     
                                                
                end
            end
        end
    end
end
% disp(size(intensity.values))
% intensity(lambda,power,r,z,d)
% potential(lambda,power,r,z,d)

%% CALCULATE: gradient force
% *************************************************************************
Fgrad.name = 'Gradient Force';
Fgrad.symbol = 'F_{grad}';
Fgrad.units_SI = 'N'; 
Fgrad.units_display = 'pN';
Fgrad.units_conversion = 1e-12;
Fgrad.values = zeros(size(intensity.values));

for i = 1:1:lambda.size
    clc, disp([num2str(i/lambda.size*100, '%.0f') '%: gradient force'])    
    for j = 1:1:power.size
        for l = 1:1:z.size
            for m = 1:1:d.size
                if units == 1 % SI
                    Fgrad.values(i,j,:,l,m) = (real(alpha(m,i))/(2*c*nmed*eps0))*...
                        gradient(squeeze(intensity.values(i,j,:,l,m)),r.values);        
                elseif units == 2 % CGS
                    Fgrad.values(i,j,:,l,m) = (real(alpha(m,i) * nmed)/(2*c)) * ...
                        gradient(squeeze(intensity.values(i,j,:,l,m)),r.values);
                end
            end
        end
    end
end

%% CALCULATE: scattering force
% *************************************************************************
Cext.name = 'Extinction Cross-Section';
Cext.symbol = 'C_{ext}';
Cext.units_SI = 'm^2'; 
Cext.units_display = 'm^2';
Cext.units_conversion = 1;
Cext.values = zeros(lambda.size,d.size);

Qext.name = 'Extinction Efficiency';
Qext.symbol = 'Q_{ext}';
Qext.units_SI = ''; 
Qext.units_display = '';
Qext.units_conversion = 1;
Qext.values = zeros(lambda.size,d.size);

Fscat.name = 'Scattering Force';
Fscat.symbol = 'F_{scat}';
Fscat.units_SI = 'N'; 
Fscat.units_display = 'pN';
Fscat.units_conversion = 1e-12;
Fscat.values = zeros(size(intensity.values));

% zeros(lambda.size,power.size,r.size,z.size,d.size);
for i = 1:1:lambda.size
% for i = 176
    clc, disp([num2str(i/lambda.size*100, '%.0f') '%: scattering'])    
    for m = 1:1:d.size
%         lambda.values(i)
%         Qext = 6; % read it from the paper
%         Cext(i) = Qext * (pi*(d.values(m)/2)^2);
%         Cext(i) = 2*pi/lambda.values(i) / eps0 / nmed^2 * imag(alpha(m,i)); % m^2, SI
        if units == 1 % alpha SI units: F.m^2
%             Cext.values(i,m) = 2*pi/lambda.values(i) / eps0 / nmed^2 * imag(alpha(m,i)); % m^2, SI
            Cext.values(i,m) = 2*pi/lambda.values(i) / eps0 * nmed^2 * imag(alpha(m,i)); % m^2, SI
        elseif units == 2 % alpha in CGS units: m^3
            Cext.values(i,m) = 2*pi/lambda.values(i) * imag(alpha(m,i)); % m^2, SI
        end
        Qext.values(i,m) = Cext.values(i,m) / (pi*(d.values(m)/2)^2);
        
        for j = 1:1:power.size
            for l = 1:1:z.size
            
%                 Fscat.values(i,j,:,l,m) = 8*pi*nmed * (2*pi/lambda.values(i))^4 * ...
%                     (d.values(m)/2)^6 / c * ...
%                     imag((n_particle(i)^2 - nmed^2) / (n_particle(i)^2 + 2*nmed^2)) *...
%                     intensity.values(i,j,:,l,m); 
                
%                 Fscat.values(i,j,:,l,m) = nmed * (2*pi/lambda.values(i))^4 /6/pi/c * ...
%                     (alpha(m,i)*conj(alpha(m,i))) * intensity.values(i,j,:,l,m); 
%                     % check that this is correct
%                     % CGS units
                    
%                 Cext(i) = 5.7e-15; % from Jain paper with Qext = 15, and reff = 11nm;
                
                Fscat.values(i,j,:,l,m) = nmed / c * intensity.values(i,j,:,l,m) * Cext.values(i,m); 
                    % check that this is correct
                    
            end
        end
    end
end

% close all
% figure
% plot(lambda.values*1e9, Cext), hold all
% plot(lambda.values*1e9, Qext.values), hold all

% alpha_0 = zeros(d.size,lambda.size);
% for m = 1:1:d.size
%     if units == 1 % SI units: F.m^2
%         alpha_0(m,:) = nmed^2*eps0*4*pi*((d.values(m))/2)^3 .* ...
%             (n_particle.^2 - nmed^2) ./ (n_particle.^2 + 2*nmed^2); 
%     elseif units == 2 % CGS units: m^3
%         alpha_0(m,:) = 4*pi*(d.values(m)/2)^3 .* ...
%             (n_particle.^2 - nmed^2) ./ (n_particle.^2 + 2*nmed^2);
%     end
%     
% end

%% CALCULATE: velocity, time
% *************************************************************************

vGrad.name = 'Particle velocity (due to Fgrad)';
vGrad.symbol = 'v';
vGrad.units_SI = 'm/s';
vGrad.units_display = '\mum/ms';
vGrad.units_conversion = 1e-6/1e-3;
vGrad.values = zeros(size(Fgrad.values));

for m = 1:1:d.size
    vGrad.values(:,:,:,:,m) = Fgrad.values(:,:,:,:,m) ./ ...
        (6*pi*viscosity*d.values(m)/2);
end

vScat.name = 'Particle velocity (due to Fscat)';
vScat.symbol = 'v';
vScat.units_SI = 'm/s';
vScat.units_display = '\mum/ms';
vScat.units_conversion = 1e-6/1e-3;
vScat.values = zeros(size(Fscat.values));

for m = 1:1:d.size
    vScat.values(:,:,:,:,m) = Fscat.values(:,:,:,:,m) ./ ...
        (6*pi*viscosity*d.values(m)/2);
end

time.name = 'Response time';
time.symbol = 't';
time.units_SI = 's';
time.units_display = 'ms';
time.units_conversion = 1e-3;
time.values = zeros(size(vGrad.values));

% MAKE SURE THERE IS AN EVEN NUMBER OF r.values
% [~,zero_index] = min(abs(r.values));
% for i = 1:1:lambda.size
%     clc, disp([num2str(i/lambda.size*100, '%.0f') '%: time'])    
%     for j = 1:1:power.size
%         for k = 1:1:r.size
%             if k ~= zero_index
%                 if k < zero_index
%                     integration_indices = k:1:zero_index;
%                 elseif k > zero_index
%                     integration_indices = zero_index:1:k;
%                 end
%                 for l = 1:1:z.size
%                     for m = 1:1:d.size
%                         time.values(i,j,k,l,m) = abs(trapz(r.values(integration_indices),...
%                             squeeze(1./vGrad.values(i,j,integration_indices,l,m))));
%                     end
%                 end
%             end
%         end
%     end
% end

%% LINEAR FIT: trap stiffness
% *************************************************************************

r_linear = r.max * 0.2;
[~, linear_indices] = intersect(find(r.values < r_linear), find(r.values > - r_linear));

kr.name = 'Radial Trap Stiffness';
kr.symbol = 'kr';
kr.units_SI = 'N/m/W'; 
kr.units_display = 'pN/\mum/mW';
kr.units_conversion = 1e-12/1e-6/1e-3;
kr.values = zeros(lambda.size,z.size,d.size);

% zeros(lambda.size,power.size,r.size,z.size,d.size);
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
            kr.values(i,l,m) = -fitting_parameters(1) / power.values(j);
            
%             if i == 1
%                 plot(r.values, ...
%                     squeeze(Fgrad.values(i,j,:,l,m)), ...
%                     '-', 'LineWidth', 2), hold all   
%                 plot_legend{end+1} = [d.name ' = ' num2str(d.values(m)) ' ' d.units_SI];
%                 plot(r.values(linear_indices), ...
%                     r.values(linear_indices)* fitting_parameters(1) + fitting_parameters(2), ...
%                     '--', 'LineWidth', 2), hold all      
%                 plot_legend{end+1} = ['kr = ' num2str(kr.values(i,l,m)) ' ' kr.units_SI];
%                 legend(plot_legend, 'Location', 'EO')
%                 xlabel([r.name ': ' r.symbol ' (' r.units_SI ')'])
%                 ylabel([Fgrad.name ': ' Fgrad.symbol ' (' Fgrad.units_SI ')'])
%                 grid on
%                 set(gca, 'FontSize', 16);     
%             end
        end
    end
end

%% Default Variables
% *************************************************************************
variables = [potential, intensity, Fgrad, kr, vGrad, time, Fscat, vScat, Cext, Qext];
variables_options = cell(size(variables));
for i = 1:1:size(variables,2)
    variables_options{i} = variables(i).name;
end
selected_variable = 2;

plot_styles = {'Line','Contour'};
selected_style = 2;
if strcmp(variables(selected_variable).symbol, 'kr')                        
    x_index = 3;
    y_index = 1;
elseif strcmp(variables(selected_variable).symbol, 'C_{ext}') % doesn't work
    x_index = 1;
    y_index = 2;
elseif strcmp(variables(selected_variable).symbol, 'Q_{ext}') % doesn't work
    x_index = 1;
    y_index = 2;
else
%     x_index = 3;
%     y_index = 4;
    x_index = 1;
    y_index = 2;
end

%% Variable Units
% *************************************************************************

Fgrad.units_display = 'fN';
Fgrad.units_conversion = 1e-15;

Fscat.units_display = 'pN';
Fscat.units_conversion = 1e-12;

time.units_display = 'ms';
time.units_conversion = 1e-3;

vGrad.units_display = 'm/s';
vGrad.units_conversion = 1;


%% Plot Options
% *************************************************************************
variables = [potential, intensity, Fgrad, kr, vGrad, time, Fscat, vScat, Cext, Qext];

[selected_style, selected_variable] = dialog_two_lists('Select plot options:', ...
                                      'Style', plot_styles, selected_style,...
                                      'Variable', variables_options, selected_variable);
if strcmp(variables(selected_variable).symbol, 'kr')                        
    parameters = [lambda,z,d];
elseif strcmp(variables(selected_variable).symbol, 'C_{ext}') 
    parameters = [lambda,d];
elseif strcmp(variables(selected_variable).symbol, 'Q_{ext}')  
    parameters = [lambda,d];
else
    parameters = [lambda,power,r,z,d];
end


parameters_options = cell(size(parameters));
i_counter = 0;
for i = parameters
    i_counter = i_counter + 1;
    parameters_options{i_counter} = [i.name ': ' i.symbol ': ' ...
        num2str(i.min / i.units_conversion) ' to ' ...
        num2str(i.max / i.units_conversion) ' ' i.units_display];
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
    default_values{i_counter} = num2str(parameters(i).selected / parameters(i).units_conversion);
end

if isempty(fixed_parameter_indices)
    selected_style = 1; % line plot
else
    input_title = 'Select fixed parameter values';
    input_data = parameters_options(fixed_parameter_indices);
    dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
    answer = inputdlg(input_data,input_title,dim,default_values,dlg_options);
    i_counter = 0;
    for i = fixed_parameter_indices
        i_counter = i_counter + 1;
        parameters(i).selected = str2double(answer{i_counter}) * parameters(i).units_conversion;
    end
end

if strcmp(variables(selected_variable).symbol, 'kr')                        
    lambda = parameters(1);
    z = parameters(2);
    d = parameters(3);
elseif strcmp(variables(selected_variable).symbol, 'C_{ext}') 
    lambda = parameters(1);
    d = parameters(2);
elseif strcmp(variables(selected_variable).symbol, 'Q_{ext}') 
    lambda = parameters(1);
    d = parameters(2);
else
    lambda = parameters(1);
    power = parameters(2);
    r = parameters(3);
    z = parameters(4);
    d = parameters(5);
end

indices = cell(size(parameters_options));
if menu_profile == 1 % bessel
    title_text = ['Au NP, Bessel, R = ' num2str(R*1e6) ' \mum'];
elseif menu_profile == 2 % gaussian
    title_text = ['Au NP, Gaussian, w0 = ' num2str(w0*1e6) ' \mum'];
end

for i = 1:1:size(parameters,2)  
    if i == x_index || i == y_index
        indices{i} = 1:1:parameters(i).size;
    else
        [~,indices{i}] = min(abs(parameters(i).values - parameters(i).selected));
        parameters(i).selected = parameters(i).values(indices{i});
        title_text = strcat(title_text, [', ' parameters(i).symbol ' = ' ...
            num2str(parameters(i).selected / parameters(i).units_conversion) ' ' ...
            parameters(i).units_display]);
    end
end


% %% Plot
% *************************************************************************

% selected_style = menu('Plot Style: ', plot_styles);
%  [selected_style, selected_variable] = dialog_two_lists('Select plot options:', ...
%                                        'Style', plot_styles, selected_style,...
%                                        'Variable', variables_options, selected_variable);

figures{end+1} = figure('Units','normalized','Position',[0.1 0.1 0.8 0.7]);

if selected_style == 1 % line plot
    for i = 1:1:size(indices,2)
    end
    h = plot(parameters(x_index).values / parameters(x_index).units_conversion, ...
             squeeze(variables(selected_variable).values(indices{:}) / variables(selected_variable).units_conversion), ...
             'LineWidth', 2); hold all

    plot_legend = cell(size(h,1),1);
    for i = indices{y_index}
        plot_legend{i} = [parameters(y_index).symbol ' = ' ...
            num2str(parameters(y_index).values(i) / parameters(y_index).units_conversion) ' ' ...
            parameters(y_index).units_display];
    end
    xlabel([parameters(x_index).name ': ' parameters(x_index).symbol ' (' parameters(x_index).units_display ')'])
    ylabel([variables(selected_variable).name ': ' variables(selected_variable).symbol ' (' variables(selected_variable).units_display ')'])

    set(gca,'FontSize', 14)
    legend(plot_legend, 'Location', 'EO', 'FontSize', 10)
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
    contour_values = squeeze(variables(selected_variable).values(indices{:})/ ...
        variables(selected_variable).units_conversion);
    
    if size(contour_values,1) == parameters(x_index).size
        contour_values = contour_values';
    end
    
    contour_levels = linspace(min(min(contour_values)), max(max(contour_values)), 100);
    
    contourf(parameters(x_index).values / parameters(x_index).units_conversion,...
             parameters(y_index).values / parameters(y_index).units_conversion,...
             contour_values,...
             'LineStyle', 'none',...
             'LevelListMode', 'manual', ...
             'LevelList', contour_levels);
%     colormap(flipud(jet))
    colormap(jet)
    colorbar
    
    xlabel([parameters(x_index).name ': ' parameters(x_index).symbol ' (' parameters(x_index).units_display ')'])
    ylabel([parameters(y_index).name ': ' parameters(y_index).symbol ' (' parameters(y_index).units_display ')'])
    title({title_text,...
        [variables(selected_variable).name ': ' variables(selected_variable).symbol ' (' variables(selected_variable).units_display ')']})
    set(gca, 'FontSize', 16);
end

%% Saving figures
% *************************************************************************
menu_save_figures = 1;
menu_save_figures = menu('Save Figures?', 'NO', 'YES');
if menu_save_figures == 2    
    for i = 1:1:max(size(figures))
        if findobj(figures{i}) ~= 0
            figure_save = figures{i};
            file_name_save = '';
            
            figure(figure_save)
            pause(0.1)
            [file_name_save,directory_save,~] = uiputfile(['.' 'png'],...
                'File to Save the Figure',[directory_save file_name_save]);
            hgexport(figure_save, [directory_save file_name_save], hgexport('factorystyle'), 'Format', 'png')
            file_name_save = strrep(file_name_save, 'png', 'fig');    
            saveas(figure_save, [directory_save file_name_save], 'fig');

        end
    end
end