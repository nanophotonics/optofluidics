clc
clear
close all

%% Parameters

reference_angle = 25; % degrees
reference_input_power = 11.5; % mW. Before the objective
reference_output_power = 4.3; % mW
fibre_length = 20; % cm. PBG HCF
wavelength = 810; % nm. Ti:Sa laser
medium = 'H2O';

input_title = 'Parameters'; 
input_data = {'Reference Ange (deg):',...
              'Reference Input Power (mW):', ...
              'Reference Output Power (mW):', ...
              'Fibre Length (cm):', ...
              'Laser Wavelength (nm):', ...
              'Medium:', ...
              };
resize = 'on'; dimensions = [1 60];
default_values = {num2str(reference_angle),...
                  num2str(reference_input_power),...
                  num2str(reference_output_power),...
                  num2str(fibre_length),...
                  num2str(wavelength),...
                  medium,...
                  };
answer = inputdlg(input_data, input_title, dimensions, default_values, resize);
reference_angle = str2double(answer{1});   
reference_input_power = str2double(answer{2});   
reference_output_power = str2double(answer{3});   
fibre_length = str2double(answer{4});   
wavelength = str2double(answer{5});   
medium = answer{6};   
    
fibre_attenuation = absorption(medium, wavelength); % dB/m. Just includes the H2O/D2O absorption.
reference_in_coupled_power = reference_output_power * 10^ (fibre_attenuation * (fibre_length / 100) / 10); % mW
in_coupling_transmission = reference_in_coupled_power / reference_input_power * 100; % percentage
fibre_transmission = reference_output_power / reference_in_coupled_power * 100; % percentage
