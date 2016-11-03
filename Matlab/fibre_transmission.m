% Created by Ana Andres-Arroyo (aa938)
% Calculates the in-coupled laser power based on 
% the Ti:S calibration and fibre loss.

function transmission = fibre_transmission(medium, wavelength, fibre_length)
% medium = 'H2O' or 'D2O: must be a string
% wavelength in nm: can be a single number or an array of numbers
% fibre_length in cm: can be a single number or an array numbers

% input_title = 'Parameters'; 
% input_data = {'Medium:', ...
%               'Laser Wavelength (nm):', ...
%               'Fibre Length (cm):', ...           
%               };
% default_values = {medium,...
%                   num2str(wavelength),...
%                   num2str(fibre_length),...
%                   };
% dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 60];
% answer = inputdlg(input_data, input_title, dim, default_values, dlg_options);
% medium = answer{1};   
% wavelength = str2double(answer{2});   
% fibre_length = str2double(answer{3});   

fibre_attenuation = absorption(medium, wavelength); % dB/m. 
% Just includes the H2O/D2O absorption and no other losses in the fibre!!!
% The total fibre attenuation should be measured experimentally with the cutback method.  

transmission = zeros(max(size(wavelength)), max(size(fibre_attenuation)));
for i = 1:1:max(size(fibre_length))
    transmission(:,i) = (10.^(- fibre_attenuation.*(fibre_length(i)/100)/10))*100; % percentage
end

end