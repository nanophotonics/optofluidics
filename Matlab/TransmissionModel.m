clc
clear 
close all

L = 0.22; % Length of the fibre in m
sigma = 2.68*10^(-16); % extinction coefficient in m^2
rho_zero = (2.6*10^10)/(3*10^-6); %concentration in m^3 from bbi: 2.6*10^10 part/ml and since we are using 2:1 D20 :Au solution
% rho_array = linspace(0.01, 10, 10) * rho_zero;
rho_array = logspace(-1, 1, 15) * rho_zero;

R = 8*10^(-6); % Radius of the fibre in m
% wvalues = linspace(0.001*R,100*R,10000); %in m
% wvalues = linspace(0.100, 100, 1000) ./ 1e6; % in m
wvalues = logspace(-1, 2, 1000) ./ 1e6; % in m

I = @(r) besselj(0, 2.405*r./R).^2 ; % intensity
I0 = 1/(integral(@(r) I(r).*2*pi.*r,0,R)); % normalisation constant so that the input power (=Integral(I(r))dA over the area of the core=1

T = zeros(size(wvalues,2),size(rho_array,2));
k = 0;
for i = 1:1:size(wvalues,2)
    w = wvalues(i);
    for j = 1:1:size(rho_array,2)
        k = k + 1;
        clc, disp([num2str(k/(size(wvalues,2)*size(rho_array,2))*100) '%'])
        rho = rho_array(j);
        
        h = @(r) -L*sigma*((rho*R^2)/(w^2*(1-exp(-R^2/w^2)))*exp(-r.^2/w^2)); % intermediate function
        f = @(r) 2*pi*(I0)*I(r).*r.*exp(h(r)); % Function we integrate

        T(i,j) = integral(f, 0, R); % Power out
    end
end

%% plotting
close all
figure('Units','normalized','Position',[0.01 0.08 0.85 0.75]);
legend_cell = {};
for j = 1:1:size(rho_array,2)
    rho = rho_array(j);
%     semilogx(wvalues*1e6, T(:,j)*100), hold all %Power out = y axis , w = x axis
    plot_handle{j} = loglog(wvalues*1e6, T(:,j)*100); hold all %Power out = y axis , w = x axis
    legend_cell{end+1} = ['\rho / \rho_0 = ' num2str(rho/rho_zero)];
end
grid on
axis tight
ylim([0,100])
ylabel('Transmission (%)')
xlabel('Gaussian concentration width (\mum)')
title(['60 nm Au NP'...
    ' // Cext = ' num2str(sigma*1e12) ' \mum^2'...
    ' // \rho_0 = ' num2str(rho_zero/1e6,'%.2e') ' particles/mL'...
    ' // PBG HCF R = ' num2str(R*1e6) ' \mum, L = ' num2str(L*100) ' cm'...
    ])
legend(legend_cell, 'Location', 'SW')

% figure('Units','normalized','Position',[0.14 0.13 0.85 0.75]);
% contourf(wvalues*1e6,rho_array/rho_zero,T')

%% colour scheme
colour_type = {'DEFAULT', ...
    'parula', 'jet', 'hsv', 'cool', ...
    'spring', 'summer', 'autumn', 'winter', ...
    'gray', 'copper',...
    'red', 'green', 'aqua', 'blue', 'purple',...
    };
menu_colour = 2;
menu_colour = menu('Colour scheme', colour_type);
if menu_colour > 1 
    for j = 1:1:size(plot_handle,2)
        colour_RGB = colour_gradient(j, size(plot_handle,2), colour_type(menu_colour));
        plot_handle{j}.Color = colour_RGB;  
        plot_handle{j}.LineStyle = '-';
        plot_handle{j}.LineWidth = 2;
    end
end