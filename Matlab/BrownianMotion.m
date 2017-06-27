clear
clc
close all
figures = {};

%% Equations
clc

kb = 1.38064852e-23; % m^2*kg/s^2/K

T = 25 + 273; % K
disp(['T = ' num2str(T,'%.4g') ' K'])

% eta = 5e-4; % Pa*s
% eta = 8.872e-4; % Pa*s
% eta = 9.55e-4; % Pa*s
% eta = 10e-4; % Pa*s
eta_log = (1.3272*(293.15-T)-0.001053*(T-293.15).^2) ./ (T-168.15) - 2.999; % for T > 293.15 K = 20 C
eta = 10.^eta_log; % Pa*s. Viscosity of water
disp(['eta = ' num2str(eta,'%.2e') ' Pa*s'])

% L = 58e-9; % m
% R = 3e-9/2; % m

L = 50.8e-9; % m
R = 12.5e-9/2; % m

% L = 50e-9; % m
% R = 3e-9/2; % m

disp(['L = ' num2str(L*1e9,'%.3g') ' nm'])
disp(['R = ' num2str(R*1e9,'%.3g') ' nm'])

% diffusion coefficient of a nanorod
g = 0.45; % depends on the detailed shape of the cylindrical rod
disp(['g = ' num2str(g,'%.3g')])

D_parallel = kb*T/2/pi/eta/L * (log(L/(R*2))-g);
disp(['D_parallel = ' num2str(D_parallel*1e12,'%.3g') ' um^2/s'])

D_perpendicular = kb*T/4/pi/eta/L * (log(L/(R*2))+g);
disp(['D_perpendicular = ' num2str(D_perpendicular*1e12,'%.3g') ' um^2/s'])

D_rotational = 3*kb*T/pi/eta/L^3 * (log(L/(R*2))-g);
% disp(['D_rotational = ' num2str(D_rotational*1e12,'%.3g') ' um^2/s'])
disp(['D_rotational = ' num2str(D_rotational,'%.3g') ' s^-1'])

% D_translational = kb*T/3/pi/eta/L * (log(2*L/(R*2))-g);
% disp(['D_translational = ' num2str(D_translational*1e12,'%.3g') ' um^2/s'])

% effective radius of a sphere with equivalent volume
% V = pi*R^2*L; % m^3
V = pi*R^2*(L-2*R)+4/3*pi*R^3; % m^3
r = (V*3/4/pi)^(1/3); % m
disp(['r = ' num2str(r*1e9,'%.3g') ' nm'])

% d = 1e-6; % m
d = 56e-9; % m
r = d/2; % m
% V = 4/3*pi*r^3; % m^3
disp(['r = ' num2str(r*1e9,'%.3g') ' nm'])

gamma = 6*pi*eta*r; % Stokes drag (SI units)
D = kb*T/gamma; % SI units: m^2/s
% disp(['D_Stokes = ' num2str(D,'%.3g') ' m^2/s'])
disp(['D_Stokes = ' num2str(D*1e12,'%.3g') ' um^2/s'])

%% other stuff

% disp(['sqrt(2*D/tau) = ' num2str(sqrt(2*D/tau),'%.3g') ' m/s'])

% ro = 19.3e3; % kg/m^3 GOLD
% ro = 2.65e3; % kg/m^3 SILICA
% m = ro*V; % kg

% tau = m/gamma;
% disp(['tau = ' num2str(tau*1e6,'%.3g') ' us'])

% c = 70e-6; % m/s
% T = c^2*m/2/kb; % K
% disp(['T = ' num2str(T,'%.3g') ' K'])

% c = sqrt(2*kb*T/m);
% disp(['c = ' num2str(c,'%.3g') ' m/s'])

% D = 11.3e-12;
% disp(['D_Stokes = ' num2str(D*1e12,'%.3g') ' um^2/s'])
% gamma = kb*T/D;
% r = gamma / 6/pi/eta;
% d = r*2;
% disp(['d = ' num2str(d*1e9,'%.3g') ' nm'])


%% Brownian motion simulation: initial parameters
% *************************************************************************

% code from: http://uk.mathworks.com/matlabcentral/fileexchange/32067-brownian-motion
% close all

D = 8e-12; % um^2/s
% tau = 0.00179e-6; % time interval (seconds)
% tau = 0.132e-6; % time interval (seconds)
frame_rate = 159.22; % fps
tau = 1/frame_rate; % time interval (seconds)
% tau = 1/140; % time interval (seconds)
% tau = 0.1; % time interval (seconds)
N = 1e2; % number of samples --> make it less than 1e7!!!!!!!!!
total_time = tau*N; % total time (seconds)
% total_time = 1; % total time (seconds)
% N = round(total_time/tau); % number of collisions
% h = sqrt(total_time/N); % scaling factor
h = sqrt(2*D*tau); % scaling factor
% h = 1; % scaling factor
np = 1e3; % Number of Particles
ndim = 2;

title_text = '';
title_text = [title_text 'tau = ' num2str(tau*1e3,'%.2g') ' ms'];
% title_text = [title_text ', t = ' num2str(total_time,'%.3g') ' s'];
title_text = [title_text ', N = ' num2str(N,'%.0f')];
title_text = [title_text ', D = ' num2str(D*1e12,'%.3g') ' \mum^2/s'];
% title_text = [title_text ', h = ' num2str(h,'%.2g')];
title_text = [title_text ', np = ' num2str(np,'%.0f')];
title_text = [title_text ', ' num2str(ndim) 'D'];
   

t_position = tau * (1:N);  % time vector for position (seconds)
t_velocity = t_position;

%% Calculating the velocity, displacement, and position
% *************************************************************************

% displacement = h*randn(N,np,ndim);
% velocity = displacement / tau;

velocity_brownian = h/tau * randn(N,np,ndim);

axis_labels = {'x','y','z'};
% velocity_axis = [0, 0, 0];
velocity_axis = [0, 8.4, 0]*1e-6;

% velocity_flow = ones(N,np,ndim);
% for i = 1:1:ndim
%     velocity_flow(:,:,i) = velocity_axis(i) * velocity_flow(:,:,i);
%     title_text = [title_text ', v' axis_labels{i} ' = ' ...
%         num2str(velocity_axis(i)*1e6, '%.2g') ' \mum/s'];
% end
% 
% velocity = velocity_brownian + velocity_flow;
% displacement = velocity * tau;
% position = cumsum(displacement);

% read velocity from figure
menu_vflow = 3;
% menu_vflow = menu('x-axis non-brownian velocity', 'Calculate', 'Read from line plot', 'Read from contour plot');

if menu_vflow > 1
    pathRead = 'R:\aa938\NanoPhotonics\Matlab\Optical Forces\2017.06.26\';
    nameRead = 'vscat-map.fig';
%     [nameRead, pathRead, ~] = uigetfile('.fig',...
%         'Select figure:',pathRead,'MultiSelect','off');

    figures{end+1} = open([pathRead nameRead]);
    hold all

    axes = get(figures{end},'Children');
    dataObjs = get(axes,'Children');
    objTypes = get(dataObjs{2}, 'Type');
    data = cell(size(objTypes,1),2);
      
    if menu_vflow == 2 % line
        selected_plot = 3;
        vflow_x = get(dataObjs{2}(selected_plot), 'XData');
        vflow_v = get(dataObjs{2}(selected_plot), 'YData');
        plot(vflow_x,vflow_v, '--k', 'LineWidth', 2)
        vflow_x = vflow_x*1e-6; % from um to m
        vflow_v = vflow_v*1e-3; % from mm/s to m/s
    elseif menu_vflow == 3 % contour
        selected_plot = 1;
        vflow_x = get(dataObjs{2}(selected_plot), 'XData');
        vflow_y = get(dataObjs{2}(selected_plot), 'YData');
        vflow_v = get(dataObjs{2}(selected_plot), 'ZData');
        vflow_x = vflow_x*1e-6; % from um to m
        vflow_y = vflow_y*1e-6; % from um to m
        vflow_v = vflow_v*1e-3; % from mm/s to m/s
        vflow_v = vflow_v * 0.1; % arbitrary scaling
    end
    close(figures{end})
end

velocity_flow = ones(N,np,ndim);
velocity = zeros(N,np,ndim);
displacement = zeros(N,np,ndim);
position = zeros(N+1,np,ndim);
% position(1,:,1) = position(1,:,1) + 10e-6;
position(1,:,1) = position(1,:,1) + rand(size(position(1,:,1)))*10e-6 + 5e-6;
position(1,:,2) = position(1,:,2) - (0.5-rand(size(position(1,:,2))))*10e-6;
for i = 1:1:N
    if menu_vflow == 1
        velocity_flow(i,:,1) = velocity_axis(1) * velocity_flow(i,:,1);
%         velocity_flow(i,:,1) = -30e-2*(80e-6 + position(i,:,1));
    elseif menu_vflow == 2 % line figure
        velocity_flow(i,:,1) = spline(vflow_x, vflow_v, position(i,:,1));
    elseif menu_vflow == 3 % contour figure
        velocity_flow(i,:,1) = interp2(vflow_x, vflow_y, vflow_v,...
            position(i,:,1), position(i,:,2));
    end    
    
    if ndim >= 2
        velocity_flow(i,:,2) = velocity_axis(2) * velocity_flow(i,:,2);
    end
    if ndim >= 3
        velocity_flow(i,:,3) = velocity_axis(3) * velocity_flow(i,:,3);
    end
    
    velocity(i,:,:) = velocity_brownian(i,:,:) + velocity_flow(i,:,:);
    displacement(i,:,:) = velocity(i,:,:) * tau;
    position(i+1,:,:) = position(i,:,:) + displacement(i,:,:);
    
%     title_text = [title_text ', v' axis_labels{i} ' = ' ...
%         num2str(velocity_axis(i)*1e6, '%.2g') ' \mum/s'];
end
position(end,:,:) = [];

displacement_squared = zeros(N,np);
position_squared = zeros(N,np);
for i = 1:1:ndim
    displacement_squared = displacement_squared + (displacement(:,:,i)).^2;
    position_squared = position_squared + (position(:,:,i)).^2;
end

%% Plot options
coordenate_labels = axis_labels(1:ndim);
if np == 1
    cmap_np = [1,0,0];
else
    cmap_np = parula(np);
end
cmap_xyz = parula(ndim);
plot_font_size = 14;

time.units = 's';
time.conversion = 1;
distance.units = '\mum';
distance.conversion = 1e-6;

plot_options = {};
plot_options{end+1} = '3D/2D Position'; 
% plot_options{end+1} = 'Position/Velocity vs. Time';
% plot_options{end+1} = 'Position Histogram';
% plot_options{end+1} = 'Position Squared vs.Time';
% plot_options{end+1} = 'Displacement Histogram';
plot_options{end+1} = 'Velocity Histogram';
plot_options{end+1} = 'Flow Velocity vs. Position';
po = 1:1:numel(plot_options);

% [po, ~] = listdlg('PromptString', 'Choose plots:',...
%     'SelectionMode', 'multiple', ...
%     'ListString', plot_options,...
%     'InitialValue', po);
    

%% Plotting the particle position
if find(strcmp(plot_options(po), '3D/2D Position'))
% if ndim > 1 && np <= 10 && N <= 1e5
    figures{end+1} = figure;
    if ndim == 3
        subplot(1,2,1)
        for k=1:np
            plot3(position(:,k,1)/distance.conversion,...
                position(:,k,2)/distance.conversion,...
                position(:,k,3)/distance.conversion,...
                'Color',cmap_np(k,:));
            hold on;
        end
        grid on
        set(gca,'FontSize', plot_font_size)
        xlabel([coordenate_labels{1} ' position (' distance.units ')'])
        ylabel([coordenate_labels{2} ' position (' distance.units ')'])
        zlabel([coordenate_labels{3} ' position (' distance.units ')'])
        title('3D position')
        % title(title_text)
        % figures{end+1} = figure;
        subplot(1,2,2)
    end
    for k=1:np
%         if ndim == 3
%             scatter(position(:,k,1),position(:,k,2),...
%                 (position(:,k,3)-min(min(position(:,:,3))))*...
%                 (max(max(position(:,k,3)))-min(min(position(:,:,3))))+1,...
%                 cmap_np(k,:)); hold on;
%         end
        plot(position(:,k,1)/distance.conversion,...
            position(:,k,2)/distance.conversion,...
            'Color',cmap_np(k,:)); hold on;
    end
    grid on
    set(gca,'FontSize', plot_font_size)
    axis square
%     xlim([-20,110])
%     ylim([-20,20])
%     title('Position top view')
    xlabel([coordenate_labels{1} ' position (' distance.units ')'])
    ylabel([coordenate_labels{2} ' position (' distance.units ')'])
    suptitle(title_text)
end

%% Plot particle position and velocity vs. time
if find(strcmp(plot_options(po), 'Position/Velocity vs. Time'))
    figures{end+1} = figure;
    position_time = cell(ndim,np);
    velocity_time = cell(ndim,np);
    j = 0;
    for i = 1:ndim
        % position
        j = j + 1;
        if np == 1
            subplot(1,2,1)
        else
            subplot(ndim,2,j)
        end

        for k=1:np
            position_time{i,k} = plot(t_position/time.conversion,...
                position(:,k,i)/distance.conversion); hold all
            position_time{i,k}.Color = cmap_np(k,:);
            position_time{i,k}.LineWidth = 2;
    %         xlim([0,600])
        end
        grid on
        set(gca,'FontSize', plot_font_size)
        xlabel('t (s)')
        ylabel([coordenate_labels{i} ' (' distance.units ')'])

        if np == 1
            position_time{i,k}.Color = cmap_xyz(i,:);
            ylabel(['position (' distance.units ')'])
            legend(coordenate_labels(1:i))
        end

        % velocity
        j = j + 1;
        if np == 1
            subplot(1,2,2)
        else
            subplot(ndim,2,j)
        end
        for k=1:np
            velocity_time{i,k} = plot(t_velocity/time.conversion,...
                velocity(:,k,i)/(distance.conversion/time.conversion)); hold all
            velocity_time{i,k}.Color = cmap_np(k,:);
            velocity_time{i,k}.LineWidth = 2;
    %         xlim([0,600])
        end
        grid on
        set(gca,'FontSize', plot_font_size)
        xlabel('t (s)')             
        ylabel(['v' coordenate_labels{i} ' (' distance.units '/' time.units ')'])

        if np == 1
            velocity_time{i,k}.Color = cmap_xyz(i,:);
            ylabel(['velocity (' distance.units '/' time.units ')'])
            legend(coordenate_labels(1:i))
        end
    end
    suptitle(title_text)
end

%% Plotting the position histograms
if find(strcmp(plot_options(po), 'Position Histogram'))
    figures{end+1} = figure;
    h_position = cell(ndim,1);
    for i = 1:ndim
        h_position{1} = histogram(position(:,:,i)); hold all
        h_position{i}.FaceColor = cmap_xyz(i,:);
        h_position{i}.Normalization = 'count';
        % h_position{i}.Normalization = 'probability';
    end
    legend(coordenate_labels)
    xlabel('position')
    suptitle(title_text)
    ylabel(h_position{1}.Normalization)
end

%% Plot particle position squared
if find(strcmp(plot_options(po), 'Position Squared vs.Time'))
    figures{end+1} = figure;
    % subplot(2,1,1)
    p_psq = cell(1,np);
    p_psq_fit = cell(1,np);
    legend_psq = {};
    text_psq = {};
    text_psq{end+1} = 'Linear fit:';
    text_psq{end+1} = 'y =  p1*x';
    % text_psq{end+1} = 'y =  p1*x + p2';
    text_psq{end+1} = '';

    linear_origin = fittype('p1*x');

    for k=1:np
        p_psq{k} = plot(t_position/time.conversion,...
            position_squared(:,k)/distance.conversion^2); hold all
        p_psq{k}.Color = cmap_np(k,:);
        p_psq{k}.LineWidth = 1;
        legend_psq{end+1} = ['particle ' num2str(k)];

    %     [fit_psq,gof_psq] = fit(t_position'/time.conversion, position_squared(:,k)/distance.conversion^2, 'poly1');
        [fit_psq,gof_psq] = fit(t_position'/time.conversion, ...
            position_squared(:,k)/distance.conversion^2, ...
            linear_origin,...
            'StartPoint', 2*ndim*D/(distance.conversion^2/time.conversion));
        confint_psq = confint(fit_psq);

        slope_error = abs((confint_psq(1,1)-confint_psq(2,1))/2/fit_psq.p1);

        p_psq_fit{k} = plot(fit_psq); hold all, 
        legend_psq{end+1} = 'linear fit';
        p_psq_fit{k}.Color = cmap_np(k,:);
        p_psq_fit{k}.LineWidth = 1;
        if np < 5
            text_psq{end+1} = [num2str(k) ': p1 = ' ...
                num2str(fit_psq.p1, '%.3g') ' \pm ' ...
                num2str(abs(confint_psq(1,1)-confint_psq(2,1))/2, '%.3g')];
        %     text_psq{end+1} = [num2str(k) ': p2 = ' ...
        %         num2str(fit_psq.p2, '%.3g') ' \pm ' ...
        %         num2str(abs(confint_psq(1,2)-confint_psq(2,2))/2, '%.3g')];
            text_psq{end+1} = [num2str(k) ': D = ' ...
                num2str(fit_psq.p1/2/ndim, '%.3g') ' \pm ' ...
                num2str(slope_error*fit_psq.p1/2/ndim, '%.3g') ...
                ' ' distance.units '^2/' time.units];
            text_psq{end+1} = [num2str(k) ': R^2 = ' ...
                num2str(gof_psq.rsquare, '%.3g')];
            text_psq{end+1} = '';
        end
    end

    if np > 1
        p_average = plot(t_position/time.conversion,...
            mean(position_squared,2)/distance.conversion^2); hold all
        p_average.LineWidth = 2;
        p_average.Color = 'r';
        p_average.LineStyle = '-';
        legend_psq{end+1} = 'average';

        [fit_psq,gof_psq] = fit(t_position'/time.conversion, ...
            mean(position_squared,2)/distance.conversion^2, ...
            linear_origin,...
            'StartPoint', 2*ndim*D);
        confint_psq = confint(fit_psq);

        slope_error = abs((confint_psq(1,1)-confint_psq(2,1))/2/fit_psq.p1);

        p_psq_fit{k+1} = plot(fit_psq); hold all, 
        legend_psq{end+1} = 'linear fit';
        p_psq_fit{k+1}.Color = 'r';
        p_psq_fit{k+1}.LineWidth = 2;
        text_psq{end+1} = ['average: p1 = ' ...
            num2str(fit_psq.p1, '%.3g') ' \pm ' ...
            num2str(abs(confint_psq(1,1)-confint_psq(2,1))/2, '%.3g')];
        text_psq{end+1} = ['average: D = ' ...
            num2str(fit_psq.p1/2/ndim, '%.3g') ' \pm ' ...
            num2str(slope_error*fit_psq.p1/2/ndim, '%.3g') ...
            ' ' distance.units '^2/' time.units];
        %     num2str(slope_error, '%.3g') '% m^2/t'];
        text_psq{end+1} = ['average: R^2 = ' ...
            num2str(gof_psq.rsquare, '%.3g')];
        text_psq{end+1} = '';
    end

    p_theory = plot(t_position/time.conversion, ...
        2*ndim*D/(distance.conversion^2/time.conversion)*t_position/time.conversion);
    p_theory.LineWidth = 2;
    p_theory.Color = 'k';
    p_theory.LineStyle = '--';
    legend_psq{end+1} = 'theory: 2*dim*D*t';

    grid on
    set(gca,'FontSize', plot_font_size)
    xlabel(['t (' time.units ')'])
    ylabel(['position squared (' distance.units '^2)'])
    suptitle(title_text)
    legend(legend_psq,'Location','EO')

    text('Units','normalized','Position',[0.08,0.95], ...
        'FontSize', 12, 'VerticalAlignment', 'top', 'String' , text_psq)
end

%% Plot particle displacement histogram
if find(strcmp(plot_options(po), 'Displacement Histogram'))
    figures{end+1} = figure;
    % subplot(2,1,2)
    h_displacement = cell(1,np);
    for k=1%:np
    %     h_displacement{k} = histogram(displacement_squared(:,k)/distance.conversion^2); hold all
        h_displacement{k} = histogram(displacement(:,k,1)/distance.conversion); hold all
        h_displacement{k}.Normalization = 'count';
        h_displacement{k}.FaceColor = cmap_np(k,:);
        disp(mean(displacement_squared(:,k))/2/ndim/tau)
    end
    ylabel(h_displacement{k}.Normalization)
    % xlabel(['displacement squared (' distance.units '^2)'])
    xlabel(['x displacement (' distance.units ')'])
    set(gca,'FontSize', plot_font_size)
    suptitle(title_text)
end

%% Plotting the velocity histogram
if find(strcmp(plot_options(po), 'Velocity Histogram'))
    figures{end+1} = figure;
    h_velocity = cell(ndim,1);
    h_fit = cell(ndim,1);
    legend_hv = {};
    text_hv = {};
    text_hv{end+1} = 'Gaussian fit:';
    text_hv{end+1} = 'y =  a1*exp(-((x-b1)/c1)^2)';
    % text_hv{end+1} = 'y =  a1*exp(-x^2/c1^2/2)';
    text_hv{end+1} = '';

    % v.units = '\mum/s';
    % v.conversion = 1e-6;

    v.units = [distance.units '/' time.units];
    v.conversion = distance.conversion / time.conversion;

    % disp(['2*sqrt(D/tau) = ' num2str(2*sqrt(D/tau)/v.conversion,'%.3g') ' ' v.units])

    for i = 1:ndim
        h_velocity{i} = histogram(velocity(:,:,i)/v.conversion,100); hold all

        h_velocity{i}.Normalization = 'count';
        % h_velocity{i}.Normalization = 'probability';
        h_velocity{i}.FaceColor = cmap_xyz(i,:);
        legend_hv{end+1} = ['v' coordenate_labels{i}];

        gaussian = fittype('a1*exp(-x^2/c1^2/2)');
        hv_x = h_velocity{i}.BinEdges(1:end-1) + h_velocity{i}.BinWidth/2;
        hv_y = h_velocity{i}.Values;
    %     fit_v = fit(hv_x', hv_y', gaussian, 'StartPoint', [7000,8e-4]);
        fit_v = fit(hv_x', hv_y', 'gauss1');
    %     fit_v = fit(hv_x', hv_y', 'gauss1', 'StartPoint', [300,0,3e4]);
        confint_v = confint(fit_v);

        c1_error = abs((confint_v(1,3)-confint_v(2,3))/2)*v.conversion;
    %     c1_error = abs((confint_v(1,2)-confint_v(2,2))/2)*v.conversion;

        h_fit{i} = plot(fit_v); hold all, 
        legend_hv{end+1} = ['v' coordenate_labels{i} ' fit'];
        h_fit{i}.Color = cmap_xyz(i,:);
        h_fit{i}.LineWidth = 2;
        text_hv{end+1} = ['v' coordenate_labels{i} ': a1 = ' ...
            num2str(fit_v.a1, '%.3g') ' \pm ' ...
            num2str(abs(confint_v(1,1)-confint_v(2,1))/2, '%.3g')];
        text_hv{end+1} = ['v' coordenate_labels{i} ': b1 = ' ...
            num2str(fit_v.b1, '%.3g') ' \pm ' ...
            num2str(abs(confint_v(1,2)-confint_v(2,2))/2, '%.3g') ' ' v.units];
        text_hv{end+1} = ['v' coordenate_labels{i} ': c1 = ' ...
            num2str(fit_v.c1, '%.3g') ' \pm ' ...
            num2str(abs(confint_v(1,3)-confint_v(2,3))/2, '%.3g') ' ' v.units];
    %         num2str(abs(confint_v(1,2)-confint_v(2,2))/2, '%.3g') ' ' v.units];
    %     text_hv{end+1} = ['v' coordenate_labels{i} ': D = ' ...
    %         num2str(tau/4*(fit_v.c1*sqrt(2)*v.conversion)^2*1e12, '%.3f') ' \pm ' ...
    %         num2str(tau*(fit_v.c1*sqrt(2)*v.conversion)/2*c1_error*1e12, '%.3f') ' \mum^2/s'];
        text_hv{end+1} = ['v' coordenate_labels{i} ': D = ' ...
            num2str(tau/4*(fit_v.c1*v.conversion)^2*1e12, '%.3f') ' \pm ' ...
            num2str(tau*(fit_v.c1*v.conversion)/2*c1_error*1e12, '%.3f') ' \mum^2/s'];
        text_hv{end+1} = '';

    end
    ylabel(h_velocity{1}.Normalization)
    xlabel(['velocity (' v.units ')'])
    legend(legend_hv)
    suptitle(title_text)
    grid on
    set(gca,'FontSize', plot_font_size)
%     xlim([-600,700])

    text('Units','normalized','Position',[0.08,0.95], ...
        'FontSize', 12, 'VerticalAlignment', 'top', 'String' , text_hv)
end

%% Plotting the flow velocity vs. position
if find(strcmp(plot_options(po), 'Flow Velocity vs. Position')) % && ndim == 2
    figures{end+1} = figure;
    suptitle(title_text)
    
    sp = 1;
    for i = 1:1:ndim
%         suptitle(title_text)
        subplot(ndim,2,sp)
        sp = sp + 1;
%         i = 1;
%         subplot(2,1,1)
                
        
%         if menu_vflow == 3 % contour
            scatter(reshape(position(:,:,1),[numel(position(:,:,1)),1])/distance.conversion,...
                reshape(position(:,:,2),[numel(position(:,:,2)),1])/distance.conversion, ...
                5, reshape(velocity_flow(:,:,i),[numel(velocity_flow(:,:,i)),1])/v.conversion, ...
                'filled');    
            ylabel(['y position (' distance.units ')'])
            cb = colorbar;  
            cb.Label.String = [axis_labels{i} ' flow v (' v.units ')'];
            cb.FontSize = plot_font_size;
%         else        
%             scatter(reshape(position(:,:,1),[numel(position(:,:,1)),1])/distance.conversion,...
%                 reshape(velocity_flow(:,:,i),[numel(velocity_flow(:,:,i)),1])/v.conversion, ...
%                 5,'filled');    
%             ylabel(['x flow v (' v.units ')'])
%         end
        
        xlabel(['x position (' distance.units ')'])
        set(gca,'FontSize', plot_font_size)
        grid on    
%         xlim([-20,110])
%         ylim([-20,20])
        
%         subplot(2,1,2)
        subplot(ndim,2,sp)
        sp = sp + 1;
        scatter(reshape(position(:,:,1),[numel(position(:,:,1)),1])/distance.conversion,...
            reshape(position(:,:,2),[numel(position(:,:,2)),1])/distance.conversion, ...
            5, reshape(velocity(:,:,i),[numel(velocity_flow(:,:,i)),1])/v.conversion, ...
            'filled');    
%         suptitle(title_text)
        xlabel(['x position (' distance.units ')'])
        ylabel(['y position (' distance.units ')'])
%         colormap(flipud(jet))
        colormap(jet)
        cb = colorbar;  
        cb.Label.String = [axis_labels{i} ' v (' v.units ')'];
        cb.FontSize = plot_font_size;
        set(gca,'FontSize', plot_font_size)
        grid on   
%         xlim([-20,110])
%         ylim([-20,20])
    end
%     suptitle(title_text)
end

%% CALCULATE velocity maps
% *************************************************************************

x = reshape(position(:,:,1),[numel(position(:,:,1)),1])/distance.conversion;
y = reshape(position(:,:,2),[numel(position(:,:,2)),1])/distance.conversion;

vx = reshape(velocity(:,:,1),[numel(velocity(:,:,1)),1])/v.conversion;
vy = reshape(velocity(:,:,2),[numel(velocity(:,:,2)),1])/v.conversion;

clear map_labels

menu_vmap = 5;
% menu_vmap = menu('Velocity:', 'Average', 'Maximum', 'Minimum', 'Mode', 'Gaussian Fit');

x_sections = 19;
y_sections = 11;

% x_sections = 1;
% y_sections = 3;

x_edges = linspace(min(x), max(x), x_sections + 1);
x_centres = (x_edges(1:end-1)+x_edges(2:end))/2;
x_bin = discretize(x,x_edges);
y_edges = linspace(min(y), max(y), y_sections + 1);
y_centres = (y_edges(1:end-1)+y_edges(2:end))/2;
y_bin = discretize(y,y_edges);

vx_map = zeros(numel(x_centres),numel(y_centres),1);
vy_map = zeros(numel(x_centres),numel(y_centres),1);

for i = 1:1:numel(x_centres)
    x_indices = find(x_bin == i);
    for j = 1:1:numel(y_centres)
        y_indices = find(y_bin == j);
        xy_indices = intersect(x_indices, y_indices);
%         disp(numel(xy_indices))

        if numel(xy_indices) < 30
            vx_map(i,j,:) = NaN;
            vy_map(i,j,:) = NaN;
        else
            if menu_vmap == 1 % average
                vx_map(i,j,1) = mean(vx(xy_indices));
                vy_map(i,j,1) = mean(vy(xy_indices));
                map_labels{1} = ['average v (' distance.units '/' time.units ')'];

            elseif menu_vmap == 2 % maximum
                if isempty(vx(xy_indices))
                    vx_map(i,j,1) = NaN;
                else           
                    vx_map(i,j,1) = max(vx(xy_indices));
                end
                if isempty(vy(xy_indices))
                    vy_map(i,j,1) = NaN;
                else           
                    vy_map(i,j,1) = max(vy(xy_indices));
                end
                map_labels{1} = ['max v (' distance.units '/' time.units ')'];

            elseif menu_vmap == 3 % mimimum
                if isempty(vx(xy_indices))
                    vx_map(i,j,1) = NaN;
                else           
                    vx_map(i,j,1) = min(vx(xy_indices));
                end
                if isempty(vy(xy_indices))
                    vy_map(i,j,1) = NaN;
                else           
                    vy_map(i,j,1) = min(vy(xy_indices));
                end
                map_labels{1} = ['min v (' distance.units '/' time.units ')'];

            elseif menu_vmap == 4 % mode
                vx_map(i,j,1) = mode(round(vx(xy_indices),0));
                vy_map(i,j,1) = mode(round(vy(xy_indices),0));
                map_labels{1} = ['mode v (' distance.units '/' time.units ')'];
            
            elseif menu_vmap == 5 % gaussian fit
                clc
                disp(['Calculating map section: (' num2str(i) ',' num2str(j) ') of (' ...
                    num2str(numel(x_centres)) ',' num2str(numel(y_centres)) ')'])
                
                number_of_bins = min(round(numel(xy_indices)/10), 100);
                
                [vx_counts,vx_edges] = histcounts(vx(xy_indices),number_of_bins);
                vx_centres = (vx_edges(1:end-1) + vx_edges(2:end)) / 2;
                
                fit_parameters = fit(vx_centres', vx_counts', 'gauss1');
                vx_map(i,j,1) = fit_parameters.b1; % gaussian centre
                vx_map(i,j,2) = fit_parameters.c1; % gaussian width
                vx_map(i,j,3) = fit_parameters.a1; % gaussian height
                vx_map(i,j,4) = 1/frame_rate/4*fit_parameters.c1^2; % D (units of velocity^2/time)
                
                [vy_counts,vy_edges] = histcounts(vy(xy_indices),number_of_bins);
                vy_centres = (vy_edges(1:end-1) + vy_edges(2:end)) / 2;
                
                fit_parameters = fit(vy_centres', vy_counts', 'gauss1');
                vy_map(i,j,1) = fit_parameters.b1; % gaussian centre
                vy_map(i,j,2) = fit_parameters.c1; % gaussian width
                vy_map(i,j,3) = fit_parameters.a1; % gaussian height               
                vy_map(i,j,4) = 1/frame_rate/4*fit_parameters.c1^2; % D (units of velocity^2/time)

                vx_map(i,j,5) = numel(xy_indices); % number of data points
                vy_map(i,j,5) = numel(xy_indices); % number of data points
                
                vx_map(i,j,6) = number_of_bins; % number of bins
                vy_map(i,j,6) = number_of_bins; % number of bins
                
%                 vx_map(i,j,5) = numel(xy_indices); % number of data points
%                 vy_map(i,j,5) = number_of_bins; % number of bins
                
                map_labels{1} = ['Gauss. centre v (' distance.units '/' time.units ')'];
                map_labels{2} = ['Gauss. width v (' distance.units '/' time.units ')'];
                map_labels{3} = 'Gauss. height (counts)';
                map_labels{4} = ['D (' distance.units '^2/' time.units ')'];
                map_labels{5} = 'No. of points';
                map_labels{6} = 'No. of bins';
                
%                 plot(vx_centres,vx_counts), hold all
%                 plot(vy_centres,vy_counts), hold all
            end
        end
    end
end
disp('Finished the map calculations.')

%% PLOT velocity color maps
% *************************************************************************
% close all

axis_labels = {'x','y'};
menu_vmap_axis = [1,2];
menu_vmap_plot = 1;

% [menu_vmap_plot, ~] = listdlg('PromptString', 'Velocity colour maps to plot:',...
%                            'SelectionMode', 'multiple', ...
%                            'ListString', map_labels,...
%                            'InitialValue', menu_vmap_plot);

% [menu_vmap_axis, ~] = listdlg('PromptString', 'Axis to plot:',...
%                            'SelectionMode', 'multiple', ...
%                            'ListString', axis_labels,...
%                            'InitialValue', menu_vmap_axis);
                       
[menu_vmap_plot, menu_vmap_axis] = dialog_two_lists('Select plot options:', ...
    'Variable:', map_labels, menu_vmap_plot,...
    'Axis:', axis_labels, menu_vmap_axis);
                       
map_variables = {vx_map,vy_map};

if min(size(vx_map(:,:,1))) > 5
    for i = menu_vmap_plot
        figures{end+1} = figure;
        k = 0;
        for j = menu_vmap_axis            
            k = k + 1;
            subplot(numel(menu_vmap_axis),1,k)

            % set NaN to blank
            imAlpha=ones(size(map_variables{j}(:,:,i)'));
            imAlpha(isnan(map_variables{j}(:,:,i)'))=0;

%             contour_levels = linspace(-36.2922, 22.1115, 100);
%             contour_levels = linspace(min(min(map_variables{j}(:,:,i))), max(max(map_variables{j}(:,:,i))), 100);
%             contour_levels = linspace(min(min(min(map_variables{1}(:,:,i))),min(min(map_variables{2}(:,:,i)))), ...
%                 max(max(max(map_variables{1}(:,:,i))),max(max(map_variables{2}(:,:,i)))), 100);
            contour_levels = linspace(min(min(min(map_variables{1}(2:end-1,2:end-1,i))), ...
                min(min(map_variables{2}(2:end-1,2:end-1,i)))), ...
                max(max(max(map_variables{1}(2:end-1,2:end-1,i))), ...
                max(max(map_variables{2}(2:end-1,2:end-1,i)))), 100);
            contour_levels = linspace(3.8, 12.5, 100);

            p_vx = imagesc(x_centres, y_centres, map_variables{j}(:,:,i)','AlphaData',imAlpha); 
%             set(gca,'YDir','normal')
    %         p_vx = pcolor(x_edges(1:end-1), y_edges(1:end-1), map_variables{j}(:,:,i)');
    %         p_vx.EdgeColor = 'none';

%             colormap(flipud(jet))
            colormap(jet)
            c = colorbar;
    %         c.Label.String = [map_labels{i} ' vx (' distance.units '/' time.units ')'];
            c.Label.String = [axis_labels{j} ' ' map_labels{i}];
            c.FontSize = 16;

            caxis([min(contour_levels), max(contour_levels)])

            xlabel(['x position (' distance.units ')'])
            ylabel(['y position (' distance.units ')'])
            title([title_text ' // average = ' ...
                num2str(nanmean(nanmean(map_variables{j}(:,:,i))), '%.3g')])
            set(gca, 'FontSize', 16)
        end
        
    end
end

%% SAVING FIGURES
% *************************************************************************
menu_save_figures = 1;
% menu_save_figures = menu('Save Figures?', 'NO', 'YES');
folder_path_save = 'R:\aa938\NanoPhotonics\Matlab\Brownian Motion\';
if menu_save_figures == 2    
    for i = 1:1:max(size(figures))
        if findobj(figures{i}) ~= 0
            figure_save = figures{i};
            file_name_save = '';
            
            figure(figure_save)
            pause(0.1)
            [file_name_save,folder_path_save,~] = uiputfile(['.' 'png'],...
                'File to Save the Figure',[folder_path_save file_name_save]);
            hgexport(figure_save, [folder_path_save file_name_save], hgexport('factorystyle'), 'Format', 'png')
            file_name_save = strrep(file_name_save, 'png', 'fig');    
            saveas(figure_save, [folder_path_save file_name_save], 'fig');

        end
    end
end

% % Calculate positions of particles
% x_increment = h*randn(N,np);
% x = cumsum(x_increment);
% x_start = x(1:end-1,:);
% x_end = x(2:end,:);
% vx = (x_end - x_start) ./ tau;
% 
% y_increment = h*randn(N,np);
% y = cumsum(y_increment);
% y_start = y(1:end-1,:);
% y_end = y(2:end,:);
% vy = (y_end - y_start) ./ tau;
% 
% z_increment = h*randn(N,np);
% z = cumsum(z_increment);
% z_start = z(1:end-1,:);
% z_end = z(2:end,:);
% vz = (z_end - z_start) ./ tau;

% position(:,:,1) = x;
% position(:,:,2) = y;
% position(:,:,3) = z;
% velocity(:,:,1) = vx;
% velocity(:,:,2) = vy;
% velocity(:,:,3) = vz;

% x = position(:,:,1);
% y = position(:,:,2);
% z = position(:,:,3);
% vx = velocity(:,:,1);
% vy = velocity(:,:,2);
% vz = velocity(:,:,3);

% x = zeros(1,np);
% y = zeros(1,np);
% z = zeros(1,np);
% for j=1:np
%     for i=1:N
% 
%         x(i+1,j)=x(i,j)+h*randn();
%         y(i+1,j)=y(i,j)+h*randn();
%         z(i+1,j)=z(i,j)+h*randn();
%     
%     end
% end
