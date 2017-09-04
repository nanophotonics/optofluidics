clc
clear 
close all

% Gaussian vector:
% (q) where 1/q = 1/R - i/n/zR 
% (1)

% Focused gaussian spot: R = infinity
% Rayleigh length:
% zR = pi * w^2 / lambda

% Once we have the total ABCD matrix:
% q2 = (A*q1 * B) / (C*q1 + D)
 
% Propagation in free space:
% (A, B) = (1 d)
% (C, D)   (0 1)

% Propagation through a thin lens:
% (A, B) = (1    0)
% (C, D)   (-1/f 1)

% ALL UNITS IN SI

% Configuration:
% PM fibre >> d1 >> Objective 1 >> d2 >> Objective 2 >> d3 >> HC fibre
%             M1        M2         M3         M4        M5

% PM fibre >> d1 >> Objective 1 >> d2 >> Lens 2 >> d3 >> Lens 3 >> d4 >> Objective 4 >> d5 >> HC fibre
%             M1        M2         M3      M4      M5     M6       M7       M8          M9        

lambda = 810e-9; 

w_TiSa = 0.2e-3;
zR_TiSa = pi * w_TiSa^2 / lambda;

w_PM = 2.45e-6;
zR_PM = pi * w_PM^2 / lambda;

q1 = 1i*zR_TiSa;
% q1 = 1i*zR_PM;

% ---
f20x = 9e-3;
d20x = 1.2e-3;

f10x = 18e-3;
d10x = 10.6e-3;

f4x = 45e-3;
d4x = 18.5e-3;

f20mm = 20e-3;

% ---

d1 = 1.75;
f1 = f20mm;
d2 = f1;
f2 = f10x;
d3 = d10x*1.5;

nx = 100;
ny = 200;

x = zeros(nx,1);
y = zeros(ny,1);
z = zeros(nx,ny);

ii = 0;
for d1 = linspace(0.5, 2.5, nx)
    ii = ii + 1;
    x(ii) = d1;
    jj = 0;
    for d2 = linspace(f1*0.8, f1*1.2, ny)
        jj = jj + 1;
        if ii == 1;
            y(jj) = d2;
        end
        

        M1 = [1,d1; 0,1];                                                                                                                                                                                                                             
        M2 = [1,0; -1/f1,1];
        M3 = [1,d2; 0,1];
        M4 = [1,0; -1/f2,1];
        M5 = [1,d3; 0,1];

%         M = M2 * M1;
        M = M3 * M2 * M1;
        % M = M5 * M4 * M3 * M2 * M1;
        % M = M1 * M2 * M3 * M4 * M5;

        % ---

        % f1 = f20x;
        % f2 = 0.25;
        % f3 = 0.10;
        % f4 = f4x;
        % 
        % d1 = d20x;
        % d2 = 0.2;
        % d3 = 0.35;
        % d4 = 1;
        % d5 = d4x;
        % 
        % M1 = [1,d1; 0,1];
        % M2 = [1,0; -1/f1,1];
        % M3 = [1,d2; 0,1];
        % M4 = [1,0; -1/f2,1];
        % M5 = [1,d3; 0,1];
        % M6 = [1,0; -1/f3,1];
        % M7 = [1,d4; 0,1];
        % M8 = [1,0; -1/f4,1];
        % M9 = [1,d5; 0,1];
        % 
        % M = M9 * M8 * M7 * M6 * M5 * M4 * M3 * M2 * M1;

        % ---

        A = M(1,1);
        B = M(1,2);
        C = M(2,1);
        D = M(2,2);

        q2 = (A*q1 + B) / (C*q1 + D);
        
        zR_PM = q2/1i;
        w_PM = real(sqrt(zR_PM *lambda / pi));
        
        % zR_HC = q2/1i;        
        % w_HC = real(sqrt(zR_HC *lambda / pi));

    z(ii,jj) = w_PM;
    % z(ii,jj) = w_HC;
    end
end

disp(['w TiSa = ' num2str(w_TiSa*1e3) ' mm'])
disp(['w PM = ' num2str(w_PM*1e6) ' um'])
% disp(['w HC = ' num2str(w_HC*1e6) ' um'])


% plot(x,y*1e6, 'LineWidth', 2)
% xlabel('d2 (m)')
% ylabel('w HC (\mum)')

contourf(x,y*1e3,z'*1e6,100,'LineColor','none')
xlabel('d1 (m)')
ylabel('d2 (mm)')
title('w PM (\mum)')
colormap(jet)
colorbar

% w_test = lambda * f1 / pi / w_PM;
% disp(['w test = ' num2str(w_test*1e6) ' um'])
