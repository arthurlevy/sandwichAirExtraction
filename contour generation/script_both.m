%This script plots the process map described in Levy, Kratz and Hubert 13.
%
%This file is part of Air Evacuation Simulation.
%
%    Inter-Strand Void Content Simulation is free software: you can
%    redistribute it and/or modify it under the terms of the GNU General
%    Public License as published by the Free Software Foundation, either
%    version 3 of the License, or (at your option) any later version.
%
%    Air Evacuation Simulation is distributed in the hope that
%    it will be useful, but WITHOUT ANY WARRANTY; without even the implied
%    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
%    the GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with Inter-Strand Void Content Simulation.  If not, see
%    <http://www.gnu.org/licenses/>.
%
% Copyright 2013 Arthur Levy

clear


%% parameters

Patm = 1.01e5;
Kx = 12e-14; % Permeability 
mu = 1.85e-5; %air viscosity
phi = 0.234; %porosity
alpha = 2; %correction factor
Kz = 5e-17; %through thickness permeability

hs = 1.2e-3;%skin thickness
hc = 20e-3; %core thickness
L = 0.02:0.001:0.3;

times = logspace(1-log10(3600),3+log10(3)-log10(3600),200);



taux = L.^2 * mu * (phi + hc/(2*hs) )  / ((1+alpha)*Kx*Patm);

tauz = ones(1,length(L)) * 2 * mu * hs * hc / Kz / Patm;

load dimensionless_solution.txt

P_trans = interp1(dimensionless_solution(:,1),dimensionless_solution(:,3), 3600*times'*(1./tauz),...
	'spline', 'extrap');

P_plane = interp1(dimensionless_solution(:,1),dimensionless_solution(:,2), 3600*times'*(1./taux),...
	'spline', 'extrap');

Mat(1,:,:) = P_trans(:,:);
Mat(2,:,:) = P_plane(:,:);

Mat(isnan(Mat)) = 1;
Mat2 = min(Mat);
Mat3(:,:) = Mat2(1,:,:);


[C,h] = contourf (L,times,P_trans>P_plane,[0.4 0.6]);
hold all
colormap([0.9 0.9 0.9]);
set(h,'color', 'none');

[C,h] = contour (L,times,Mat3,0.1:0.1:0.9);
xlabel('L [m]')
ylabel ('time [h]')
set(gca,'FontSize', 14)
set(h,'ShowText','on','LabelSpacing', 400, 'color', [0 0 0]);
set(gca,'FontSize', 14, 'Yscale','log')
%colormap([0 0 0]);
text_handle = clabel(C,h);
set(text_handle,'FontSize', 14)




% Lc = sqrt(4*hs*hs*(1+alpha) *Kx/Kz);

%vline(Lc)

