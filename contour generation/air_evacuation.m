%This function plots the pressure at end of line for a sandwich structure
%under vacuum hold.
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


function air_evacuation

%% parameters

Patm = 1.01e5; %(Pa)
Kx = 6e-14; % Permeability (m^2)
mu = 1.85e-5; %air viscosity (Pa.s)
phi = 0.234; %porosity
alpha = 0.32; %correction factor
L = 2;%length (m)
hc = 20e-3; %core thickness (m)
hs = 2e-3; %skin thckness (m)


load dimensionless_solution.txt

%times:
t = 0:3600*72;

tau = L^2 * mu * (phi+hc/2/hs)  / (  (1+alpha)*Kx*Patm  ); 

tstar =t/tau;

Pressure = interp1(...
	dimensionless_solution(:,1), dimensionless_solution(:,2),...
	tstar);

plot(t/3600, Pressure);
xlabel('time (h)')
ylabel('Pressure (atm)')

end

