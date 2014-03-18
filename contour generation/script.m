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




%% parameters

Patm = 1.01e5;
Kx = 6e-14; % Permeability 
mu = 1.85e-5; %air viscosity
phi = 0.234; %porosity
alpha = 0.32; %correction factor


char_time_span = logspace(3, 6+log10(2), 400);

%% %%%%%%%%%%%%%%%%%%%%%%%% The end of line pressure contour
second = subplot('position', [0.075 0.06 0.86 0.32 ]);

t = 0:60*10:3600*72;

tau = t'*(1./char_time_span);

load dimensionless_solution.txt

Mat = interp1(dimensionless_solution(:,1),dimensionless_solution(:,2), tau);

[C,h] = contour(char_time_span,t/3600,Mat,0.1:0.1:0.9,'LineWidth',2);


% formatting
xlabel('Characteristic time \tau [s]','FontSize', 14);

set(second, 'Xscale', 'log', 'FontSize', 14);

ylabel('time (h)');
colormap([0,0,0]);
set(h,'ShowText','on','LabelSpacing', 400);
set(gca,'FontSize', 14)
text_handle = clabel(C,h);
set(text_handle,'FontSize', 14)
h = title('Highest pressure in part [Atm].');
set(h,'Position',get(h,'Position') - [0 5 0])



%% %%%%%%%%%%%%%%%The iso L^2/K contour , general case
general = subplot('position', [0.075 0.72  0.86 0.23 ]);

L_K_span = [1e12 2e12 5e12 1e13 2e13 5e13 1e14 2e14 5e14 1e15 2e15];
ratio_span = 0:0.25:30;

denom = (char_time_span).^(-1)'*mu*(phi+ratio_span/2);

Mat_L2_K = (1.5*Patm./denom);


[C,h] = contour(char_time_span, ratio_span,...
    Mat_L2_K'/1e14,   L_K_span/1e14);

set(h,'ShowText','on','LabelSpacing', 500');
text_handle = clabel(C,h);
set(text_handle,'FontSize', 14)
% formatting
set(general,'XTick',[],'xscale','log','XAxisLocation', 'top')
set(general,'XAxisLocation', 'top',...
    'Xscale', 'log',...
    'FontSize', 14,...
	'XTick', interp1(Mat_L2_K(:,end),char_time_span,L_K_span),...
    'XTickLabel', [],...num2str(L_K_span',1),...
	'XMinorTick', 'off',...
    'box', 'off','TickLength',[0 0]);

ylabel('h_c/h_s');
xlabel('General case. Dimensionless number L^2/K (\times10^{14})');
%xlabh = get(general,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 7 0])
colormap([0,0,0]);

%additional axes for ticks
h = axes('position', get(general,'Position'));
set(h, 'XTick', get(second,'XTick'),...
    'XLim',get(second,'XLim'), 'XTickLabel' , '',...
	'YTickLabel' , [],...
	'YLim', [1e15 2.9e15],...
	'YTick', 2e15,'FontSize', 14,...
	'YAxisLocation', 'right',...
    'Color', 'none',...
    'Xscale', 'log',...
    'TickLength', [0 0])

h2 =  axes('position', get(general,'Position'));
set(h2, 'YTick', get(general,'YTick'),...
	'YLim',get(general,'YLim'),...
	'XTickLabel' , '',...
	'XTick', [],...
	'Color', 'none',...
	'FontSize', 14, 'box', 'on')

h3 =  axes('position', get(general,'Position'));
set(h3, 'XTick', get(second,'XTick'),...
	'XLim',get(second,'XLim'),...
	'XTickLabel' , '',...
	'XScale', 'log',...
	'YTick', [],...
	'Color', 'none',...
	'FontSize', 14, 'box', 'off')

%% %%%%%%%% The iso L contour, 5320 case
cycom = subplot('position', [0.075 0.43  0.86 0.23 ]);

L = [0.2 0.25 0.3 0.35 0.4 0.5 0.6 0.7 0.8 1 1.2 1.5 2 2.5 3 4 5 7 10];

denom = (char_time_span).^(-1)'*mu*(phi+ratio_span/2);
Mat_L = (Kx*(1+alpha)*Patm./denom).^(1/2);

[C,h] = contour(char_time_span,ratio_span,Mat_L',L);
set(h,'ShowText','on','LabelSpacing', 500');
text_handle = clabel(C,h);
set(text_handle,'FontSize', 14)

%formatting
set(cycom,'XTick',[]);
set(cycom,'XAxisLocation', 'top',...
    'Xscale', 'log',...
    'FontSize', 14,...
	'XTick', interp1(Mat_L(:,end),char_time_span,L),...
    'XTickLabel', [],...
	'XMinorTick', 'off',...
    'box', 'off','TickLength',[0 0]);
ylabel('h_c/h_s');
colormap([0,0,0]);
xlabel('Cycom 5320 PW case. Distance from edge breathing L [m]')
%xlabh = get(cycom,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 7 0])


% additional grid
h = axes('position', get(cycom,'Position'));
set(h, 'XTick', get(second,'XTick'),...
    'XLim',get(second,'XLim'), 'XTickLabel' , '',...
	'YTickLabel' , [],...
	'YLim', [ 0,1],...
	'YTick', 0.58,'FontSize', 14,...
	'YAxisLocation', 'right',...
    'Color', 'none',...
    'Xscale', 'log',...
    'TickLength', [0 0])

h2 =  axes('position', get(cycom,'Position'));
set(h2, 'YTick', get(cycom,'YTick'),...
	'YLim',get(cycom,'YLim'),...
	'XTickLabel' , '',...
	'XTick', [],...
	'Color', 'none',...
	'FontSize', 14, 'box', 'on')

h3 =  axes('position', get(cycom,'Position'));
set(h3, 'XTick', get(second,'XTick'),...
	'XLim',get(second,'XLim'),...
	'XTickLabel' , '',...
	'XScale', 'log',...
	'YTick', [],...
	'Color', 'none',...
	'FontSize', 14, 'box', 'off')

%% Background grid
back = axes('position', [0.075 0.06  0.86 0.89 ]);

set(back, 'XTick', get(second,'XTick'),...
    'XLim',get(second,'XLim'),...
	'XTickLabel' , '',...
	'YTickLabel' , '',...
	'Xgrid', 'on',...
    'Color', 'none',...
	'Box', 'on',...
    'Xscale', 'log',...
    'TickLength', [0 0])


set(gcf,'Pointer','fullcross')


