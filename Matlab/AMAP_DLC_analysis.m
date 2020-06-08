% Author: Daniele, la Cecilia (The University of Sydney @ daniele.lacecilia@sydney.edu.eu, now at eawag)
% Author: Giovanni M., Porta (Politecnico di Milano @ giovanni.porta@polimi.it)

% Should you use the AMAP function in your analysis, please cite the following publication
% Data and script were used in the journal article:
% Daniele la Cecilia, Giovanni M. Porta, Fiona H.M. Tang, Monica Riva, Federico Maggi
% Probabilistic indicators for soil and groundwater contamination risk assessment
% Ecological Indicators
% Volume 115
% 2020
% 106424
% https://doi.org/10.1016/j.ecolind.2020.106424
% http://www.sciencedirect.com/science/article/pii/S1470160X20303617

% Should you use the AMA family (AMAE, AMAV, AMAskew, AMAkurt) in your analysis, please cite the following publication
% The AMA family of indicators was introduced in:
% Aronne Dell'Oca, Monica Riva, Alberto Guadagnini
% Moment-based metrics for global sensitivity analysis of hydrological systems
% Hydrology and Earth System Sciences
% Volume 21
% 2017
% 6219-6234
% https://doi.org/10.5194/hess-21-6219-2017
% https://www.hydrol-earth-syst-sci.net/21/6219/2017/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function computes AMA sensitivity indices and principal Sobol Indices (Si) and AMAP
% starting from a set of Monte Carlo realizations. It also gives
% conditional moments with respect to one parameter at a time up to order 4
% inputs:

% The function needs a matrix called: Sampling_points , containing MC realizations of the parameters [n simulations x n parameters]
%                                     Output_Mat: Full model evaluations [n simulations x n time points]
% The user defines: nclass , the number of classes used to compute conditional statistics
%                   thr: threshold value
% The script automatically get: name_var , the name of each output variable (text string)

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT

% which indicator to be plotted? replace "AMAP" with any indicator you wish to plot
% indicators: AMAP, AMAE, AMAV, AMAskew, AMAkurt, Si
var_out = 'AMAP'; % chosen by the user

dotime = 0; % chosen by the user
% dotime = 0: calculate indicators at a fixed time
% dotime = 1: calculate indicators over time

%% THIS IS YOUR MATRIX Output_mat (for the example without analysis over time)
% For simplicity, instead of having the whole time series for all simulations, we reported 2 cases:
file = 'Matrice_REF_C_GLP_15years.csv'; % predicted concentrations in microg/l at a fixed time = 15 years
% organized as: each row is a different simulation, each column is a different time
opts = detectImportOptions(file, 'NumHeaderLines', 1, 'Delimiter',';');
opts.VariableNamesLine = 1;
% opts.VariableNames
GLP_aq_BRZ_at15y = readmatrix(file, opts);

%open file
fid = fopen(file,'r');
%Read the first line of the file containing the header information
headerline = fgetl(fid);
%close file
fclose(fid);
headers = textscan(headerline,'%s','Delimiter',';');
labout_at5 = headers{1,1};

%% OR THIS IS YOUR MATRIX Output_mat (for the example with analysis over time)
file = 'Matrice_REF_C_GLP_by5years.csv'; % predicted concentrations in microg/l over time
% In this example, predicted concentrations every 5 years
% organized as: each row is a different simulation, each column is a different time
opts = detectImportOptions(file, 'NumHeaderLines', 1, 'Delimiter',';');
opts.VariableNamesLine = 1;
GLP_aq_BRZ_by5y = readmatrix(file, opts);

%open file
fid = fopen(file,'r');
%Read the first line of the file containing the header information
headerline = fgetl(fid);
%close file
fclose(fid);
headers = textscan(headerline,'%s','Delimiter',';');
labout_by5 = headers{1,1};

%% THIS IS YOUR MATRIX sampling_points
file = 'Matrice_parametri.csv'; % values of uncertain parameters sampled from a uniform distribution using a Quasi Monte-Carlo technique
opts = detectImportOptions(file, 'NumHeaderLines', 1, 'Delimiter',';');
opts.VariableNamesLine = 1;
sampling_points = readmatrix(file, opts);
% organized as: each row is a different simulation, each column is a different parameter
% names(sampling_points) <- c("k","b","psi_s","phi","Sr_l","Sr_g")
% name_var = {'k','b','\psi_s','\phi','S_{lr}','S_{gr}'};

%open file
fid = fopen(file,'r');
%Read the first line of the file containing the header information
headerline = fgetl(fid);
%close file
fclose(fid);
headers = textscan(headerline,'%s','Delimiter',';');
name_var = headers{1,1};

nclass = 20; % to be chosen by the user, in theory, the more the better

thr = 0.1; % to be chosen by the user
% In this example, thr corresponds to 0.1 microg/l, the safety threshold concentration for a concerning plant protection product in Europe
% For a mixture of concerning plant protection products and transformation products, thr <- 0.5
% For fluxes of concerning PPP, we used thr = 0.2*1000 * 0.01/100; % mg/m2/year, that is 0.01% of a gross application rate of 2 kg/ha/year
% For fluxes of mixtures of concerning PPP and transformation product, we used thr = 0.2*1000 * 0.05/100; % mg/m2/year, that is 0.05% of a gross application rate of 2 kg/ha/year given that only 1 active ingredient and 1 transformation product were accounted for in this analysis

%%% END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% define plotting parameters
mycol = [0 0 0;    % black
         0.66 0.66 0.66; % dark grey
         1 0.84 0;   % gold
         1 0 0 ; %firebrick
         0 0 0.8; % medium blue
         0.53  0.81 0.92];  % sky blue
         
mylnwidth = [1 1.5 1 1 1 1.5]; 
labelfsize = 20;
mytitle = {'Managed (REF)'};
mysort = 1;
mylegend = mytitle(mysort);
nscenario = length(mysort);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (dotime == 1)
    Output_Mat = GLP_aq_BRZ_by5y; % Output_Mat goes in the function to calculate the indicators
else
    Output_Mat = GLP_aq_BRZ_at15y;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RUN FUNCTION
[AMAP, AMAE, AMAV, AMAskew, AMAkurt, Si, cond_mean, cond_var, cond_skew, cond_kurt ,cond_PThr, number, lower, upper] = gsa_amaP_DLC(sampling_points, Output_Mat, nclass,thr, name_var);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOT
Npar = size(sampling_points,2);

tmp = AMAP;
if (dotime == 0)
    tmp = [tmp,zeros(Npar,1)];
    tmp = tmp';
    x_pos = 1;
    labout = labout_at5;
else
    x_pos = 1:size(tmp,2);
    tmp = tmp'; % for plotting
    labout = labout_by5;
end
bar(tmp,'stacked' );
grid on
box on
ylabel(var_out,'fontname','Times New Roman','fontsize',labelfsize+2)
legend(name_var,'Location','northeast')

xticks(x_pos)
xticklabels(labout)


