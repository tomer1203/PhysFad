
% clear all
close all

% figure
disp('Example of a Channel Matrix Evaluation');
%% DEFINITION OF PHYSFAD PARAMETERS

freq = linspace(0.9,1.1,120);

%% Configurable Dipole Properties
%## Transmitters ##
% locations
x_tx = [0 0 0];
y_tx = [4 4.5 5];
% dipole properties
fres_tx = [1 1 1];
chi_tx = [0.5 0.5 0.5];
gamma_tx = [0 0 0];

%##  Receivers ##
% locations
x_rx = [15 15 15 15];
y_rx = [11 11.5 12 12.5];
% properties
fres_rx = [1 1 1 1];
chi_rx = [0.5 0.5 0.5 0.5];
gamma_rx = [0 0 0 0];

%## Scattering Environment ##
% locations
load('ComplexEnclosure2.mat')
% properties
fres_env = 10*ones(size(x_env));%10*ones maybe change this to 2
chi_env = 50*ones(size(x_env));% 50*ones
gamma_env = 0*ones(size(x_env));% 100*ones

%## RIS ##
% locations
load('ExampleRIS3.mat','x_ris','y_ris');
% x_ris = x_ris(1,1:1);%1,1:25
% y_ris = y_ris(1,1:1);%1,1:25
% properties
ris_num_samples = 15000;
if length(x_ris)~=length(y_ris)
    disp('Error: x_ris and y_ris do not have the same length.');
else
    N_RIS = length(x_ris);
end



%--- Read Values from python simulation---%
RISConfiguration = csvread("Physfad_optimal_parameters.txt");
% RISConfiguration = csvread("optimal_parameters.txt");
% H_python = csvread("Physfad_H_mat.txt");
% H_python = H_python.reshape(120,4,3);
% full_ris_configuration = csvread("full_RISConfiguration.txt");
% full_H_realizations = csvread("full_H_realizations.txt");
% Y1 = reshape(full_H_realizations(1,:),120,4,3);
% RISConfiguration = full_ris_configuration(1,:)';
% load('H_python_mat.mat','H_python_mat');
% load('estOptInp.mat','estOptInp');
RISConfiguration = RISConfiguration.^2;
fres_ris = RISConfiguration(1:length(x_ris))';
chi_ris = RISConfiguration(length(x_ris)+1:2*length(x_ris))'+0.00000001;
gamma_ris = RISConfiguration(2*length(x_ris)+1:3*length(x_ris))';


%% Testing Parameter Validity
if length(x_tx)~=length(y_tx)
    disp('Error: x_tx and y_tx do not have the same length.');
else
    N_T = length(x_tx);
end
if length(fres_tx)~=N_T
    disp('Error: x_tx and fres_tx do not have the same length.');
end
if length(chi_tx)~=N_T
    disp('Error: x_tx and chi_tx do not have the same length.');
end
if length(gamma_tx)~=N_T
    disp('Error: x_tx and gamma_tx do not have the same length.');
end
if length(x_rx)~=length(y_rx)
    disp('Error: x_rx and y_rx do not have the same length.');
else
    N_R = length(x_rx);
end
if length(fres_rx)~=N_R
    disp('Error: x_rx and fres_rx do not have the same length.');
end
if length(chi_rx)~=N_R
    disp('Error: x_rx and chi_rx do not have the same length.');
end
if length(gamma_rx)~=N_R
    disp('Error: x_rx and gamma_rx do not have the same length.');
end
if length(x_env)~=length(y_env)
    disp('Error: x_env and y_env do not have the same length.');
else
    N_E = length(x_env);
end
if length(fres_env)~=N_E
    disp('Error: x_env and fres_env do not have the same length.');
end
if length(chi_env)~=N_E
    disp('Error: x_env and chi_env do not have the same length.');
end
if length(gamma_env)~=N_E
    disp('Error: x_env and gamma_env do not have the same length.');
end

if length(chi_ris(1,:))~=N_RIS
    disp('Error: x_ris and chi_ris do not have the same length.');
end
if length(gamma_ris(1,:))~=N_RIS
    disp('Error: x_ris and gamma_ris do not have the same length.');
end


%% EVALUATE CHANNEL MATRIX
sampled_Hs = zeros(ris_num_samples,length(freq),length(x_rx),length(x_tx));
% [fres_ris,chi_ris,gamma_ris] = RISConfiguration;%TODO: FIX THIS%
[~,H] = getH2(freq,...
                x_tx,y_tx,fres_tx,chi_tx,gamma_tx,...
                x_rx,y_rx,fres_rx,chi_rx,gamma_rx,...
                x_env,y_env,fres_env,chi_env,gamma_env,...
                x_ris,y_ris,fres_ris,chi_ris,gamma_ris);
% avg_sampled_Hs(sample,:) = mean(mean(abs(H),2),3);
avg_sampled_H = abs(H(:,1,2));
SumRate = getSumRate(H,ones(length(H)),1)
plot(avg_sampled_H)
