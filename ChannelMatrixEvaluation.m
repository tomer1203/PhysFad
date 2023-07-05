
% clear all
close all
clc

disp('Example of a Channel Matrix Evaluation');
%% DEFINITION OF PHYSFAD PARAMETERS

freq = linspace(0.9,1.1,101);

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
load('ComplexEnclosure.mat')
% properties
fres_env = 10*ones(size(x_env));
chi_env = 50*ones(size(x_env));
gamma_env = 100*ones(size(x_env));

%## RIS ##
% locations
load('ExampleRIS.mat','x_ris','y_ris');
% properties
ris_num_samples = 5;

fres_ris_ON = 1*ones(ris_num_samples,1);
fres_ris_OFF = 5*ones(ris_num_samples,1);
chi_phase=2*pi*rand(1,ris_num_samples);
chi_frequency=0.1*rand(1,ris_num_samples);
chi_ris = 25+25*sin(chi_frequency'*[1:size(x_ris,2)]+chi_phase'*ones(size(x_ris)));
% chi_ris = 10*rand(cat(2,ris_num_samples,size(x_ris)));%50
% gamma_ris = 1*rand(cat(2,ris_num_samples,size(x_ris)));
gamma_ris = 1*sin(chi_frequency'*[1:size(x_ris,2)]+chi_phase'*ones(size(x_ris)));


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
if length(x_ris)~=length(y_ris)
    disp('Error: x_ris and y_ris do not have the same length.');
else
    N_RIS = length(x_ris);
end
if length(chi_ris)~=N_RIS
    disp('Error: x_ris and chi_ris do not have the same length.');
end
if length(gamma_ris)~=N_RIS
    disp('Error: x_ris and gamma_ris do not have the same length.');
end


%% EVALUATE CHANNEL MATRIX
sampled_Hs = zeros(ris_num_samples,length(freq));
for sample=1:ris_num_samples
    % RIS Configuration
    config_ris = round(rand(1,N_RIS));
    clear fres_ris;
    for cc=1:length(config_ris)
        if config_ris(cc)==0
            fres_ris(cc) = fres_ris_OFF(sample);
        elseif config_ris(cc)==1
            fres_ris(cc) = fres_ris_ON(sample);
        end
    end
    if length(fres_ris)~=N_RIS
        disp('Error: x_ris and fres_ris do not have the same length.');
    end
    disp(sample);
    [freq,H] = getH(freq,...
                    x_tx,y_tx,fres_tx,chi_tx,gamma_tx,...
                    x_rx,y_rx,fres_rx,chi_rx,gamma_rx,...
                    x_env,y_env,fres_env,chi_env,gamma_env,...
                    x_ris(sample),y_ris(sample),fres_ris(sample),chi_ris(sample),gamma_ris(sample));
    sampled_Hs(sample,:) = mean(mean(abs(H),2),3);
%     sampled_Hs(sample,:) = abs(H(:,1,1));
end
plot(sampled_Hs')
similarity_matrix = (1/101) * sampled_Hs * sampled_Hs'/mean(sampled_Hs(:).^2);
disp(min(similarity_matrix(:)));