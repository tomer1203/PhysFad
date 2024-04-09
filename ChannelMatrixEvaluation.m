
% clear all
close all
clc
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
ris_num_samples = 200;
if length(x_ris)~=length(y_ris)
    disp('Error: x_ris and y_ris do not have the same length.');
else
    N_RIS = length(x_ris);
end



%--- Yotam and Ron's Values ---%
% chi_ris = 6*ones(size(x_ris));
% gamma_ris = unifrnd(0,100,size(x_ris));
% fres_ris   = unifrnd(1,5,ris_num_samples,N_RIS);

%--- Sin distribution for chi and gamma ---%
% chi_phase=2*pi*rand(1,ris_num_samples);
% chi_phase=2*pi*zeros(1,ris_num_samples);
chi_phase=2*pi*linspace(0,1,ris_num_samples);
chi_frequency=0.3*rand(1,ris_num_samples);
% chi_frequency=0.3*ones(1,ris_num_samples);
% chi_ris = 50+50*sin(chi_frequency'*[1:size(x_ris,2)]+chi_phase'*ones(size(x_ris)));
% chi_ris = 25*rand(ris_num_samples,1)+25*rand(ris_num_samples,1).*sin(chi_frequency'*[1:size(x_ris,2)]+chi_phase'*ones(size(x_ris)));
% gamma_ris = 0.01+0.01*sin(chi_frequency'*[1:size(x_ris,2)]+chi_phase'*ones(size(x_ris)));
% fres_ris   = 1+0.4*sin(chi_frequency'*[1:size(x_ris,2)]+chi_phase'*ones(size(x_ris)));

%--- Normal distribution ---%
% gamma_ris = zeros(cat(2,ris_num_samples,size(x_ris)));
% chi_ris = zeros(cat(2,ris_num_samples,size(x_ris)));
% gamma_mu_linespace=linspace(-1.5,1.5,ris_num_samples);
% gamma_std_linespace=linspace(0.3,3,ris_num_samples);
% for J=1:ris_num_samples
%     gamma_ris(J,:) = normpdf(linspace(-1.5,1.5,N_RIS),0,gamma_std_linespace(J))';
%     chi_ris(J,:) = normpdf(linspace(-1.5,1.5,N_RIS),0,gamma_std_linespace(J))';
% end
% gamma_ris(isnan(gamma_ris))=0;
% chi_ris(isnan(chi_ris))=0;

%--- uniform random ---%
chi_ris = unifrnd(0.001,4,cat(2,ris_num_samples,size(x_ris))); %0-100
gamma_ris = unifrnd(0,4,cat(2,ris_num_samples,size(x_ris))); %0-100
% gamma_ris = unifrnd(0,100,cat(2,ris_num_samples,size(x_ris))); %0-100
% fres_ris   = unifrnd(0.001,5,ris_num_samples,N_RIS);%1-5

logical_fres_ris = reshape(randsample([true,false],ris_num_samples*N_RIS,true),[ris_num_samples,N_RIS]);
logical_chi_ris = reshape(randsample([true,false],ris_num_samples*N_RIS,true),[ris_num_samples,N_RIS]);
logical_gamma_ris = reshape(randsample([true,false],ris_num_samples*N_RIS,true,[0.2,0.8]),[ris_num_samples,N_RIS]);


resonant_freq = unifrnd(0.92,1.08,ris_num_samples,N_RIS);
non_resonant_freq = unifrnd(0.05,5,ris_num_samples,N_RIS);

resonant_chi = unifrnd(0.001,1,ris_num_samples,N_RIS);
non_resonant_chi = unifrnd(0.001,5,ris_num_samples,N_RIS);
resonant_gamma = unifrnd(0,0.1,ris_num_samples,N_RIS);
non_resonant_gamma = unifrnd(0,5,ris_num_samples,N_RIS);


fres_ris = logical_fres_ris.*resonant_freq+not(logical_fres_ris).*non_resonant_freq;
chi_ris = logical_fres_ris.*resonant_chi+not(logical_fres_ris).*non_resonant_chi;
gamma_ris = logical_fres_ris.*resonant_gamma+not(logical_fres_ris).*non_resonant_gamma;
% fres_ris = reshape(randsample([1,5],ris_num_samples*N_RIS,true),[ris_num_samples,N_RIS]);
%--- constants ---%
% chi_ris = 6*ones(cat(2,ris_num_samples,size(x_ris)));
% gamma_ris = 0.2*ones(cat(2,ris_num_samples,size(x_ris)));%50
% chi_ris = 0.2*ones(cat(2,ris_num_samples,size(x_ris)));
% gamma_ris = 0*ones(cat(2,ris_num_samples,size(x_ris)));%50
% fres_ris = ones(cat(2,ris_num_samples,size(x_ris)));



[Validity,N_T,N_R,N_E] = TestParameters(x_tx,y_tx,fres_tx,chi_tx,gamma_tx, ...
                                                 x_rx,y_rx,fres_rx,chi_rx,gamma_rx, ...
                                                 x_env,y_env,fres_env,chi_env,gamma_env, ...
                                                 chi_ris,gamma_ris,N_RIS);
if Validity == false
    error("One of the parameters provided is invalid");
end

%% EVALUATE CHANNEL MATRIX
avg_sampled_Hs = zeros(ris_num_samples,length(freq));
sampled_Hs = zeros(ris_num_samples,length(freq),length(x_rx),length(x_tx));
RISConfiguration = zeros(ris_num_samples,length([fres_ris(1,:),chi_ris(1,:),gamma_ris(1,:)]));
%% Debug
% load('new_full_RISConfiguration.mat','RISConfiguration');
% fres_ris = RISConfiguration(:,1:length(x_ris));
% chi_ris = RISConfiguration(:,length(x_ris)+1:2*length(x_ris));
% gamma_ris = RISConfiguration(:,2*length(x_ris)+1:3*length(x_ris));
%% Debug End
SumRateVec = zeros(ris_num_samples,1);
ppm = ParforProgressbar(ris_num_samples);
W = Calc_W(freq,x_tx,y_tx,x_rx,y_rx, x_env,y_env,x_ris,y_ris);
parfor sample=1:ris_num_samples
    
    disp(sample);
    % [~,H] = getH(freq,...
    %                 x_tx,y_tx,fres_tx,chi_tx,gamma_tx,...
    %                 x_rx,y_rx,fres_rx,chi_rx,gamma_rx,...
    %                 x_env,y_env,fres_env,chi_env,gamma_env,...
    %                 x_ris,y_ris,fres_ris(sample,:),chi_ris(sample,:),gamma_ris(sample,:));
    % [~,H2] = getH2(freq,...
    %                 x_tx,y_tx,fres_tx,chi_tx,gamma_tx,...
    %                 x_rx,y_rx,fres_rx,chi_rx,gamma_rx,...
    %                 x_env,y_env,fres_env,chi_env,gamma_env,...
    %                 x_ris,y_ris,fres_ris(sample,:),chi_ris(sample,:),gamma_ris(sample,:));
    [~,H] = getH3(freq,W,...
        fres_tx,chi_tx,gamma_tx,...
        fres_rx,chi_rx,gamma_rx,...
        fres_env,chi_env,gamma_env,...
        fres_ris(sample,:),chi_ris(sample,:),gamma_ris(sample,:));
    % disp(H3==H);
    % avg_sampled_Hs(sample,:) = mean(mean(abs(H),2),3);
    avg_sampled_Hs(sample,:) = abs(H(:,1,2));
    sampled_Hs(sample,:,:,:) = abs(H);
    SumRateVec(sample) = getSumRate(H,ones(length(H)),1);
    RISConfiguration(sample,:) = [fres_ris(sample,:),chi_ris(sample,:),gamma_ris(sample,:)];
    ppm.increment();
%     sampled_Hs(sample,:) = abs(H(:,1,1));
end
delete(ppm);
plot(avg_sampled_Hs')
title("Channel H between first transmitor to second reciever")
xlabel("freq");
% figure;
% plot(std(avg_sampled_Hs,0,1))
% xlabel("freq");
% ylabel("std(|H_{12}|)");
% figure;
% plot(SumRateVec);
% xlabel("test x");
% ylabel("channel sum rate");
% similarity_matrix = (1/101) * avg_sampled_Hs * avg_sampled_Hs'/mean(avg_sampled_Hs(:).^2);
% disp(min(similarity_matrix(:)));
[dataDiversity,~,dataSampleDiv]= getDataDiversity(sampled_Hs);
figure;
plot(dataSampleDiv);
disp('diversity '+string(dataDiversity));
figure();
plot(SumRateVec)
disp('average rate'+string(mean(SumRateVec)))
disp('max rate'+string(max(SumRateVec)))
% mean_H = mean(avg_sampled_Hs);
% mean_similarity = (1/101) * mean_H * avg_sampled_Hs'/mean(mean_H(:).^2);
% figure;
% plot(mean_similarity*RISConfiguration);
% plot(RISConfiguration(:,70)',mean_similarity,'*')
% plot(avg_sampled_Hs(abs(mean_similarity-1)>0.2,:)','LineWidth',1,'LineStyle',':')
% hold on
% plot(mean_H,'LineWidth',2.5)
% csvwrite("H_realizations.txt",avg_sampled_Hs(abs(mean_similarity-1)>0.3,:));
% csvwrite("H_residuals.txt",avg_sampled_Hs-mean_H);
% csvwrite("RISConfiguration.txt",RISConfiguration(abs(mean_similarity-1)>0.3,:));
% save("new_full_H_realizations.mat",'sampled_Hs');
% csvwrite("H_residuals.txt",avg_sampled_Hs-mean_H);
% save("new_full_RISConfiguration.mat",'RISConfiguration');