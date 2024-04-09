
% clear all
close all
clc

disp('Example of a Channel Matrix Evaluation');
%% DEFINITION OF PHYSFAD PARAMETERS

freq = linspace(0.9,1.1,120);
% freq = [1];

    %% Dipole Properties

        %% Transmitters

        % locations

        x_tx = [0 0 0];
        y_tx = [4 4.5 5];
        % x_tx = [0];
        % y_tx = [4];

        if length(x_tx)~=length(y_tx)
            disp('Error: x_tx and y_tx do not have the same length.');
        else
            N_T = length(x_tx);
        end

        % dipole properties

        fres_tx = [1 1 1];
        chi_tx = [0.5 0.5 0.5];
        gamma_tx = [0 0 0];
        % fres_tx = [1];
        % chi_tx = [0.5];
        % gamma_tx = [0];

        if length(fres_tx)~=N_T
            disp('Error: x_tx and fres_tx do not have the same length.');
        end
        if length(chi_tx)~=N_T
            disp('Error: x_tx and chi_tx do not have the same length.');
        end
        if length(gamma_tx)~=N_T
            disp('Error: x_tx and gamma_tx do not have the same length.');
        end


        %% Receivers
        x_sample_density = 26*5;
        y_sample_density = 22*5;
        % locations
        x_diff = repmat(linspace(-20,5,x_sample_density)',1,y_sample_density)';
        y_diff = repmat(linspace(-14.5,6.5,y_sample_density),x_sample_density,1)';
        x_rx = [15 15 15 15];
        y_rx = [11 11.5 12 12.5];

        if length(x_rx)~=length(y_rx)
            disp('Error: x_rx and y_rx do not have the same length.');
        else
            N_R = length(x_rx);
        end

        % properties

        fres_rx = [1,1,1,1];
        chi_rx = [0.5,0.5,0.5,0.5];
        gamma_rx = [0,0,0,0];

        if length(fres_rx)~=N_R
            disp('Error: x_rx and fres_rx do not have the same length.');
        end
        if length(chi_rx)~=N_R
            disp('Error: x_rx and chi_rx do not have the same length.');
        end
        if length(gamma_rx)~=N_R
            disp('Error: x_rx and gamma_rx do not have the same length.');
        end


        %% Scattering Environment

        % locations
        load('ComplexEnclosure2.mat')

        if length(x_env)~=length(y_env)
            disp('Error: x_env and y_env do not have the same length.');
        else
            N_E = length(x_env);
        end

        % properties

        fres_env = 10*ones(size(x_env));
        chi_env = 50*ones(size(x_env));
        gamma_env = 100*ones(size(x_env));

        if length(fres_env)~=N_E
            disp('Error: x_env and fres_env do not have the same length.');
        end
        if length(chi_env)~=N_E
            disp('Error: x_env and chi_env do not have the same length.');
        end
        if length(gamma_env)~=N_E
            disp('Error: x_env and gamma_env do not have the same length.');
        end

        %% RIS

        % locations
        load('ExampleRIS3.mat','x_ris','y_ris');
%         x_ris = x_ris(1,1:25);
%         y_ris = y_ris(1,1:25);
%         save('ExampleRIS3.mat','x_ris','y_ris');


        if length(x_ris)~=length(y_ris)
            disp('Error: x_ris and y_ris do not have the same length.');
        else
            N_RIS = length(x_ris);
        end

        % properties

        % fres_ris_ON = 1;
        % fres_ris_OFF = 5;
        % chi_ris = 50*ones(size(x_ris));%50
        % gamma_ris = 1*zeros(size(x_ris));
        % 
        % if length(chi_ris)~=N_RIS
        %     disp('Error: x_ris and chi_ris do not have the same length.');
        % end
        % if length(gamma_ris)~=N_RIS
        %     disp('Error: x_ris and gamma_ris do not have the same length.');
        % end

        %% RIS Configuration
        
        % config_ris = round(rand(1,N_RIS));
        % clear fres_ris;
        % for cc=1:length(config_ris)
        %     if config_ris(cc)==0
        %         fres_ris(cc) = fres_ris_OFF;
        %     elseif config_ris(cc)==1
        %         fres_ris(cc) = fres_ris_ON;
        %     end
        % end
        RISConfiguration = csvread("Physfad_optimal_parameters.txt");
        RISConfiguration = RISConfiguration.^2;
        fres_ris = RISConfiguration(1:length(x_ris))';
        chi_ris = RISConfiguration(length(x_ris)+1:2*length(x_ris))'+0.00000001;
        gamma_ris = RISConfiguration(2*length(x_ris)+1:3*length(x_ris))';
        if length(fres_ris)~=N_RIS
            disp('Error: x_ris and fres_ris do not have the same length.');
        end

    %% Visualize Dipole Locations
    % 
    % figure, hold on,box on,
    % plot(x_tx,y_tx,'bo','displayname','TX');
    % plot(x_rx,y_rx,'ro','displayname','RX');
    % plot(x_env,y_env,'k.','displayname','Scat. Env.');
    % plot(x_ris,y_ris,'g.','markersize',7.5,'displayname','RIS');
    % axis equal;
    % xlabel('x [a.u.]');
    % ylabel('y [a.u.]');
    % set(gca,'fontsize',15);
    % xlim([min([x_tx x_rx x_env x_ris])-1 max([x_tx x_rx x_env x_ris])+1]);
    % ylim([min([y_tx y_rx y_env y_ris])-1 max([y_tx y_rx y_env y_ris])+1]);
    % legend('show','location','eastoutside');
    % drawnow;
    figure
    imagesc(spatial_field_dist)
    colormap(gca, turbo)
    hold on
    plot(5.2*(x_tx+5),5.238*(y_tx+2.5),'bo','MarkerSize',8,'displayname','TX','LineWidth',1.5);
    plot(5.2*(x_rx+5),5.238*(y_rx+2.5),'ro','MarkerSize',8,'displayname','RX','LineWidth',1.5);
    plot(5.2*(x_env+5),5.238*(y_env+2.5),'k.','MarkerSize',12,'displayname','Scat. Env.');
    plot(5.2*(x_ris+5),5.238*(y_ris+2.5),'g.','markersize',12,'displayname','RIS');
    xlabel('x [a.u.]');
    ylabel('y [a.u.]');
    lgd = legend('show','location','southwest');
    fontsize(lgd,26,'points')
    drawnow;
    figure
    imagesc(spatial_field_dist_optimized)
    colormap(gca, turbo)
    hold on
    plot(5.2*(x_tx+5),5.238*(y_tx+2.5),'bo','MarkerSize',8,'displayname','TX','LineWidth',1.5);
    plot(5.2*(x_rx+5),5.238*(y_rx+2.5),'ro','MarkerSize',8,'displayname','RX','LineWidth',1.5);
    plot(5.2*(x_env+5),5.238*(y_env+2.5),'k.','MarkerSize',12,'displayname','Scat. Env.');
    plot(5.2*(x_ris+5),5.238*(y_ris+2.5),'g.','markersize',12,'displayname','RIS');
    xlabel('x [a.u.]');
    ylabel('y [a.u.]');
    lgd = legend('show','location','southwest');
    fontsize(lgd,26,'points')
    drawnow;
    figure
    imagesc((spatial_field_dist_optimized-spatial_field_dist).^2)
    colormap(gca, turbo)
    hold on
    plot(5.2*(x_tx+5),5.238*(y_tx+2.5),'bo','MarkerSize',8,'displayname','TX','LineWidth',1.5);
    plot(5.2*(x_rx+5),5.238*(y_rx+2.5),'ro','MarkerSize',8,'displayname','RX','LineWidth',1.5);
    plot(5.2*(x_env+5),5.238*(y_env+2.5),'k.','MarkerSize',12,'displayname','Scat. Env.');
    plot(5.2*(x_ris+5),5.238*(y_ris+2.5),'g.','markersize',12,'displayname','RIS');
    xlabel('x [a.u.]');
    ylabel('y [a.u.]');
    lgd = legend('show','location','southwest');
    fontsize(lgd,26,'points')
    drawnow;

    plot(x_tx,y_tx,'bo','displayname','TX');
    plot(x_env,y_env,'k.','displayname','Scat. Env.');
    plot(x_ris,y_ris,'g.','markersize',7.5,'displayname','RIS');
    axis equal;
    xlabel('x [a.u.]');
    ylabel('y [a.u.]');
    set(gca,'fontsize',15);
    % xlim([min([x_tx x_rx x_env x_ris])-1 max([x_tx x_rx x_env x_ris])+1]);
    % ylim([min([y_tx y_rx y_env y_ris])-1 max([y_tx y_rx y_env y_ris])+1]);
    legend('show','location','eastoutside');
    drawnow;
    hold off
    
%% EVALUATE CHANNEL MATRIX
disp("starting loop")
spatial_field_dist = zeros([y_sample_density,x_sample_density]);
parfor x=1:x_sample_density
    for y=1:y_sample_density
        disp([string(y)," ",string(x)])
        [~,H] = getH2(freq,...
                x_tx,y_tx,fres_tx,chi_tx,gamma_tx,...
                x_rx+x_diff(y,x),y_rx+y_diff(y,x),fres_rx,chi_rx,gamma_rx,...
                x_env,y_env,fres_env,chi_env,gamma_env,...
                x_ris,y_ris,fres_ris,chi_ris,gamma_ris);
        spatial_field_dist(y,x) = getSumRate(H,ones(length(H)),1);
        % spatial_field_dist(x,y) = abs(H(1,1));
    end
end

save("spatial_field_optimized.mat","spatial_field_dist")
% plot(mean(mean(abs(H),2),3));
% hold on;
% plot(mean(mean(abs(H2),2),3));
% H_avg = mean(mean(abs(H),2),3);
% H2_avg = mean(mean(abs(H2),2),3);
% similarity_matrix = (1/101) * H_avg' * H2_avg/mean((H_avg(:).^2+H2_avg(:).^2)/2);
% disp(similarity_matrix(:));
