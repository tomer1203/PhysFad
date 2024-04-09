function [freq,H] = getH3(freq,W_freq,...
    fres_tx,chi_tx,gamma_tx,...
    fres_rx,chi_rx,gamma_rx,...
    fres_env,chi_env,gamma_env,...
    fres_ris,chi_ris,gamma_ris)

% constants of proportionality and constants of unity value are neglected for better readability

k       = 2*pi*freq;

%% collect dipole parameters


fres  = [fres_tx fres_rx fres_env fres_ris];
chi   = [chi_tx chi_rx chi_env chi_ris];
gamma = [gamma_tx gamma_rx gamma_env gamma_ris];

N_T   = length(fres_tx);
N_R   = length(fres_rx);
N_E   = length(fres_env);
N_RIS = length(fres_ris);
N     = N_T+N_R+N_E+N_RIS;

%% loop over frequencies
% inv_alpha2 = ((2*pi*fres).^2-(2*pi*freq').^2)...
%         ./(chi.^2) + 1i*(((k'.^2)/4) + 2*pi*freq'*gamma./(chi.^2));
for ff=1:length(freq)
    W = squeeze(W_freq(ff,:,:));
    % inv_alpha1 = ((2*pi*fres).^2-(2*pi*freq(ff)).^2)...
    %     ./(chi.^2) + 1i*(((k(ff).^2)/4) + 2*pi*freq(ff)*gamma./(chi.^2));
    % W = W+inv_alpha2(ff,:);

    for ii=1:N
        % diagonal entries of W are the inverse polarizabilities
        inv_alpha = ((2*pi*fres(ii))^2-(2*pi*freq(ff))^2)...
                    /(chi(ii)^2) + 1i*(((k(ff)^2)/4) + 2*pi*freq(ff)*gamma(ii)/(chi(ii)^2));

        W(ii,ii)  = inv_alpha;
        
    end
    
    %% Invert W and extract H
    % load('W_mat.mat','W_mat');
    Winv = W^(-1);
    
    V = diag(diag(W))*Winv;
    
    H(ff,:,:) = V((N_T+1):(N_T+N_R),1:N_T);
    
end

end

