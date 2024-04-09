function [freq,H] = getH4(freq,...
    x_tx,y_tx,fres_tx,chi_tx,gamma_tx,...
    x_rx,y_rx,fres_rx,chi_rx,gamma_rx,...
    x_env,y_env,fres_env,chi_env,gamma_env,...
    x_ris,y_ris,fres_ris,chi_ris,gamma_ris)

% constants of proportionality and constants of unity value are neglected for better readability

k       = 2*pi*freq;

%% collect dipole parameters

x     = [x_tx x_rx x_env x_ris];
y     = [y_tx y_rx y_env y_ris];
fres  = [fres_tx fres_rx fres_env fres_ris];
chi   = [chi_tx chi_rx chi_env chi_ris];
gamma = [gamma_tx gamma_rx gamma_env gamma_ris];

N_T   = length(x_tx);
N_R   = length(x_rx);
N_E   = length(x_env);
N_RIS = length(x_ris);
N     = N_T+N_R+N_E+N_RIS;
W_freq = zeros(length(freq),N,N);

%% loop over frequencies
for ff=1:length(freq)
    X_diff = zeros(N,N);
    Y_diff = zeros(N,N);
    for l=1:N
        xl_vec = x(l)*ones(1,N);
        yl_vec = y(l)*ones(1,N);
        X_diff(l,:) = x - xl_vec;
        Y_diff(l,:) = y - yl_vec;
    end
    
    BesselInp= k(ff)*sqrt(X_diff.^2+Y_diff.^2);
    % BesselInp= sqrt(X_diff.^2+Y_diff.^2);
    BesselOut = besselh(0,2,BesselInp);
    W_freq(ff,:,:) = 1i*(k(ff)^2/4)*BesselOut;
end
for ff=1:length(freq)
    

    %% Assemble W
    % X_diff = zeros(N,N);
    % Y_diff = zeros(N,N);
    % for l=1:N
    %     xl_vec = x(l)*ones(1,N);
    %     yl_vec = y(l)*ones(1,N);
    %     X_diff(l,:) = x - xl_vec;
    %     Y_diff(l,:) = y - yl_vec;
    % end
    
    % BesselInp= k(ff)*sqrt(X_diff.^2+Y_diff.^2);
    % BesselInp= sqrt(X_diff.^2+Y_diff.^2);
    % BesselOut = besselh(0,2,BesselInp);
    % off-diagonal entries of W are the negative free-space Green's functions between ith and jth dipoles
    % W = 1i*(k(ff)^2/4)*BesselOut;
    W2 = squeeze(W_freq(ff,:,:));
    
    for ii=1:N
        % diagonal entries of W are the inverse polarizabilities
        inv_alpha = ((2*pi*fres(ii))^2-(2*pi*freq(ff))^2)...
                    /(chi(ii)^2) + 1i*(((k(ff)^2)/4) + 2*pi*freq(ff)*gamma(ii)/(chi(ii)^2));
        % inv_alpha = 0;
        % W(ii,ii)  = inv_alpha;
        W2(ii,ii) = inv_alpha;
    end
    
    %% Invert W and extract H
    % load('W_mat.mat','W_mat');
    % Winv = W^(-1);
    Winv2 = W2^(-1);
    % V = diag(diag(W))*Winv;
    V2 = diag(diag(W2))*Winv2;
    % H(ff,:,:) = V((N_T+1):(N_T+N_R),1:N_T);
    H(ff,:,:) = V2((N_T+1):(N_T+N_R),1:N_T);
    % if isequal(H,H2)
    %     disp("equal")
    % else
    %     disp("not equal")
    % end
end

end

