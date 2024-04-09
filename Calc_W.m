function W = Calc_W(freq,...
    x_tx,y_tx,...
    x_rx,y_rx,...
    x_env,y_env,...
    x_ris,y_ris)

% constants of proportionality and constants of unity value are neglected for better readability

k       = 2*pi*freq;

%% collect dipole parameters

x     = [x_tx x_rx x_env x_ris];
y     = [y_tx y_rx y_env y_ris];


N_T   = length(x_tx);
N_R   = length(x_rx);
N_E   = length(x_env);
N_RIS = length(x_ris);
N     = N_T+N_R+N_E+N_RIS;
W = zeros(length(freq),N,N);
%% loop over frequencies
X_diff = zeros(N,N);
Y_diff = zeros(N,N);
for l=1:N
    xl_vec = x(l)*ones(1,N);
    yl_vec = y(l)*ones(1,N);
    X_diff(l,:) = x - xl_vec;
    Y_diff(l,:) = y - yl_vec;
end
for ff=1:length(freq)
    
%     disp(['Currently evaluating frequency point ',num2str(ff),' / ',num2str(length(freq)),'.']);
    BesselInp= k(ff)*sqrt(X_diff.^2+Y_diff.^2);
    % BesselInp= sqrt(X_diff.^2+Y_diff.^2);
    BesselOut = besselh(0,2,BesselInp);
    % jbessel = besselj(0,BesselInp);
    % epsilon = 0.000000001;
    % ybessel = (besselj(epsilon,BesselInp)*cos(epsilon*pi)-besselj(-epsilon,BesselInp))/sin(epsilon*pi);
    % BesselOut = besselj(0,BesselInp) - 1i*ybessel;

    
    %% Assemble W
    
    % off-diagonal entries of W are the negative free-space Green's functions between ith and jth dipoles
    W(ff,:,:) = 1i*(k(ff)^2/4)*BesselOut;
    
    
end

end

