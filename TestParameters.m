function [Validity,N_T,N_R,N_E] = TestParameters(x_tx,y_tx,fres_tx,chi_tx,gamma_tx, ...
                                                 x_rx,y_rx,fres_rx,chi_rx,gamma_rx, ...
                                                 x_env,y_env,fres_env,chi_env,gamma_env, ...
                                                 chi_ris,gamma_ris,N_RIS)
    Validity = true;
    N_T = 0;
    N_R = 0;
    N_E = 0;
    %% Testing Parameter Validity
    if length(x_tx)~=length(y_tx)
        error('Error: x_tx and y_tx do not have the same length.');
        Validity = false;
    else
        N_T = length(x_tx);
    end
    if length(fres_tx)~=N_T
        error('Error: x_tx and fres_tx do not have the same length.');
        Validity = false;
    end
    if length(chi_tx)~=N_T
        error('Error: x_tx and chi_tx do not have the same length.');
        Validity = false;
    end
    if length(gamma_tx)~=N_T
        error('Error: x_tx and gamma_tx do not have the same length.');
        Validity = false;
    end
    if length(x_rx)~=length(y_rx)
        error('Error: x_rx and y_rx do not have the same length.');
        Validity = false;
    else
        N_R = length(x_rx);
    end
    if length(fres_rx)~=N_R
        error('Error: x_rx and fres_rx do not have the same length.');
        Validity = false;
    end
    if length(chi_rx)~=N_R
        error('Error: x_rx and chi_rx do not have the same length.');
        Validity = false;
    end
    if length(gamma_rx)~=N_R
        error('Error: x_rx and gamma_rx do not have the same length.');
        Validity = false;
    end
    if length(x_env)~=length(y_env)
        error('Error: x_env and y_env do not have the same length.');
        Validity = false;
    else
        N_E = length(x_env);
    end
    if length(fres_env)~=N_E
        error('Error: x_env and fres_env do not have the same length.');
        Validity = false;
    end
    if length(chi_env)~=N_E
        error('Error: x_env and chi_env do not have the same length.');
        Validity = false;
    end
    if length(gamma_env)~=N_E
        error('Error: x_env and gamma_env do not have the same length.');
        Validity = false;
    end
    
    if length(chi_ris(1,:))~=N_RIS
        error('Error: x_ris and chi_ris do not have the same length.');
        Validity = false;
    end
    if length(gamma_ris(1,:))~=N_RIS
        error('Error: x_ris and gamma_ris do not have the same length.');
        Validity = false;
    end
end
