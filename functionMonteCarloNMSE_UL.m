function [NMSE, NMSE_dB] = functionMonteCarloNMSE_UL(Hhat, H, A, tau_c, tau_p, nbrOfRealizations, M, K, p)
%--------------------------------------------------------------------------
% Computes the normalized mean square error (NMSE) of the estimated channels
% via Monte Carlo simulations.
%
% INPUT:
% Hhat              = M x nbrOfRealizations x K matrix of estimated channels
% H                 = M x nbrOfRealizations x K matrix of true channels
% A                 = M x M x K matrix of LSFD coefficients (not used here,
%                     kept for structural consistency with SE function)
% tau_c, tau_p      = coherence block length and pilot length (unused)
% nbrOfRealizations = number of channel realizations
% M                 = number of APs
% K                 = number of UEs
% p                 = 1 x K vector of UE uplink powers (unused here)
%
% OUTPUT:
% NMSE              = K x 1 vector of normalized MSE per UE
% NMSE_dB           = NMSE in dB (10*log10)
%
%--------------------------------------------------------------------------

NMSE = zeros(K,1);

for k = 1:K
    
    % True and estimated channel realizations for UE k
    H_true = H(:,:,k);
    H_est  = Hhat(:,:,k);
    
    % Mean-square error across all APs and realizations
    mse = mean( abs(H_est(:) - H_true(:)).^2 );
    
    % Average true channel power
    powerH = mean( abs(H_true(:)).^2 );
    
    % Normalized MSE
    NMSE(k) = mse / powerH;
end

% Convert to dB for plotting or comparison
NMSE_dB = 10*log10(NMSE);


end
