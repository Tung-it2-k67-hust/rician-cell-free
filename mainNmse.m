%Empty workspace and close figures
clear
close all
clc

%Define the range of number of Access Points (APs)
M=60:20:100;

%Number of UEs
K=40; 
%Pilot length
tau_p=5;
%Select length of coherence block
tau_c=200;

%The cell size 1km x 1km
cellRange=1000;

%Compute the noise power
noiseFigure = 7;
B=20e6; %Communication bandwidth
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;
sigma2=db2pow(noiseVariancedBm-30); %(in Watt)

%Uplink transmit power per UE (W)
p=0.2; %200 mW
%Create the power vector for all UEs (The uplink power is the same
%(p)at each UE)
pv=p*ones(1,K);



%Select the number of setups with random AP/UE locations
nbrOfSetups=2;
%Select the number of channel realizations per setup
nbrOfRealizations=5e2;

% Khai báo mảng lưu NMSE trung bình theo từng số lượng AP (M)
NMSE_MMSE_mean = zeros(1, length(M));
NMSE_LMMSE_mean = zeros(1, length(M));
NMSE_LS_mean    = zeros(1, length(M));


%Go through the number of APs
for m=1:length(M)
    % Biến tạm lưu kết quả cho từng setup
    NMSE_MMSE_all = zeros(1, nbrOfSetups);
    NMSE_LMMSE_all = zeros(1, nbrOfSetups);
    NMSE_LS_all = zeros(1, nbrOfSetups);

    %Deploy APs randomly
    APpositions=cellRange*(rand(M(m),1) + 1i*rand(M(m),1));
    
    %For single layer decoding set all large-scale fading coefficients to 1
    A_singleLayer=reshape(repmat(eye(M(m)),1,K),M(m),M(m),K);
    
    %Go through all setups
    for n=1:nbrOfSetups
      
       %Deploy UEs and generate the covariance and mean matrices
       [R,HMeanWithoutPhase,channelGain] = functionCellFreeSetup( M(m),K,cellRange,APpositions,sigma2,1);
       
       
       %Create channel generations for each UE-AP pair
       [H,HMean] = functionChannelGeneration( R,HMeanWithoutPhase,M(m),nbrOfRealizations,K );
       
       %Pilot allocation 
       [Pset] = functionPilotAllocation( R,HMeanWithoutPhase,A_singleLayer,K,M(m),pv,tau_p);
       
       %Second Layer Decoding
       %Generate Large-scale fading coefficients for phase-aware MMSE,
       %Linear MMSE and LS estimators
       %[ A_MMSE,A_LMMSE,A_LS] = functionLSFD( R,HMeanWithoutPhase,M(m),K,pv,tau_p,Pset);
       
       %Channel estimation 
       %Channel estimation with phase-aware MMSE estimator
       [Hhat_MMSE] =functionCellFreeMMSE(R,HMean,H,nbrOfRealizations,M(m),K,pv,tau_p,Pset);
       %Channel estimation with LMMSE estimator
       [Hhat_LMMSE] =functionCellFreeLMMSE(R,HMeanWithoutPhase,H,nbrOfRealizations,M(m),K,pv,tau_p,Pset);
       %Channel estimation with LS estimator
       [Hhat_LS] = functionCellFreeLS( H,nbrOfRealizations,M(m),K,pv,tau_p,Pset);
       
  

        %-------------------------------------------------------------
        % Monte-Carlo NMSE computation for each estimator
        %-------------------------------------------------------------
        [NMSE_MMSE, NMSE_MMSE_dB]   = functionMonteCarloNMSE_UL(Hhat_MMSE, H, [], tau_c, tau_p, nbrOfRealizations, M(m), K, pv);
        [NMSE_LMMSE, NMSE_LMMSE_dB] = functionMonteCarloNMSE_UL(Hhat_LMMSE, H, [], tau_c, tau_p, nbrOfRealizations, M(m), K, pv);
        [NMSE_LS, NMSE_LS_dB]       = functionMonteCarloNMSE_UL(Hhat_LS, H, [], tau_c, tau_p, nbrOfRealizations, M(m), K, pv);
        fprintf('\n===== M = %d | Setup %d =====\n', M(m), n);
        fprintf('--- MMSE ---\n');
        disp(NMSE_MMSE_dB.');  % hiển thị theo hàng, dB
        fprintf('--- LMMSE ---\n');
        disp(NMSE_LMMSE_dB.');
        fprintf('--- LS ---\n');
        disp(NMSE_LS_dB.');

        %-------------------------------------------------------------
        % In ra kết quả trung bình NMSE (theo dòng)
        %-------------------------------------------------------------
        % fprintf('\n[M = %d | Setup %d/%d]\n', M(m), n, nbrOfSetups);
        % fprintf('  ➤ NMSE_MMSE   = %.4f (%.2f dB)\n', mean(NMSE_MMSE), mean(NMSE_MMSE_dB));
        % fprintf('  ➤ NMSE_LMMSE  = %.4f (%.2f dB)\n', mean(NMSE_LMMSE), mean(NMSE_LMMSE_dB));
        % fprintf('  ➤ NMSE_LS     = %.4f (%.2f dB)\n', mean(NMSE_LS), mean(NMSE_LS_dB));
        NMSE_MMSE_all(n) = mean(NMSE_MMSE_dB);
        NMSE_LMMSE_all(n) = mean(NMSE_LMMSE_dB);
        NMSE_LS_all(n) = mean(NMSE_LS_dB);
    end
    NMSE_MMSE_mean(m) = mean(NMSE_MMSE_all);
    NMSE_LMMSE_mean(m) = mean(NMSE_LMMSE_all);
    NMSE_LS_mean(m) = mean(NMSE_LS_all);
end
figure;
plot(M, NMSE_MMSE_mean, '-o', 'LineWidth', 1.5); hold on;
plot(M, NMSE_LMMSE_mean, '-s', 'LineWidth', 1.5);
plot(M, NMSE_LS_mean, '-^', 'LineWidth', 1.5);

xlabel('Number of Access Points (M)');
ylabel('Average NMSE [dB]');
title('NMSE performance vs Number of APs');
legend('MMSE', 'LMMSE', 'LS', 'Location', 'best');
grid on;
grid on;
