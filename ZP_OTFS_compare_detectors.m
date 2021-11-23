%
% Copyright (c) 2021, Tharaj Thaj, Emanuele Viterbo, and  Yi Hong, Monash University
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
% 3. The reference listed below should be cited if the corresponding codes are used for
%   publication..
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%    - Latest version of this code may be downloaded from: https://ecse.monash.edu/staff/eviterbo/
%    - Freely distributed for educational and research purposes
%References

%  [R1]. T. Thaj and E. Viterbo, "Low Complexity Iterative Rake Decision Feedback Equalizer for Zero-Padded OTFS Systems," in IEEE Transactions on Vehicular Technology, vol. 69, no. 12, pp. 15606-15622, Dec. 2020, doi: 10.1109/TVT.2020.3044276.
%  [R2]. T. Thaj and E. Viterbo,``Low Complexity Iterative Rake Detector for Orthogonal Time Frequency Space Modulation’’ 2020 IEEE Wireless Communications and Networking Conference (WCNC), 2020, pp. 1-6, doi: 10.1109/WCNC45663.2020.9120526.
%  [R3]. Y. Hong, T. Thaj, E. Viterbo, ``Delay-Doppler Communications: Principles and Applications'', Academic Press, 2022, ISBN:9780323850285






close all
clear all
rng(1)
%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time
N = 16;
% M: number of subcarriers in frequency
M = 64;
% M_mod: size of QAM constellation
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));

%% delay-Doppler grid symbol placement
% ZP length  should be set to greater than or equal to maximum value
% of delay_taps
length_ZP = M/16;
% data positions of OTFS delay-Doppler domain data symbols  in the 2-D grid
M_data = M-length_ZP;
data_grid=zeros(M,N);
data_grid(1:M_data,1:N)=1;
% number of symbols per frame
N_syms_perfram = sum(sum(data_grid));
% number of bits per frame
N_bits_perfram = N_syms_perfram*M_bits;



% Time and frequency resources
car_fre = 4*10^9;% Carrier frequency
delta_f = 15*10^3; % subcarrier spacing: 15 KHz
T = 1/delta_f; %one time symbol duration in OTFS frame

% SNR and variance of the noise
% SNR = P/\sigma^2; P: avg. power of albhabet transmitted
SNR_dB = 5:5:30;
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;



%% Initializing simulation error count variables

N_fram = 1000;
Initialize_error_count_variables;

%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix
current_frame_number=zeros(1,length(SNR_dB));

for iesn0 = 1:length(SNR_dB)
    for ifram = 1:N_fram
        current_frame_number(iesn0)=ifram;
        %% random input bits generation%%%%%
        trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);
        %%2D QAM symbols generation %%%%%%%%
        data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');
        X = Generate_2D_data_grid(N,M,data,data_grid);
        
        
        %% OTFS modulation%%%%
        X_tilda=X*Fn';
        s = reshape(X_tilda,N*M,1);
        
        
        %% OTFS channel generation%%%% The user can use either the synthetic channel model or the 3GPP channel by uncommenting the corresonding piece of code. 
        %% synthetic channel model with equal power paths with delays [0,l_max] and Dopplers [-k_max,k_max]
%         taps=4;
%         l_max=delay_spread;
%         k_max=4;
%         chan_coef=1/sqrt(2)*(randn(1,taps)+1i.*randn(1,taps));
%         delay_taps=randi(taps,[1,l_max+1]);  
%         delay_taps=sort(delay_taps-min(delay_taps));  %% random delay shifts in the range [0,l_max] 
%         Doppler_taps=k_max-2*k_max*rand(1,taps);   %% uniform Doppler profile [-k_max,k_max]
%         L_set=unique(delay_taps);
        %% channel model following 3GPP standard
        max_speed=500;  % km/hr
        [chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(N,M,car_fre,delta_f,T,max_speed);
        L_set=unique(delay_taps);

        
              
        %% channel output%%%%%
        [G,gs]=Gen_time_domain_channel(N,M,taps,delay_taps,Doppler_taps,chan_coef);
        [H,H_tilda,P]= Gen_DD_and_DT_channel_matrices(N,M,G,Fn);
        r=zeros(N*M,1);
        noise= sqrt(sigma_2(iesn0)/2)*(randn(size(s)) + 1i*randn(size(s)));
        l_max=max(delay_taps);
        for q=1:N*M
            for l=(L_set+1)
                if(q>=l)
                    r(q)=r(q)+gs(l,q)*s(q-l+1);
                end
            end
        end
        r=r+noise;
        
        %% OTFS demodulation%%%%
        Y_tilda=reshape(r,M,N);
        Y = Y_tilda*Fn;
        y=reshape(Y.',N*M,1);
        
        %% test: the received time domain signal can be generated element in the matrix form (using r=G.s+noise).
        %         r_test=G*s+noise;
        %         test_delay_time_matrix_error=norm(r_test-r)
        %% test: the received DD domain signal can be generated in the matrix form (using y=H.x+noise).
        %         noise_DD=kron(eye(M),Fn)*P'*noise;
        %         x=reshape(X.',N*M,1);
        %         y_test=H*x_vec+noise_DD;
        %         text_delay_Doppler_matrix_error=norm(y_test-y)
        
        %% Generate delay-time and delay-Doppler channel vectors from the time domain channel.
%         [nu_ml_tilda]=Gen_delay_time_channel_vectors(N,M,l_max,gs);
        [nu_ml_tilda,nu_ml,K_ml]=Gen_DT_and_DD_channel_vectors(N,M,L_set,gs); 
        
        %% Generate block-wise time-frequency domain channel
        [H_tf]=Generate_time_frequency_channel_ZP(N,M,gs,L_set);
        
        %% Different detection methods
        
        n_ite_MRC=15; % maximum number of MRC detector iterations  (Algorithm 2 in [R1])
        n_ite_algo3=15; % maximum number of matched filtered Guass Seidel (Algorithm 3 in [R1])
        n_ite_MPA=15; % maximum number of MPA detector iterations
        %damping parameter - reducing omega improves error performance at the cost of increased detector iterations
        omega=1;
        if(M_mod>=64)
            omega=0.25;     % set omega to a smaller value (for example: 0.05) for modulation orders greater than 64-QAM
        end
        decision=1; %1-hard decision, 0-soft decision
        init_estimate=1; %1-use the TF single tap estimate as the initial estimate for MRC detection, 0-initialize the symbol estimates to 0 at the start of MRC iteration
        %(Note: it is recommended to set init_estimate to 0 for higher order modulation schemes like 64-QAM as the single tap equalizer estimate may be less accurate)
        
        %MRC detectors in [R1]
        
        [est_info_bits_MRC,det_iters_MRC,data_MRC] = MRC_delay_time_detector(N,M,M_data,M_mod,sigma_2(iesn0),data_grid,r,H_tf,nu_ml_tilda,L_set,omega,decision,init_estimate,n_ite_MRC);
        [est_info_bits_Algo1,det_iters_Algo1,data_Algo1]= Algorithm1_detector(N,M,M_data,M_mod,sigma_2(iesn0),data_grid,r,H_tf,L_set,omega,decision,init_estimate,n_ite_MRC,K_ml,Y);
        [est_info_bits_Algo3,det_iters_Algo3,data_Algo3] = Algorithm3_detector(N,M,M_mod,sigma_2(iesn0),data_grid,r,G,omega,decision,init_estimate,n_ite_algo3);  
        [est_info_bits_Algo3_low_complexity,det_iters_Algo3_low_complexity,data_Algo3_low_complexity] = Algorithm3_low_complexity_detector(N,M,M_mod,sigma_2(iesn0),data_grid,r,H_tf,gs,L_set,omega,decision,init_estimate,n_ite_algo3);
        
        % Other detectors in the literature
        
        [est_info_bits_MPA,det_iters_MPA,data_MPA] = MPA_detector(N,M,M_mod,sigma_2(iesn0),data_grid,y,H,n_ite_MPA);
        [est_info_bits_1tap,data_1tap] = TF_single_tap_equalizer(N,M,M_mod,sigma_2(iesn0),data_grid,Y,H_tf);
        [est_info_bits_LMMSE,data_LMMSE] = Block_LMMSE_detector(N,M,M_mod,sigma_2(iesn0),data_grid,r,gs,L_set);
        
        
        %% errors count%%%%%
        count_errors_per_frame;
        
        
        %% DISP error performance details
        display_errors_per_frame;
        
        
        
    end
    
end

figure(1)
semilogy(SNR_dB,avg_ber_Algo1,'-x','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,avg_ber_MRC,'-x','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,avg_ber_Algo3,'-x','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,avg_ber_Algo3_low_complexity,'-x','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,avg_ber_MPA,'-x','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,avg_ber_1tap,'-x','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,avg_ber_LMMSE,'-s','LineWidth',2,'MarkerSize',8)
legend('MRC (Algorithm 1 in [R1])','MRC (Algorithm 2 in [R1])','Algorithm 3 in [R1]','low-complexity implementation of Algorithm 3 in [R1]','MPA','time-freq single tap in [R3]','block-wise time-domain LMMSE in [R3]')
grid on
xlabel('SNR(dB)')
ylabel('BER')

