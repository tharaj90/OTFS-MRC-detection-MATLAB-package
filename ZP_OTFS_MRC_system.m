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
%   publication.
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

%  [R1]. T. Thaj and E. Viterbo, "Low Complexity Iterative Rake Decision Feedback Equalizer for Zero-Padded OTFS Systems," in IEEE Transactions on Vehicular Technology, vol. 69, no. 12, pp. 15606-15622, Dec. 2020, doi: 10.1109/TVT.2020.3044276.
%  [R2]. T. Thaj and E. Viterbo,``Low Complexity Iterative Rake Detector for Orthogonal Time Frequency Space Modulation’’ 2020 IEEE Wireless Communications and Networking Conference (WCNC), 2020, pp. 1-6, doi: 10.1109/WCNC45663.2020.9120526.
%  [R3]. Y. Hong, T. Thaj, E. Viterbo, ``Delay-Doppler Communications: Principles and Applications'', Academic Press, 2022, ISBN:9780323850285






close all
clear all
rng(1)
%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time
N = 64;
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
SNR_dB = 10:2.5:20;
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;



%% Initializing simulation error count variables

N_fram = 1000;
est_info_bits_MRC=zeros(N_bits_perfram,1);
err_ber_MRC = zeros(1,length(SNR_dB));
no_of_detetor_iterations_MRC= zeros(length(SNR_dB),1);
avg_no_of_iterations_MRC=zeros(1,length(SNR_dB));


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
        X_tilda=X*Fn';                     %equation (2) in [R1]
        s = reshape(X_tilda,N*M,1);        %equation (4) in [R1]
        
        
        %% OTFS channel generation%%%% The user can use either the synthetic channel model or the 3GPP channel by uncommenting the corresonding piece of code.        
        %% synthetic channel model with equal power paths with delays [0,l_max] and Dopplers [-k_max,k_max]
%                 taps=4;
%                 l_max=delay_spread;
%                 k_max=4;
%                 chan_coef=1/sqrt(2)*(randn(1,taps)+1i.*randn(1,taps));
%                 delay_taps=randi(taps,[1,l_max+1]);  
%                 delay_taps=sort(delay_taps-min(delay_taps));  %% random delay shifts in the range [0,l_max] 
%                 Doppler_taps=k_max-2*k_max*rand(1,taps);   %% uniform Doppler profile [-k_max,k_max]
        
        % channel model following 3GPP standard
        max_speed=500;  % km/hr
        [chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(N,M,car_fre,delta_f,T,max_speed);
        L_set=unique(delay_taps);       
        gs=Gen_discrete_time_channel(N,M,taps,delay_taps,Doppler_taps,chan_coef);  %equation (14) in [R1]
        
        %% channel output%%%%%             
        r=zeros(N*M,1);        
        l_max=max(delay_taps);
        for q=1:N*M
            for l=(L_set+1)
                if(q>=l)
                    r(q)=r(q)+gs(l,q)*s(q-l+1);  %equation (18) in [R1]
                end
            end
        end
        noise= sqrt(sigma_2(iesn0)/2)*(randn(size(s)) + 1i*randn(size(s)));
        r=r+noise;
                
        %% OTFS demodulation%%%%
        Y_tilda=reshape(r,M,N);           %equation (5) in [R1]
        Y = Y_tilda*Fn;                   %equation (6) in [R1]
               
        %% Generate delay-time channel vectors from the time domain channel.
        [nu_ml_tilda]=Gen_delay_time_channel_vectors(N,M,l_max,gs);  % equation (42) in [R1]
        
        %% Generate block-wise time-frequency domain channel
        [H_tf]=Generate_time_frequency_channel_ZP(N,M,gs,L_set);  
                
        %% MRC delay-time detection in [R1,R2,R3]
        
        n_ite_MRC=50; % maximum number of MRC detector iterations
        %damping parameter - reducing omega improves error performance at the cost of increased detector iterations
        omega=1;
        if(M_mod==64)
            omega=0.25;     % set omega to a smaller value (for example: 0.05) for modulation orders greater than 64-QAM
        end
        decision=1;         %1-hard decision, 0-soft decision
        init_estimate=1;    %1-use the TF single tap estimate as the initial estimate for MRC detection, 0-initialize the symbol estimates to 0 at the start of MRC iteration 
        %(Note: it is recommended to set init_estimate to 0 for higher order modulation schemes like 64-QAM or 256-QAM as the single tap equalizer estimate may be less accurate)
        
        [est_info_bits_MRC,det_iters_MRC,data_MRC] = MRC_delay_time_detector(N,M,M_data,M_mod,sigma_2(iesn0),data_grid,r,H_tf,nu_ml_tilda,L_set,omega,decision,init_estimate,n_ite_MRC);               

        %% errors count%%%%%
        errors_MRC = sum(xor(est_info_bits_MRC,trans_info_bit));                
        err_ber_MRC(iesn0) = err_ber_MRC(iesn0) + errors_MRC;        
        avg_ber_MRC(iesn0)=err_ber_MRC(iesn0).'/length(trans_info_bit)/ifram;
        
               
        %%  iterations count        
        no_of_detetor_iterations_MRC(iesn0)=no_of_detetor_iterations_MRC(iesn0)+det_iters_MRC;
        avg_no_of_iterations_MRC(iesn0)=no_of_detetor_iterations_MRC(iesn0)/ifram;
        
               
        %% DISP error performance details        
         clc
        fprintf('ZP-OTFS-(N,M,QAM size)');disp([N,M,M_mod]);
        display(current_frame_number,'Number of frames');
        display(avg_ber_MRC,'Average BER - Delay-time domain Maximal Ratio Combining (MRC)');        
        display(avg_no_of_iterations_MRC,'Average number of iterations for the delay-time domain MRC detector');       
               
    end
    
end

figure(1)
semilogy(SNR_dB,avg_ber_MRC,'-x','LineWidth',2,'MarkerSize',8)
legend('MRC (Algorithm 2 in [R1])')
grid on
xlabel('SNR(dB)')
ylabel('BER')
