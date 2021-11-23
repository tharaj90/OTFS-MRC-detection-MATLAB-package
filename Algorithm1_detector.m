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





function [est_bits,ite,x_data] = Algorithm1_detector(N,M,M_data,M_mod,no,data_grid,r,H_tf,L_set,omega,decision,init_estimate,n_ite_MRC,K_ml,Y)
%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix
%% Initial assignments
%Number of symbols per frame
N_syms_perfram=sum(sum((data_grid>0)));
%Arranging the delay-Doppler grid symbols into an array
data_array=reshape(data_grid,1,N*M);
%finding position of data symbols in the array
[~,data_index]=find(data_array>0);
%number of bits per QAM symbol
M_bits=log2(M_mod);
%number of bits per frame
N_bits_perfram = N_syms_perfram*M_bits;
%received delay-time samples
Y_tilda=reshape(r,M,N);
M_prime=M_data;
L_set=L_set+1; % Since matlab cannot work with zero indices
%% initial estimate using single tap TF equalizer
if(init_estimate==1)
    Y_tf=fft(Y_tilda).'; % delay-time to frequency-time domain                                                                      % equation (63) in [R1]
    X_tf=conj(H_tf).*Y_tf./(H_tf.*conj(H_tf)+no); % single tap equalizer                                                            % equation (64) in [R1]
    X_est = ifft(X_tf.')*Fn; % SFFT                                                                                                 % equation (65) in [R1]
    X_est=qammod(qamdemod(X_est,M_mod,'gray'),M_mod,'gray');
    X_est=X_est.*data_grid;
else
    X_est=zeros(M,N);
end
x_m_hat=X_est.';

%% MRC detector    %% Algorithm 1 in [R1]
%% initial computation
D_m=zeros(N,N,M);
D_m_inv=zeros(N,N,M);
y_m=Y.';
% RNPI_error=y_m;
for m=1:M_prime
    for l=L_set
        D_m(:,:,m)=D_m(:,:,m)+K_ml(:,:,m,l)'*K_ml(:,:,m,l);                                                             % equation (47) in [R1]
    end
    D_m_inv(:,:,m)=inv(D_m(:,:,m));
end
c_m=x_m_hat;
b_m_l=zeros(N,M,length(L_set));

%% iterative computation
for ite=1:n_ite_MRC
    g_m=zeros(N,M);
    RNPI_error=y_m;
    for m=1:M_prime
        for l=L_set
            b_m_l(:,m,l)=y_m(:,m+l-1);
            for p=L_set
                if(l~=p)
                    if(m+(l-p)>0)
                        b_m_l(:,m,l)=b_m_l(:,m,l)-K_ml(:,:,m+l-1,p)*x_m_hat(:,m+(l-p));                                % Line 5 of Algorithm 1 in [R1] 
                    end
                end
            end
            g_m(:,m)=g_m(:,m)+K_ml(:,:,m+(l-1),l)'*b_m_l(:,m,l);                                                       % Line 7 of Algorithm 1 in [R1] 
        end
        c_m(:,m)=D_m_inv(:,:,m)*g_m(:,m);                                                                              % Line 8 of Algorithm 1 in [R1] 
        if(decision==1)
            x_m_hat(:,m)=(1-omega)*c_m(:,m)+omega*qammod(qamdemod((c_m(:,m)),M_mod,'gray'),M_mod,'gray');              % Line 9 of Algorithm 1 in [R1] 
        else
            x_m_hat(:,m)=c_m(:,m);
        end
        
        for l=L_set
            if(m>=l)
                RNPI_error(:,m)=RNPI_error(:,m)-reshape(K_ml(:,:,m,l),N,N)*x_m_hat(:,m-(l-1));                         % equation (50) in [R1] - residual interference
            end
        end
    end
    %% convergence criteria
    error(ite)=norm(RNPI_error);
    if(ite>1)
        if(error(ite)>=error(ite-1))
            break;
        end
    end
end
if(n_ite_MRC==0)
    ite=0;
end
%% detector output bits
X_est=x_m_hat.';
x_est=reshape(X_est,1,N*M);
x_data=x_est(data_index);
est_bits=reshape(qamdemod(x_data,M_mod,'gray','OutputType','bit'),N_bits_perfram,1);

end
