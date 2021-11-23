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

function [est_bits,ite,x_data] = Algorithm3_detector(N,M,M_mod,no,data_grid,r,G,omega,decision,init_estimate,n_ite)
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
[Gn_block_matrix,Tn_block_matrix,zn_block_vector,H_tf]=Generate_Algo3_GS_iteration_matrices(N,M,G,r);
%% initial time-frequency low complexity estimate assuming ideal pulses
if(init_estimate==1)
    Y_tf=fft(Y_tilda).'; % delay-time to frequency-time domain                                                                      % equation (63) in [R1]                                   
    X_tf=conj(H_tf).*Y_tf./(H_tf.*conj(H_tf)+no); % single tap equalizer                                                            % equation (64) in [R1]
    X_est = ifft(X_tf.')*Fn; % SFFT                                                                                                 % equation (65) in [R1]
    X_est=qammod(qamdemod(X_est,M_mod,'gray'),M_mod,'gray');
    X_est=X_est.*data_grid;
    X_tilda_est=X_est*Fn';
else
    X_est=zeros(M,N);
    X_tilda_est=X_est*Fn';
end
X_tilda_est=X_tilda_est.*data_grid;

               
%% Matched Filter Gauss Siedel algorithm (Algorithm 3 in [R1])
error=zeros(n_ite);
x_soft=zeros(M,N);

for ite=1:n_ite
    for i=1:N
        Tn=Tn_block_matrix(:,:,i);
        zn=zn_block_vector(:,i);
        sn_prev=X_tilda_est(:,i);
        sn_next=-Tn*sn_prev+zn;                                              %step 12 of Algorithm 3 in [R1]
        x_soft(:,i)=sn_next.*data_grid(:,i);
        X_tilda_est(:,i)=(x_soft(:,i));        
    end
    if(decision==1)
        x_m=x_soft*Fn;
        X_tilda_est=(1-omega)*X_tilda_est+omega*((qammod(qamdemod(x_m,M_mod,'gray'),M_mod,'gray').*data_grid)*Fn');  %equation (50) in [R1] (or equivalently equation (27) in [R2])
    end    
    for i=1:N
        Gn=Gn_block_matrix(:,:,i);
        Y_tilda_est(1:M,i)=Gn*(X_tilda_est(1:M,i));
    end
    error(ite)=sum(sum(abs(Y_tilda_est-Y_tilda).^2))./N_syms_perfram;        %step 14 of Algorithm 3 in [R1]
    if(ite>1)
        if(error(ite)>=error(ite-1))
            break;
        end
    end    
end
if(n_ite==0)
    ite=0;
end
%% detector output likelihood calculations for turbo decode
X_est=X_tilda_est*Fn;
x_est=reshape(X_est,1,N*M);
x_data=x_est(data_index);
est_bits=reshape(qamdemod(x_data,M_mod,'gray','OutputType','bit'),N_bits_perfram,1);

end
