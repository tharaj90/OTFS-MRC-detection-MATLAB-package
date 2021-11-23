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




function [est_bits,ite,x_data] = Algorithm3_low_complexity_detector(N,M,M_mod,no,data_grid,r,H_tf,gs,L_set,omega,decision,init_estimate,n_ite)
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
data_location=reshape(data_grid,N*M,1);
%number of bits per QAM symbol
M_bits=log2(M_mod);
%number of bits per frame
N_bits_perfram = N_syms_perfram*M_bits;
%received delay-time samples 
Y_tilda=reshape(r,M,N);
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

%% Matched Filter Gauss Siedel algorithm
error=zeros(n_ite);
s_est=reshape(X_tilda_est,N*M,1);
delta_r=r;
d=zeros(N*M,1);
for q=1:N*M
    for l=(L_set+1)
        if(data_location(q)==1)
            d(q)=d(q)+abs(gs(l,q+(l-1))).^2;      % time-domain equivalent of equation (59) in [R1]
        end
    end
end
for q=1:N*M         
    for l=(L_set+1)
        if(q>=l)
            delta_r(q)=delta_r(q)-gs(l,q)*s_est(q-(l-1)); % time-domain equivalent of line 3 of Algorithm 2 in [R1]
        end
    end
end
for ite=1:n_ite    
    delta_g=zeros(N*M,1);
    s_est_old=s_est;
    for q=data_index
            for l=(L_set+1)                
                delta_g(q)=delta_g(q)+conj(gs(l,q+(l-1)))*delta_r(q+(l-1));  % time-domain equivalent of line 8 of Algorithm 2 in [R1]              
            end
            s_est(q)=s_est_old(q)+delta_g(q)/d(q);                           % time-domain equivalent of line 9 of Algorithm 2 in [R1]   
            for l=(L_set+1)                                                             % time-domain equivalent of line 11 of Algorithm 2 in [R1]  
                delta_r(q+(l-1))=delta_r(q+(l-1))-gs(l,q+(l-1))*(s_est(q)-s_est_old(q)); % time-domain equivalent of line 12 of Algorithm 2 in [R1]   
            end                                                                          % time-domain equivalent of line 13 of Algorithm 2 in [R1]  
    end
    s_est_old=s_est;
    if(decision==1)
        X_est=reshape(s_est,M,N)*Fn;
        X_tilda_est=((qammod(qamdemod(X_est,M_mod,'gray'),M_mod,'gray').*data_grid)*Fn'); 
        s_est=(1-omega)*s_est+omega*reshape(X_tilda_est,N*M,1);
    end     
    for q=data_index
            for l=(L_set+1)                                                                % time-domain equivalent of line 11 of Algorithm 2 in [R1]         
                delta_r(q+(l-1))=delta_r(q+(l-1))-gs(l,q+(l-1))*(s_est(q)-s_est_old(q));   % time-domain equivalent of line 12 of Algorithm 2 in [R1] 
            end                                                                             % time-domain equivalent of line 13 of Algorithm 2 in [R1] 
    end
    error(ite)=norm(delta_r);
    if(ite>1)
        if(error(ite)>=error(ite-1))                                                        % time-domain equivalent of line 15 of Algorithm 2 in [R1]
            break;
        end
    end    
end
if(n_ite==0)
    ite=0;
end
%% detector output likelihood calculations for turbo decode
X_tilda_est=reshape(s_est,M,N);
X_est=X_tilda_est*Fn;
x_est=reshape(X_est,1,N*M);
x_data=x_est(data_index);
est_bits=reshape(qamdemod(x_data,M_mod,'gray','OutputType','bit'),N_bits_perfram,1);

end
