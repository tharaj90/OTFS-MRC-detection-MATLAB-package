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


function [Gn_block_matrix,Tn_block_matrix,zn_block_vector,H_t_f]=Generate_Algo3_GS_iteration_matrices(N,M,G,r)

% Generate time-domain GS iteration matrices for low complexity iterative detection
Gn_block_matrix=zeros(M,M,N);
Tn_block_matrix=zeros(M,M,N);
Qn_block_matrix=zeros(M,M,N);
zn_block_vector=zeros(M,N);

H_t_f=zeros(N,M); % Time-frequency single tap channel matrix
Fn=dftmtx(M);
Fn=Fn./norm(Fn);
for n=1:N
    rn=r((n-1)*M+1:n*M);
    Gn_block_matrix(:,:,n)=G((n-1)*M+1:n*M,(n-1)*M+1:n*M);
    Gn=Gn_block_matrix(:,:,n);
    H_t_f(n,1:M)=diag(Fn*Gn*Fn').';  % Generate time-frequency channel matrix for low complexity initial estimate in [R1,R3]
    Rn=Gn'*Gn;    
    Dn=diag(diag(Rn));
    Ln=tril(Rn,-1);
    Un=triu(Rn,1);
    Qn=(Dn+Ln)^(-1);
    Tn=Qn*Un;
    Tn_block_matrix(:,:,n)=Tn;
    Qn_block_matrix(:,:,n)=Qn;        
    zn_block_vector(:,n)=Qn*Gn'*rn;
end
end