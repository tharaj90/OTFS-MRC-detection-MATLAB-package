est_info_bits_MRC=zeros(N_bits_perfram,1);
est_info_bits_Algo1=zeros(N_bits_perfram,1);
est_info_bits_Algo3=zeros(N_bits_perfram,1);
est_info_bits_Algo3_low_complexity=zeros(N_bits_perfram,1);
est_info_bits_MPA=zeros(N_bits_perfram,1);
est_info_bits_1tap=zeros(N_bits_perfram,1);
est_info_bits_LMMSE=zeros(N_bits_perfram,1);

err_ber_MRC = zeros(1,length(SNR_dB));
err_ber_Algo1 = zeros(1,length(SNR_dB));
err_ber_Algo3 = zeros(1,length(SNR_dB));
err_ber_Algo3_low_complexity = zeros(1,length(SNR_dB));
err_ber_MPA = zeros(1,length(SNR_dB));
err_ber_1tap = zeros(1,length(SNR_dB));
err_ber_LMMSE = zeros(1,length(SNR_dB));

avg_ber_MRC=zeros(1,length(SNR_dB));
avg_ber_Algo1=zeros(1,length(SNR_dB));
avg_ber_Algo3=zeros(1,length(SNR_dB));
avg_ber_Algo3_low_complexity=zeros(1,length(SNR_dB));
avg_ber_MPA=zeros(1,length(SNR_dB));
avg_ber_1tap=zeros(1,length(SNR_dB));
avg_ber_LMMSE=zeros(1,length(SNR_dB));

det_iters_MRC=0;
no_of_detetor_iterations_MRC= zeros(length(SNR_dB),1);
avg_no_of_iterations_MRC=zeros(1,length(SNR_dB));

det_iters_Algo1=0;
no_of_detetor_iterations_Algo1= zeros(length(SNR_dB),1);
avg_no_of_iterations_Algo1=zeros(1,length(SNR_dB));

det_iters_Algo3=0;
no_of_detetor_iterations_Algo3= zeros(length(SNR_dB),1);
avg_no_of_iterations_Algo3=zeros(1,length(SNR_dB));

det_iters_Algo3_low_complexity=0;
no_of_detetor_iterations_Algo3_low_complexity= zeros(length(SNR_dB),1);
avg_no_of_iterations_Algo3_low_complexity=zeros(1,length(SNR_dB));

det_iters_MPA=0;
no_of_detetor_iterations_MPA= zeros(length(SNR_dB),1);
avg_no_of_iterations_MPA=zeros(1,length(SNR_dB));
