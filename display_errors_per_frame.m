clc
fprintf('ZP-OTFS-(N,M,QAM size)');disp([N,M,M_mod]);
display(current_frame_number,'Number of frames');
display(SNR_dB,'SNR (dB)');
display(avg_ber_MRC,'Average BER - Delay-time domain Maximal Ratio Combining (MRC)');
display(avg_ber_MRC,'Average BER - Algorithm 1 in [R1]');
display(avg_ber_Algo3,'Average BER - Algorithm 3 in [R1]');
display(avg_ber_Algo3_low_complexity,'Average BER -low-complexity Algorithm 3 in [R1]');
display(avg_ber_MPA,'Average BER - Message passing algorithm (MPA)');
display(avg_ber_1tap,'Average BER - Single tap TF equalizer');
display(avg_ber_LMMSE,'Average BER - LMMSE equalizer');
display(avg_no_of_iterations_MRC,'Average number of iterations for the delay-time domain MRC detector');
display(avg_no_of_iterations_Algo3,'Average number of iterations for Algorithm 3 in [R1]');
display(avg_no_of_iterations_Algo3_low_complexity,'Average number of iterations for low-complexity Algorithm 3 in [R1]');
display(avg_no_of_iterations_MPA,'Average number of iterations for the MPA detector');