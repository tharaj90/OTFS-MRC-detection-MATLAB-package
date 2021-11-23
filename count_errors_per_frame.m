     %% errors count%%%%%
        errors_MRC = sum(xor(est_info_bits_MRC,trans_info_bit));
        errors_Algo1 = sum(xor(est_info_bits_Algo1,trans_info_bit));
        errors_Algo3 = sum(xor(est_info_bits_Algo3,trans_info_bit));
        errors_Algo3_low_complexity = sum(xor(est_info_bits_Algo3_low_complexity,trans_info_bit));
        errors_MPA = sum(xor(est_info_bits_MPA,trans_info_bit));
        errors_1tap = sum(xor(est_info_bits_1tap,trans_info_bit));
        errors_LMMSE = sum(xor(est_info_bits_LMMSE,trans_info_bit));
        
        err_ber_MRC(1,iesn0) = err_ber_MRC(1,iesn0) + errors_MRC;
        err_ber_Algo1(1,iesn0) = err_ber_Algo1(1,iesn0) + errors_Algo1;
        err_ber_Algo3(1,iesn0) = err_ber_Algo3(1,iesn0) + errors_Algo3;
        err_ber_Algo3_low_complexity(1,iesn0) = err_ber_Algo3_low_complexity(1,iesn0) + errors_Algo3_low_complexity;
        err_ber_MPA(1,iesn0) = err_ber_MPA(1,iesn0) + errors_MPA;
        err_ber_1tap(1,iesn0) = err_ber_1tap(1,iesn0) + errors_1tap;
        err_ber_LMMSE(1,iesn0) = err_ber_LMMSE(1,iesn0) + errors_LMMSE;
        
        avg_ber_MRC(1,iesn0)=err_ber_MRC(1,iesn0).'/length(trans_info_bit)/ifram;
        avg_ber_Algo1(1,iesn0)=err_ber_Algo1(1,iesn0).'/length(trans_info_bit)/ifram;
        avg_ber_Algo3(1,iesn0)=err_ber_Algo3(1,iesn0).'/length(trans_info_bit)/ifram;
        avg_ber_Algo3_low_complexity(1,iesn0)=err_ber_Algo3_low_complexity(1,iesn0).'/length(trans_info_bit)/ifram;
        avg_ber_MPA(1,iesn0)=err_ber_MPA(1,iesn0).'/length(trans_info_bit)/ifram;
        avg_ber_1tap(1,iesn0)=err_ber_1tap(1,iesn0).'/length(trans_info_bit)/ifram;
        avg_ber_LMMSE(1,iesn0)=err_ber_LMMSE(1,iesn0).'/length(trans_info_bit)/ifram;
        
        
        %%  iterations count
        
        no_of_detetor_iterations_MRC(iesn0)=no_of_detetor_iterations_MRC(iesn0)+det_iters_MRC;
        avg_no_of_iterations_MRC(iesn0)=no_of_detetor_iterations_MRC(iesn0)/ifram;
        no_of_detetor_iterations_Algo1(iesn0)=no_of_detetor_iterations_Algo1(iesn0)+det_iters_Algo1;
        avg_no_of_iterations_Algo1(iesn0)=no_of_detetor_iterations_Algo1(iesn0)/ifram;
        no_of_detetor_iterations_Algo3(iesn0)=no_of_detetor_iterations_Algo3(iesn0)+det_iters_Algo3;
        avg_no_of_iterations_Algo3(iesn0)=no_of_detetor_iterations_Algo3(iesn0)/ifram;
        no_of_detetor_iterations_Algo3_low_complexity(iesn0)=no_of_detetor_iterations_Algo3_low_complexity(iesn0)+det_iters_Algo3_low_complexity;
        avg_no_of_iterations_Algo3_low_complexity(iesn0)=no_of_detetor_iterations_Algo3_low_complexity(iesn0)/ifram;
        no_of_detetor_iterations_MPA(iesn0)=no_of_detetor_iterations_MPA(iesn0)+det_iters_MPA;
        avg_no_of_iterations_MPA(iesn0)=no_of_detetor_iterations_MPA(iesn0)/ifram;