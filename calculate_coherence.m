
function [coherence_matrix, coherence_mean_matrix] = calculate_coherence(eeg,numChan,fs)
    
    eeg_coherences = cell(numChan,numChan); %we create a matrix where we we store a matrix in each position. That matrix will store the coh and f of all the possible combinations of channels
    coh_mean = zeros(numChan,numChan);
    
    for ch_first= 1:numChan
        for ch_second= 1:numChan
            try
                [coh, f] = mscohere(eeg(ch_first,:), eeg(ch_second,:), [], [], [], fs);
                eeg_coherences{ch_first, ch_second}= [f,coh];
                coh_mean(ch_first,ch_second) = mean(coh); %Matrix that stores the mean of the coerence between channels
            
            catch ME
                warning(strcat("Erro na coerência entre canal", num2str(ch_first), " e ",  num2str(ch_second), ME.message)); %Catch any errors that might occur and stop the loop
    
            end
        end
    end
    
    coherence_matrix= eeg_coherences;
    coherence_mean_matrix= coh_mean;


end