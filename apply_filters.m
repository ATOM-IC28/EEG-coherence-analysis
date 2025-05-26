%apply filters function
function filtered_eeg = apply_filters(eeg,varargin) %we use varargin so that the function can accept infinite numbers of filters. It is important to note that varargin is a cell
    numChan = size(eeg,1);
    eeg_filtered = zeros(size(eeg)); %first we create an empty matrix with the same size as the original matrix
    eeg_to_filter = eeg;

    for i = 1:length(varargin)
        for ch = 1:numChan
            Hd = varargin{i}; %Accessing the name of the filter
            eeg_filtered(ch, :) = filter(Hd, eeg_to_filter(ch, :)); %with the for loop we "fill" the matrix with the filtered values
        end

        eeg_to_filter = eeg_filtered; %Updating the EEG to filter, to make sure that the next filter will be applied on the already filtered EEG

    end

    filtered_eeg = eeg_filtered;

end

