%Extracting the data from the EEG files

%EEG - ACTIVE, EYES CLOSED
[data1,numChan,labels,txt,fs,gain,prefiltering,ChanDim] = eeg_read_bdf("CPFRO_calculo_OF_02.bdf","all","no");

%EEG - RESTING EYES CLOSED
[data2,~,~,~,~,~,~,~] = eeg_read_bdf("CPFRO_repouso_OF_01.bdf","all","no");

%EEG - RESTING, EYES OPEN   
[data3,~,~,~,~,~,~,~] = eeg_read_bdf("CPFRO_repouso_OA_01.bdf","all","no");

%%
%here we can see it is the 33rd channel that is the marker
%It appears empty
figure;
for ch = 1:33
    subplot(7,5,ch)
    plot(data1(ch,1:20*fs))
    title("Channel ",ch)

end

%We ignore the 33rd channel then, as it is a marker channel
numChan = numChan - 1;

%%
%XCreating a sort of dictionary to store the names of the EEG channels
%associated with the channel number
channel_labels = {'Fp1','AF3','F7','F3','FC1','FC5','T7','C3','CP1','CP5','P7','P3','Pz','PO3','O1','Oz','O2','PO4','P4','P8','CP6','CP2','C4','T8','FC6','FC2','F4','F8','AF4','Fp2','Fz','Cz'};
indexes = 1:length(channel_labels);
channels = containers.Map(channel_labels, indexes);

%Use channels("name of the channel") to access the channel index
%Use channels_label("number of channel") to access the channel label

%%
%Initialize eeglab
eeglab;

%%
%Visualize the electrodes
figure;
chanlocs = pop_readlocs('channel_locs.ced');
topoplot([], chanlocs, 'style', 'blank', 'electrodes', 'labels');
%%
%Calculating the duration of the signal
raw_eeg_1= data1(1:numChan,:);
raw_eeg_2= data2(1:numChan,:);
raw_eeg_3= data3(1:numChan,:);

signal_size = size(raw_eeg_1);
samples= signal_size(2);
duration= (samples/fs)/60; %duration of the signal in minutes

%%
%We remove the mean through time ("2")
%Removes the Offset of the data
raw_eeg_nooffset_01 = raw_eeg_1 - mean(raw_eeg_1,2);
raw_eeg_nooffset_02 = raw_eeg_2 - mean(raw_eeg_2,2); 
raw_eeg_nooffset_03 = raw_eeg_3 - mean(raw_eeg_3,2); 

%%
%Removing the trend in the data
%Using "detrend" we remove any linear trend in the EEG data

raw_eeg_nooffset_1 = detrend(raw_eeg_nooffset_01,"linear");
raw_eeg_nooffset_2 = detrend(raw_eeg_nooffset_02,"linear"); 
raw_eeg_nooffset_3 = detrend(raw_eeg_nooffset_03,"linear");

%%
%band_stop_50Hz filter (removing the frequency from the EU eletrical grid, 50Hz)
%Removing the noise from slow movements (< 1Hz)
%Removing the noise from high frequencies (> 100Hz)

eeg_filtered_1 = apply_filters(raw_eeg_nooffset_1,band_stop_50Hz,high_pass_filt,low_pass_filt);
eeg_filtered_2 = apply_filters(raw_eeg_nooffset_2,band_stop_50Hz,high_pass_filt,low_pass_filt);
eeg_filtered_3 = apply_filters(raw_eeg_nooffset_3,band_stop_50Hz,high_pass_filt,low_pass_filt);

%%
%Theta band (4-8Hz)
theta_band_1 = apply_filters(eeg_filtered_1, theta_band_filt);
theta_band_2 = apply_filters(eeg_filtered_2, theta_band_filt);
theta_band_3 = apply_filters(eeg_filtered_3, theta_band_filt);

%%
%Alpha band (8-13Hz)
alpha_band_1 = apply_filters(eeg_filtered_1, alpha_band_filt);
alpha_band_2 = apply_filters(eeg_filtered_2, alpha_band_filt);
alpha_band_3 = apply_filters(eeg_filtered_3, alpha_band_filt);

%%
%Beta band (13-30Hz)
beta_band_1 = apply_filters(eeg_filtered_1, beta_band_filt);
beta_band_2 = apply_filters(eeg_filtered_2, beta_band_filt);
beta_band_3 = apply_filters(eeg_filtered_3, beta_band_filt);

%%
%Gamma band (30-100Hz)
gamma_band_1 = apply_filters(eeg_filtered_1, gamma_band_filt);
gamma_band_2 = apply_filters(eeg_filtered_2, gamma_band_filt);
gamma_band_3 = apply_filters(eeg_filtered_3, gamma_band_filt);

%%
%Calculating the coherences and mean coherences matrixes (ALL FREQUENCIES)

[eeg_coherences_all_1, coh_all_mean_1] = calculate_coherence(eeg_filtered_1, numChan,fs);
[eeg_coherences_all_2, coh_all_mean_2] = calculate_coherence(eeg_filtered_2, numChan,fs);
[eeg_coherences_all_3, coh_all_mean_3] = calculate_coherence(eeg_filtered_3, numChan,fs);

%%
%Calculating the coherences and mean coherences matrixes (Theta band)
[theta_coh_1, theta_coh_mean_1] = calculate_coherence(theta_band_1,numChan,fs);
[theta_coh_2, theta_coh_mean_2] = calculate_coherence(theta_band_2,numChan,fs);
[theta_coh_3, theta_coh_mean_3] = calculate_coherence(theta_band_3,numChan,fs);

%%
%Calculating the coherences and mean coherences matrixes (Alpha band)
[alpha_coh_1, alpha_coh_mean_1] = calculate_coherence(alpha_band_1,numChan,fs);
[alpha_coh_2, alpha_coh_mean_2] = calculate_coherence(alpha_band_2,numChan,fs);
[alpha_coh_3, alpha_coh_mean_3] = calculate_coherence(alpha_band_3,numChan,fs);

%%
%Calculating the coherences and mean coherences matrixes (Beta band)
[beta_coh_1, beta_coh_mean_1] = calculate_coherence(beta_band_1,numChan,fs);
[beta_coh_2, beta_coh_mean_2] = calculate_coherence(beta_band_2,numChan,fs);
[beta_coh_3, beta_coh_mean_3] = calculate_coherence(beta_band_3,numChan,fs);

%%
%Calculating the coherences and mean coherences matrixes (Gamma band)
[gamma_coh_1, gamma_coh_mean_1] = calculate_coherence(gamma_band_1,numChan,fs);
[gamma_coh_2, gamma_coh_mean_2] = calculate_coherence(gamma_band_2,numChan,fs);
[gamma_coh_3, gamma_coh_mean_3] = calculate_coherence(gamma_band_3,numChan,fs);

%%
%Testing if the filter worked
ch=6;
filtered_signal=alpha_band_1;

subplot(2,1,1)
[Pxx,F] = pwelch(eeg_filtered_1(ch,:),[],[],[],fs);
plot(F,Pxx);
title(strcat("Before Filter - Channel ", num2str(ch)))
xlabel("Frequency (Hz)")
ylabel("PSD")
xlim([5,105])

subplot(2,1,2)
[Pxx,F] = pwelch(filtered_signal(ch,:),[],[],[],fs);
plot(F,Pxx);
title(strcat("After Filter - Channel ", num2str(ch)))
xlabel("Frequency (Hz)")
ylabel("PSD")
xlim([5,105])

%%
%Testing if the previous step worked
figure('Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8]);
ch1= 15;
ch2= 16;
inf_lim = 8;
sup_lim = 13;

coh_band_1 = alpha_coh_1;
coh_band_2 = alpha_coh_2;
coh_band_3 = alpha_coh_3;

coh_1 = coh_band_1{ch1,ch2}(:,2);
f_1 = coh_band_1{ch1,ch2}(:,1);
coh_smoothed_1 = movmean(coh_1, 15); %smoothed coeherences

subplot(3,1,1)
plot(f_1,coh_smoothed_1)
xlabel("Frequency (Hz)")
ylabel("Coherence (0-1)")
xlim([inf_lim,sup_lim])
title(strcat("Coherence between Channel ", channel_labels{ch1}, " and ", channel_labels{ch2}, " Active, Eyes closed"));

coh_2 = coh_band_2{ch1,ch2}(:,2);
f_2 = coh_band_2{ch1,ch2}(:,1);
coh_smoothed_2 = movmean(coh_2, 15); %smoothed coeherences

subplot(3,1,2)
plot(f_2,coh_smoothed_2)
xlabel("Frequency (Hz)")
ylabel("Coherence (0-1)")
xlim([inf_lim,sup_lim])
title(strcat("Coherence between Channel ", channel_labels{ch1}, " and ", channel_labels{ch2}, " Resting, Eyes closed"));

coh_3 = coh_band_3{ch1,ch2}(:,2);
f_3 = coh_band_3{ch1,ch2}(:,1);
coh_smoothed_3 = movmean(coh_3, 15); %smoothed coeherences

subplot(3,1,3)
plot(f_3,coh_smoothed_3)
xlabel("Frequency (Hz)")
ylabel("Coherence (0-1)")
xlim([inf_lim,sup_lim])
title(strcat("Coherence between Channel ", channel_labels{ch1}, " and ", channel_labels{ch2}, " Active, Eyes open"));

%%
%Analizing the mean of the coherence between all channels
figure('Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8]);
band = "Alpha band";
band_width = "(8-13Hz) (MEAN COHERENCE)";

coh_mean_1 = alpha_coh_mean_1;
coh_mean_2 = alpha_coh_mean_2;
coh_mean_3 = alpha_coh_mean_3;

coh_no_diag_1 = coh_mean_1;
coh_no_diag_1(tril(true(numChan), 0)) = 0;  % Removing the part of the matrix we do not need (the diagonal and one of the symetrical triangles)

coh_no_diag_2 = coh_mean_2;
coh_no_diag_2(tril(true(numChan), 0)) = 0;

coh_no_diag_3 = coh_mean_3;
coh_no_diag_3(tril(true(numChan), 0)) = 0;

% Transforming the matrix in a sorted vector (from biggest to lowest) and storing the indexes of the
% original values in the matrix
[~, linear_idx_1] = sort(coh_no_diag_1(:), "descend");
[~, linear_idx_2] = sort(coh_no_diag_2(:), "descend");
[~, linear_idx_3] = sort(coh_no_diag_3(:), "descend");

% Getting the top 10 values of coherence between channels
top_n = 10;

%top_vals_1 = sorted_vals_1(1:top_n);
top_idx_1 = linear_idx_1(1:top_n);

%top_vals_2 = sorted_vals_2(1:top_n);
top_idx_2 = linear_idx_2(1:top_n);

%top_vals_3 = sorted_vals_3(1:top_n);
top_idx_3 = linear_idx_3(1:top_n);

% Getting the position of the top coherence values in the original matrix
[row_idx_1, col_idx_1] = ind2sub(size(coh_mean_1), top_idx_1); %we use this function to return the position of the top coherence values in the original matrix
[row_idx_2, col_idx_2] = ind2sub(size(coh_mean_2), top_idx_2);
[row_idx_3, col_idx_3] = ind2sub(size(coh_mean_3), top_idx_3);


imagesc(coh_mean_1);
colorbar;
title(strcat(band, " Coherence ", band_width, " Active, Eyes closed"));
xticks(1:numChan);
yticks(1:numChan);
xticklabels(channel_labels);
yticklabels(channel_labels);
xlabel("Channels");
ylabel("Channels");
xtickangle(90);
hold on;

for i = 1:top_n
    plot(col_idx_1(i), row_idx_1(i), "ro", "Markersize", 12, "LineWidth", 1.5);
end
hold off;

figure('Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8]);
imagesc(coh_mean_2);
colorbar;
title(strcat(band, " Coherence ", band_width, " Resting, Eyes closed"));
xticks(1:numChan);
yticks(1:numChan);
xticklabels(channel_labels);
yticklabels(channel_labels);
xlabel("Channels");
ylabel("Channels");
xtickangle(90);
hold on;

for i = 1:top_n
    plot(col_idx_2(i), row_idx_2(i), "ro", "Markersize", 12, "LineWidth", 1.5);
end
hold off;

figure('Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8]);
imagesc(coh_mean_3);
colorbar;
title(strcat(band, " Coherence ", band_width, " Resting, Eyes open"));
xticks(1:numChan);
yticks(1:numChan);
xticklabels(channel_labels);
yticklabels(channel_labels);
xlabel("Channels");
ylabel("Channels");
xtickangle(90);
hold on;

for i = 1:top_n
    plot(col_idx_3(i), row_idx_3(i), "ro", "Markersize", 12, "LineWidth", 1.5);
end
hold off;

%%
%Analizing the difference in mean coherence between 
%the three different states on the alpha band
figure('Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8]);
band = "Alpha Band";
band_width = "(MEAN COHERENCE)";

coh_mean_1 = alpha_coh_mean_1;
coh_mean_2 = alpha_coh_mean_2;
coh_mean_3 = alpha_coh_mean_3;

diff_1_2 = coh_mean_1 - coh_mean_2;
diff_1_3 = coh_mean_1 - coh_mean_3;
diff_2_3 = coh_mean_2 - coh_mean_3;

% Making all the coherences in the diagonal 0 so it doesnt affect the
% calculus of the top coherences between channels
coh_no_diag_1 = diff_1_2;
coh_no_diag_1(tril(true(numChan), 0)) = 0;  % Removing the part of the matrix we do not need (the diagonal and one of the symetrical triangles)

coh_no_diag_2 = diff_1_3;
coh_no_diag_2(tril(true(numChan), 0)) = 0;

coh_no_diag_3 = diff_2_3;
coh_no_diag_3(tril(true(numChan), 0)) = 0;

% Transforming the matrix in a sorted vector (from biggest to lowest) and storing the indexes of the
% original values in the matrix
[~, linear_idx_1] = sort(coh_no_diag_1(:), "descend");
[~, linear_idx_2] = sort(coh_no_diag_2(:), "descend");
[~, linear_idx_3] = sort(coh_no_diag_3(:), "descend");

% Getting the top 10 values of coherence between channels
top_n = 10;

%top_vals_1 = sorted_vals_1(1:top_n);
top_idx_1 = linear_idx_1(1:top_n);

%top_vals_2 = sorted_vals_2(1:top_n);
top_idx_2 = linear_idx_2(1:top_n);

%top_vals_3 = sorted_vals_3(1:top_n);
top_idx_3 = linear_idx_3(1:top_n);

% Getting the position of the top coherence values in the original matrix
[row_idx_1, col_idx_1] = ind2sub(size(coh_mean_1), top_idx_1); %we use this function to return the position of the top coherence values in the original matrix
[row_idx_2, col_idx_2] = ind2sub(size(coh_mean_2), top_idx_2);
[row_idx_3, col_idx_3] = ind2sub(size(coh_mean_3), top_idx_3);

max_diff = max(abs([diff_1_2(:); diff_1_3(:); diff_2_3(:)]));

imagesc(diff_1_2);
clim([-max_diff, max_diff]); %So we can more easily compare between the 3 states
colorbar;
title(strcat(band, " Coherence ", band_width), (" Active, Eyes closed VS Resting, Eyes Closed"));
xticks(1:numChan);
yticks(1:numChan);
xticklabels(channel_labels);
yticklabels(channel_labels);
xlabel("Channels");
ylabel("Channels");
xtickangle(90);

hold on;

for i = 1:top_n
    plot(col_idx_1(i), row_idx_1(i), "ro", "Markersize", 12, "LineWidth", 1.5); 
end
hold off;

figure('Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8]);
imagesc(diff_1_3);
clim([-max_diff, max_diff]);
colorbar;
title(strcat(band, " Coherence ", band_width),(" Active, Eyes closed VS Resting, Eyes open"));
xticks(1:numChan);
yticks(1:numChan);
xticklabels(channel_labels);
yticklabels(channel_labels);
xlabel("Channels");
ylabel("Channels");
xtickangle(90);
hold on;

for i = 1:top_n
    plot(col_idx_2(i), row_idx_2(i), "ro", "Markersize", 12, "LineWidth", 1.5);
end
hold off;

figure('Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8]);
imagesc(diff_2_3);
clim([-max_diff, max_diff]);
colorbar;
title(strcat(band, " Coherence ", band_width),( " Resting, Eyes closed VS Resting, Eyes open "));
xticks(1:numChan);
yticks(1:numChan);
xticklabels(channel_labels);
yticklabels(channel_labels);
xlabel("Channels");
ylabel("Channels");
xtickangle(90);

hold on;

for i = 1:top_n
    plot(col_idx_3(i), row_idx_3(i), "ro", "Markersize", 12, "LineWidth", 1.5); 
end
hold off;
 

%%
%Analyzing the channels with the most coherence between them (base code)

coh_mean_1 = alpha_coh_mean_1;
coh_mean_2 = alpha_coh_mean_2;
coh_mean_3 = alpha_coh_mean_3;

% Making all the coherences in the diagonal 0 so it doesnt affect the
% calculus of the top coherences between channels
coh_no_diag_1 = coh_mean_1;
coh_no_diag_1(tril(true(numChan), 0)) = 0;  % Removing the part of the matrix we do not need (the diagonal and one of the symetrical triangles)

coh_no_diag_2 = coh_mean_2;
coh_no_diag_2(tril(true(numChan), 0)) = 0;

coh_no_diag_3 = coh_mean_3;
coh_no_diag_3(tril(true(numChan), 0)) = 0;

% Transforming the matrix in a sorted vector (from biggest to lowest) and storing the indexes of the
% original values in the matrix
[sorted_vals_1, linear_idx_1] = sort(coh_no_diag_1(:), "descend");
[sorted_vals_2, linear_idx_2] = sort(coh_no_diag_2(:), "descend");
[sorted_vals_3, linear_idx_3] = sort(coh_no_diag_3(:), "descend");

% Getting the top 10 values of coherence between channels
top_n = 10;

top_vals_1 = sorted_vals_1(1:top_n);
top_idx_1 = linear_idx_1(1:top_n);

top_vals_2 = sorted_vals_2(1:top_n);
top_idx_2 = linear_idx_2(1:top_n);

top_vals_3 = sorted_vals_3(1:top_n);
top_idx_3 = linear_idx_3(1:top_n);

% Getting the position of the top coherence values in the original matrix
[row_idx_1, col_idx_1] = ind2sub(size(coh_mean_1), top_idx_1); %we use this function to return the position of the top coherence values in the original matrix
[row_idx_2, col_idx_2] = ind2sub(size(coh_mean_2), top_idx_2);
[row_idx_3, col_idx_3] = ind2sub(size(coh_mean_3), top_idx_3);

%Active, Eyes closed
figure;
imagesc(coh_mean_1);
colorbar;
title(strcat(band, " Coherence ", band_width, " Active, Eyes closed"),("Top 10 values highlighted"));
xticks(1:numChan);
yticks(1:numChan);
xticklabels(channel_labels);
yticklabels(channel_labels);
xtickangle(90);
hold on;

for i = 1:top_n
    plot(col_idx_1(i), row_idx_1(i), "ro", "Markersize", 12, "LineWidth", 1.5);
    plot(row_idx_1(i), col_idx_1(i), "ro", "Markersize", 12, "LineWidth", 1.5); 
end
hold off;

%Resting, eyes closed
figure;
imagesc(coh_mean_2);
colorbar;
title(strcat(band, " Coherence ", band_width, " Resting, Eyes closed"),("Top 10 values highlighted"));
xticks(1:numChan);
yticks(1:numChan);
xticklabels(channel_labels);
yticklabels(channel_labels);
xtickangle(90);
hold on;

for i = 1:top_n
    plot(col_idx_2(i), row_idx_2(i), "ro", "Markersize", 12, "LineWidth", 1.5);
    plot(row_idx_2(i), col_idx_2(i), "ro", "Markersize", 12, "LineWidth", 1.5); 
end
hold off;

%Resting, eyes open
figure;
imagesc(coh_mean_3);
colorbar;
title(strcat(band, " Coherence ", band_width, " Resting, Eyes open"),("Top 10 values highlighted"));
xticks(1:numChan);
yticks(1:numChan);
xticklabels(channel_labels);
yticklabels(channel_labels);
xtickangle(90);
hold on;

for i = 1:top_n
    plot(col_idx_3(i), row_idx_3(i), "ro", "Markersize", 12, "LineWidth", 1.5);
    plot(row_idx_3(i), col_idx_3(i), "ro", "Markersize", 12, "LineWidth", 1.5);
end
hold off;

%%
%Ploting the mean coherence between a fixed channel and all the other
%channels in the three states in a specific frequency band

coh_mean_1 = alpha_coh_mean_1;
coh_mean_2 = alpha_coh_mean_2;
coh_mean_3 = alpha_coh_mean_3;

band = "Alpha band";

fixed_chan = channels('Oz');
chanlocs = pop_readlocs('channel_locs.ced');

coh_vec_1 = coh_mean_1(fixed_chan,:);
coh_vec_2 = coh_mean_2(fixed_chan,:);
coh_vec_3 = coh_mean_3(fixed_chan,:);


all_vals = [coh_vec_1(:); coh_vec_2(:); coh_vec_3(:)]; %Preparing the values as a 1D vector so we can use the min and max functions
clim = [min(all_vals, [], "omitnan"), max(all_vals, [], "omitnan")];

figure;
%Active, Eyes Closed
subplot(3,1,1)
topoplot(coh_vec_1, chanlocs, 'style', 'map', 'electrodes', 'on', 'maplimits',clim);
title(strcat("Active, Eyes Closed - ", band, " - Coeherences relative to channel ", channel_labels{fixed_chan}))
colorbar;

%Resting, Eyes Closed
subplot(3,1,2)
topoplot(coh_vec_2, chanlocs, 'style', 'map', 'electrodes', 'on','maplimits',clim);
title(strcat("Resting, Eyes Closed - ", band, " - Coeherences relative to channel ", channel_labels{fixed_chan}))
colorbar;

%Resting, Eyes Open
subplot(3,1,3)
topoplot(coh_vec_3, chanlocs, 'style', 'map', 'electrodes', 'on', 'maplimits',clim);
title(strcat("Resting, Eyes Open - ", band, " - Coeherences relative to channel ", channel_labels{fixed_chan}))
colorbar;

%%
%Ploting the mean coherence between a fixed channel and all the other
%channels for all 32 channels in one specific state and frequency band

coh_mean_1 = theta_coh_mean_1;
coh_mean_2 = theta_coh_mean_2;
coh_mean_3 = theta_coh_mean_3;

band = "Theta band";

fixed_chan = channels('Fp1');
chanlocs = pop_readlocs('channel_locs.ced');

coh_vec_1 = coh_mean_1(fixed_chan,:);
coh_vec_2 = coh_mean_2(fixed_chan,:);
coh_vec_3 = coh_mean_3(fixed_chan,:);


all_vals = [coh_vec_1(:); coh_vec_2(:); coh_vec_3(:)]; %Preparing the values as a 1D vector so we can use the min and max functions
clim = [min(all_vals, [], "omitnan"), max(all_vals, [], "omitnan")];

%Active, Eyes Closed

figure('Name', strcat('Active, Eyes Closed - ', band), 'NumberTitle', 'off');

for ch= 1:numChan
    subplot(8,4,ch)
    coh_vec_1 = coh_mean_1(ch,:);
    coh_vec_1(ch) = 0;
    topoplot(coh_vec_1, chanlocs, ...
    'style', 'map', ...
    'electrodes', 'on', ...
    'maplimits',clim, ...
    'emarker2', {ch, 'o', 'w', 3, 1}); 
    title(strcat("Fixed Channel ", channel_labels{ch}))
    colorbar;
end

%Resting, Eyes Closed

figure('Name', strcat('Resting, Eyes Closed - ', band), 'NumberTitle', 'off');

for ch= 1:numChan
    subplot(8,4,ch)
    coh_vec_2 = coh_mean_2(ch,:);
    coh_vec_2(ch) = 0;
    topoplot(coh_vec_2, chanlocs, ...
    'style', 'map', ...
    'electrodes', 'on', ...
    'maplimits',clim, ...
    'emarker2', {ch, 'o', 'w', 3, 1});  
    title(strcat("Fixed Channel ", channel_labels{ch}))
    colorbar;

end

%Resting, Eyes Open

figure('Name', strcat('Resting, Eyes Open - ', band), 'NumberTitle', 'off');

for ch= 1:numChan
    subplot(8,4,ch)
    coh_vec_3 = coh_mean_3(ch,:);
    coh_vec_3(ch) = 0;
    topoplot(coh_vec_3, chanlocs, ...
    'style', 'map', ...
    'electrodes', 'on', ...
    'maplimits',clim, ...
    'emarker2', {ch, 'o', 'w', 3, 1});  
    title(strcat("Fixed Channel ", channel_labels{ch}))
    colorbar;

end

