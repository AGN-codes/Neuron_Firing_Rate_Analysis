%% mgs
clc;
clear all;
close all;

%% loading the data
load("data.mat");

%% neuron selection
neuron_no = 15;
SoM = 1; % 1 = multi, 2 = single

%% Q2 - a
trialno_led51 = (Event.mgs.codes(:,5)==51)';
trialno_led52 = (Event.mgs.codes(:,5)==52)';
trialno_led53 = (Event.mgs.codes(:,5)==53)';
trialno_led54 = (Event.mgs.codes(:,5)==54)';
trialno_led55 = (Event.mgs.codes(:,5)==55)';
trialno_led56 = (Event.mgs.codes(:,5)==56)';
trialno_led57 = (Event.mgs.codes(:,5)==57)';
trialno_led58 = (Event.mgs.codes(:,5)==58)';

% ===> select multi or single here
% 1 for multi, 2 for single
trials_led51 = spike{neuron_no,SoM}.mgs(trialno_led51, :);
trials_led52 = spike{neuron_no,SoM}.mgs(trialno_led52, :);
trials_led53 = spike{neuron_no,SoM}.mgs(trialno_led53, :);
trials_led54 = spike{neuron_no,SoM}.mgs(trialno_led54, :);
trials_led55 = spike{neuron_no,SoM}.mgs(trialno_led55, :);
trials_led56 = spike{neuron_no,SoM}.mgs(trialno_led56, :);
trials_led57 = spike{neuron_no,SoM}.mgs(trialno_led57, :);
trials_led58 = spike{neuron_no,SoM}.mgs(trialno_led58, :);

%% %% Q2 - a - raster process
raster_trials_led51 = imresize(trials_led51,[2001,4001], 'nearest');
raster_trials_led52 = imresize(trials_led52,[2001,4001], 'nearest');
raster_trials_led53 = imresize(trials_led53,[2001,4001], 'nearest');
raster_trials_led54 = imresize(trials_led54,[2001,4001], 'nearest');
raster_trials_led55 = imresize(trials_led55,[2001,4001], 'nearest');
raster_trials_led56 = imresize(trials_led56,[2001,4001], 'nearest');
raster_trials_led57 = imresize(trials_led57,[2001,4001], 'nearest');
raster_trials_led58 = imresize(trials_led58,[2001,4001], 'nearest');

raster_trials_led51 = imresize(bwmorph(raster_trials_led51, 'thicken', 4),0.5);
raster_trials_led52 = imresize(bwmorph(raster_trials_led52, 'thicken', 4),0.5);
raster_trials_led53 = imresize(bwmorph(raster_trials_led53, 'thicken', 4),0.5);
raster_trials_led54 = imresize(bwmorph(raster_trials_led54, 'thicken', 4),0.5);
raster_trials_led55 = imresize(bwmorph(raster_trials_led55, 'thicken', 4),0.5);
raster_trials_led56 = imresize(bwmorph(raster_trials_led56, 'thicken', 4),0.5);
raster_trials_led57 = imresize(bwmorph(raster_trials_led57, 'thicken', 4),0.5);
raster_trials_led58 = imresize(bwmorph(raster_trials_led58, 'thicken', 4),0.5);

raster_trials_led51 = 1 - raster_trials_led51;
raster_trials_led52 = 1 - raster_trials_led52;
raster_trials_led53 = 1 - raster_trials_led53;
raster_trials_led54 = 1 - raster_trials_led54;
raster_trials_led55 = 1 - raster_trials_led55;
raster_trials_led56 = 1 - raster_trials_led56;
raster_trials_led57 = 1 - raster_trials_led57;
raster_trials_led58 = 1 - raster_trials_led58;

%% Q2 - a - raster plot
figure;
tiledlayout(4,2);

nexttile;
imshow(raster_trials_led51);
ylabel('led 51');
nexttile;
imshow(raster_trials_led52);
ylabel('led 52');
nexttile;
imshow(raster_trials_led53);
ylabel('led 53');
nexttile;
imshow(raster_trials_led54);
ylabel('led 54');
nexttile;
imshow(raster_trials_led55);
ylabel('led 55');
nexttile;
imshow(raster_trials_led56);
ylabel('led 56');
nexttile;
imshow(raster_trials_led57);
ylabel('led 57');
nexttile;
imshow(raster_trials_led58);
ylabel('led 58');

%% Q2 - a - psth proccess
psth_trials_led51 = sum(trials_led51, 1)/size(trials_led51,1);
psth_trials_led52 = sum(trials_led52, 1)/size(trials_led52,1);
psth_trials_led53 = sum(trials_led53, 1)/size(trials_led53,1);
psth_trials_led54 = sum(trials_led54, 1)/size(trials_led54,1);
psth_trials_led55 = sum(trials_led55, 1)/size(trials_led55,1);
psth_trials_led56 = sum(trials_led56, 1)/size(trials_led56,1);
psth_trials_led57 = sum(trials_led57, 1)/size(trials_led57,1);
psth_trials_led58 = sum(trials_led58, 1)/size(trials_led58,1);

all_psth = [psth_trials_led51;psth_trials_led52;psth_trials_led53; ...
    psth_trials_led54;psth_trials_led55;psth_trials_led56; ...
    psth_trials_led57; psth_trials_led58];
max_all_psth = max(all_psth, [], 'all');

%% Q2 - a - psth plot
figure;
tiledlayout(4,2);

nexttile;
bar(psth_trials_led51,2);
ylabel('led 51');
xlim([1,4001]);
ylim([0, max_all_psth]);
nexttile;
bar(psth_trials_led52,2);
ylabel('led 52');
xlim([1,4001]);
ylim([0, max_all_psth]);
nexttile;
bar(psth_trials_led53,2);
ylabel('led 53');
xlim([1,4001]);
ylim([0, max_all_psth]);
nexttile;
bar(psth_trials_led54,2);
ylabel('led 54');
xlim([1,4001]);
ylim([0, max_all_psth]);
nexttile;
bar(psth_trials_led55,2);
ylabel('led 55');
xlim([1,4001]);
ylim([0, max_all_psth]);
nexttile;
bar(psth_trials_led56,2);
ylabel('led 56');
xlim([1,4001]);
ylim([0, max_all_psth]);
nexttile;
bar(psth_trials_led57,2);
ylabel('led 7');
xlim([1,4001]);
ylim([0, max_all_psth]);
nexttile;
bar(psth_trials_led58,2);
ylabel('led 58');
xlim([1,4001]);
ylim([0, max_all_psth]);

%% Q2 - a - mean firing rate
'mean firing rate selected neuron'
mean_trials_led51 = sum(psth_trials_led51)/4001
mean_trials_led52 = sum(psth_trials_led52)/4001
mean_trials_led53 = sum(psth_trials_led53)/4001
mean_trials_led54 = sum(psth_trials_led54)/4001
mean_trials_led55 = sum(psth_trials_led55)/4001
mean_trials_led56 = sum(psth_trials_led56)/4001
mean_trials_led57 = sum(psth_trials_led57)/4001
mean_trials_led58 = sum(psth_trials_led58)/4001

% ===> IN for 15'th neuron is led no. 51
% ===> OUT for 15'th neuron is led no. 55

%% Q2 - a - IN&OUT psth plot
figure;
hold on;

bar(psth_trials_led51,2);

xlim([1,4001]);
ylim([0, max_all_psth]);

bar(psth_trials_led55,2);

xlim([1,4001]);
ylim([0, max_all_psth]);

ylabel('psth for in and out');
legend('IN', 'OUT');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Q2 - b - gaussian kernel parameters
var_gaus = 10000;
fil_length = 1+250;
fil_gaus = exp( - ( (1:fil_length) - round(fil_length/2) ).^2 / var_gaus) / sqrt(2*pi*var_gaus);

% figure;
% plot(fil_gaus);
% title('gaussian kernel');
% 
% psth_trials_signle = zeros(1,4001);
% psth_trials_signle(2000) = 1;
% figure;
% plot(conv(psth_trials_signle, fil_gaus, 'same'));


%% Q2 - b - gaussian kernel
gaus_psth_trials_led51 = conv(psth_trials_led51, fil_gaus, 'same');
gaus_psth_trials_led52 = conv(psth_trials_led52, fil_gaus, 'same');
gaus_psth_trials_led53 = conv(psth_trials_led53, fil_gaus, 'same');
gaus_psth_trials_led54 = conv(psth_trials_led54, fil_gaus, 'same');
gaus_psth_trials_led55 = conv(psth_trials_led55, fil_gaus, 'same');
gaus_psth_trials_led56 = conv(psth_trials_led56, fil_gaus, 'same');
gaus_psth_trials_led57 = conv(psth_trials_led57, fil_gaus, 'same');
gaus_psth_trials_led58 = conv(psth_trials_led58, fil_gaus, 'same');

max_all_gaus = [gaus_psth_trials_led51;gaus_psth_trials_led52;gaus_psth_trials_led53; ...
    gaus_psth_trials_led54;gaus_psth_trials_led55;gaus_psth_trials_led56; ...
    gaus_psth_trials_led57; gaus_psth_trials_led58];
max_all_gaus = max(max_all_gaus, [], 'all');

%% Q2 - b - gaussian kernel plot
figure;
tiledlayout(4,2);
sgtitle('Gaussian Kernel Result');

nexttile;
plot(gaus_psth_trials_led51);
ylabel('led 51');
xlim([1,4001]);
ylim([0, max_all_gaus]);
nexttile;
plot(gaus_psth_trials_led52);
ylabel('led 52');
xlim([1,4001]);
ylim([0, max_all_gaus]);
nexttile;
plot(gaus_psth_trials_led53);
ylabel('led 53');
xlim([1,4001]);
ylim([0, max_all_gaus]);
nexttile;
plot(gaus_psth_trials_led54);
ylabel('led 54');
xlim([1,4001]);
ylim([0, max_all_gaus]);
nexttile;
plot(gaus_psth_trials_led55);
ylabel('led 55');
xlim([1,4001]);
ylim([0, max_all_gaus]);
nexttile;
plot(gaus_psth_trials_led56);
ylabel('led 56');
xlim([1,4001]);
ylim([0, max_all_gaus]);
nexttile;
plot(gaus_psth_trials_led57);
ylabel('led 57');
xlim([1,4001]);
ylim([0, max_all_gaus]);
nexttile;
plot(gaus_psth_trials_led58);
ylabel('led 58');
xlim([1,4001]);
ylim([0, max_all_gaus]);

figure;
sgtitle('Gaussian Kernel Result');
hold on;
plot(gaus_psth_trials_led51, 'LineWidth', 3);
plot(gaus_psth_trials_led52, 'LineWidth', 1);
plot(gaus_psth_trials_led53, 'LineWidth', 1);
plot(gaus_psth_trials_led54, 'LineWidth', 1);
plot(gaus_psth_trials_led55, 'LineWidth', 3);
plot(gaus_psth_trials_led56, 'LineWidth', 1);
plot(gaus_psth_trials_led57, 'LineWidth', 1);
plot(gaus_psth_trials_led58, 'LineWidth', 1);
xlim([1,4001]);
ylim([0, max_all_gaus]);
hold off;
legend('LED51','LED52','LED53','LED54','LED55','LED56','LED57','LED58');


%% Q2 - b - window kernel paramters
fil_length = 1+300;
fil_wind = ones(1,fil_length);

% figure;
% plot(fil_wind);
% title('window kernel');
% 
% psth_trials_signle = zeros(1,4001);
% psth_trials_signle(2000) = 1;
% figure;
% plot(conv(psth_trials_signle, fil_wind, 'same'));


%% Q2 - b - winow kernel proccess
wind_psth_trials_led51 = conv(psth_trials_led51, fil_wind, 'same');
wind_psth_trials_led52 = conv(psth_trials_led52, fil_wind, 'same');
wind_psth_trials_led53 = conv(psth_trials_led53, fil_wind, 'same');
wind_psth_trials_led54 = conv(psth_trials_led54, fil_wind, 'same');
wind_psth_trials_led55 = conv(psth_trials_led55, fil_wind, 'same');
wind_psth_trials_led56 = conv(psth_trials_led56, fil_wind, 'same');
wind_psth_trials_led57 = conv(psth_trials_led57, fil_wind, 'same');
wind_psth_trials_led58 = conv(psth_trials_led58, fil_wind, 'same');

max_all_wind = [wind_psth_trials_led51;wind_psth_trials_led52;wind_psth_trials_led53; ...
    wind_psth_trials_led54;wind_psth_trials_led55;wind_psth_trials_led56; ...
    wind_psth_trials_led57; wind_psth_trials_led58];
max_all_wind = max(max_all_wind, [], 'all');

%% Q2 - b - window kernel plot
figure;
tiledlayout(4,2);

nexttile;
plot(wind_psth_trials_led51);
ylabel('led 51');
xlim([1,4001]);
ylim([0, max_all_wind]);
nexttile;
plot(wind_psth_trials_led52);
ylabel('led 52');
xlim([1,4001]);
ylim([0, max_all_wind]);
nexttile;
plot(wind_psth_trials_led53);
ylabel('led 53');
xlim([1,4001]);
ylim([0, max_all_wind]);
nexttile;
plot(wind_psth_trials_led54);
ylabel('led 54');
xlim([1,4001]);
ylim([0, max_all_wind]);
nexttile;
plot(wind_psth_trials_led55);
ylabel('led 55');
xlim([1,4001]);
ylim([0, max_all_wind]);
nexttile;
plot(wind_psth_trials_led56);
ylabel('led 56');
xlim([1,4001]);
ylim([0, max_all_wind]);
nexttile;
plot(wind_psth_trials_led57);
ylabel('led 57');
xlim([1,4001]);
ylim([0, max_all_wind]);
nexttile;
plot(wind_psth_trials_led58);
ylabel('led 58');
xlim([1,4001]);
ylim([0, max_all_wind]);

figure;
hold on;
plot(wind_psth_trials_led51, 'LineWidth', 3);
plot(wind_psth_trials_led52, 'LineWidth', 1);
plot(wind_psth_trials_led53, 'LineWidth', 1);
plot(wind_psth_trials_led54, 'LineWidth', 1);
plot(wind_psth_trials_led55, 'LineWidth', 3);
plot(wind_psth_trials_led56, 'LineWidth', 1);
plot(wind_psth_trials_led57, 'LineWidth', 1);
plot(wind_psth_trials_led58, 'LineWidth', 1);
hold off;
legend('LED51','LED52','LED53','LED54','LED55','LED56','LED57','LED58');

%% Q2 - b - window kernel 2 parameter and proccess
fil_length = 1+49;

win2_psth_trials_led51 = window(psth_trials_led51, fil_length);
win2_psth_trials_led52 = window(psth_trials_led52, fil_length);
win2_psth_trials_led53 = window(psth_trials_led53, fil_length);
win2_psth_trials_led54 = window(psth_trials_led54, fil_length);
win2_psth_trials_led55 = window(psth_trials_led55, fil_length);
win2_psth_trials_led56 = window(psth_trials_led56, fil_length);
win2_psth_trials_led57 = window(psth_trials_led57, fil_length);
win2_psth_trials_led58 = window(psth_trials_led58, fil_length);

%% Q2 - b - window kernel 2 plot
figure;
hold on;
plot(win2_psth_trials_led51, 'LineWidth', 3);
plot(win2_psth_trials_led52, 'LineWidth', 1);
plot(win2_psth_trials_led53, 'LineWidth', 1);
plot(win2_psth_trials_led54, 'LineWidth', 1);
plot(win2_psth_trials_led55, 'LineWidth', 3);
plot(win2_psth_trials_led56, 'LineWidth', 1);
plot(win2_psth_trials_led57, 'LineWidth', 1);
plot(win2_psth_trials_led58, 'LineWidth', 1);
hold off;
legend('LED51','LED52','LED53','LED54','LED55','LED56','LED57','LED58');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Q2 - c - selecting in and out for fef
% ===> select multi or single here
% 1 for multi, 2 for single
feftrials_led51 = spike{17,SoM}.mgs(trialno_led51, :);
feftrials_led52 = spike{17,SoM}.mgs(trialno_led52, :);
feftrials_led53 = spike{17,SoM}.mgs(trialno_led53, :);
feftrials_led54 = spike{17,SoM}.mgs(trialno_led54, :);
feftrials_led55 = spike{17,SoM}.mgs(trialno_led55, :);
feftrials_led56 = spike{17,SoM}.mgs(trialno_led56, :);
feftrials_led57 = spike{17,SoM}.mgs(trialno_led57, :);
feftrials_led58 = spike{17,SoM}.mgs(trialno_led58, :);
fefpsth_trials_led51 = sum(feftrials_led51, 1)/size(feftrials_led51,1);
fefpsth_trials_led52 = sum(feftrials_led52, 1)/size(feftrials_led51,1);
fefpsth_trials_led53 = sum(feftrials_led53, 1)/size(feftrials_led51,1);
fefpsth_trials_led54 = sum(feftrials_led54, 1)/size(feftrials_led51,1);
fefpsth_trials_led55 = sum(feftrials_led55, 1)/size(feftrials_led51,1);
fefpsth_trials_led56 = sum(feftrials_led56, 1)/size(feftrials_led51,1);
fefpsth_trials_led57 = sum(feftrials_led57, 1)/size(feftrials_led51,1);
fefpsth_trials_led58 = sum(feftrials_led58, 1)/size(feftrials_led51,1);
'mean firing rate selected neuron'
fefmean_trials_led51 = sum(fefpsth_trials_led51)/4001
fefmean_trials_led52 = sum(fefpsth_trials_led52)/4001
fefmean_trials_led53 = sum(fefpsth_trials_led53)/4001
fefmean_trials_led54 = sum(fefpsth_trials_led54)/4001
fefmean_trials_led55 = sum(fefpsth_trials_led55)/4001
fefmean_trials_led56 = sum(fefpsth_trials_led56)/4001
fefmean_trials_led57 = sum(fefpsth_trials_led57)/4001
fefmean_trials_led58 = sum(fefpsth_trials_led58)/4001

% ===> IN for 15'th neuron is led no. 53
% ===> OUT for 15'th neuron is led no. 55

%% Q2 - c - average firing rate - multi
fr_in_fixation = sum(psth_trials_led51(300:900))/(601)
fr_out_fixation = sum(psth_trials_led55(300:900))/(601)

fr_in_visual = sum(psth_trials_led51(1000:1600))/(1600-900+1)
fr_out_visual = sum(psth_trials_led55(1000:1600))/(1600-900+1)

fr_in_memory = sum(psth_trials_led51(2400:3000))/(3000-2400+1)
fr_out_memory = sum(psth_trials_led55(2400:3000))/(3000-2400+1)

fr_in_all = sum(psth_trials_led51)/4001
fr_out_all = sum(psth_trials_led55)/4001



fef_fr_in_fixation = sum(fefpsth_trials_led51(300:900))/(601)
fef_fr_out_fixation = sum(fefpsth_trials_led55(300:900))/(601)

fef_fr_in_visual = sum(fefpsth_trials_led51(1000:1600))/(1600-900+1)
fef_fr_out_visual = sum(fefpsth_trials_led55(1000:1600))/(1600-900+1)

fef_fr_in_memory = sum(fefpsth_trials_led51(2400:3000))/(3000-2400+1)
fef_fr_out_memory = sum(fefpsth_trials_led55(2400:3000))/(3000-2400+1)

fef_fr_in_all = sum(fefpsth_trials_led51)/4001
fef_fr_out_all = sum(fefpsth_trials_led55)/4001

%% Q2 - c - plotting average firing rate
figure;
tiledlayout(2,1);

nexttile;
bary = [fr_in_fixation,fr_out_fixation,fr_in_visual,fr_out_visual,fr_in_memory,fr_out_memory,fr_in_all,fr_out_all];
hold on;
for i = 1:length(bary)
    bar(i,bary(i));
end
hold off;
grid on;
ylabel('firing rate');
title('firing rate for selected neuron');
legend('fixation in','fixation out','visual in','visual out','memory in','memory out','all in', 'all out', 'Location','eastoutside');

nexttile;
baryf = [fef_fr_in_fixation,fef_fr_out_fixation,fef_fr_in_visual,fef_fr_out_visual,fef_fr_in_memory,fef_fr_out_memory,fef_fr_in_all,fef_fr_out_all];
hold on;
for i = 1:length(baryf)
    bar(i,baryf(i));
end
hold off;
grid on;
ylabel('firing rate');
title('firing rate for FEF neuron');
legend('fixation in','fixation out','visual in','visual out','memory in','memory out','all in', 'all out', 'Location','eastoutside');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Q3 - ISI process
% ===> select multi or single here
% 1 for multi, 2 for single
isi_trials_all = [];
for i = 1:421
    this_isi = ISIcalc(spike{neuron_no,SoM}.mgs(i, :), 1,4001);
    isi_trials_all = cat(2,isi_trials_all, this_isi);
end
isi_trials_fix = [];
for i = 1:421
    this_isi = ISIcalc(spike{neuron_no,SoM}.mgs(i, :), 300,900);
    isi_trials_fix = cat(2,isi_trials_fix, this_isi);
end
isi_trials_vis = [];
for i = 1:421
    this_isi = ISIcalc(spike{neuron_no,SoM}.mgs(i, :), 1000,1600);
    isi_trials_vis = cat(2,isi_trials_vis, this_isi);
end
isi_trials_mem = [];
for i = 1:421
    this_isi = ISIcalc(spike{neuron_no,SoM}.mgs(i, :), 2400,3000);
    isi_trials_mem = cat(2,isi_trials_mem, this_isi);
end

fefisi_trials_all = [];
for i = 1:421
    this_isi = ISIcalc(spike{17,SoM}.mgs(i, :), 1,4001);
    fefisi_trials_all = cat(2,fefisi_trials_all, this_isi);
end
fefisi_trials_fix = [];
for i = 1:421
    this_isi = ISIcalc(spike{17,SoM}.mgs(i, :), 300,900);
    fefisi_trials_fix = cat(2,fefisi_trials_fix, this_isi);
end
fefisi_trials_vis = [];
for i = 1:421
    this_isi = ISIcalc(spike{17,SoM}.mgs(i, :), 1000,1600);
    fefisi_trials_vis = cat(2,fefisi_trials_vis, this_isi);
end
fefisi_trials_mem = [];
for i = 1:421
    this_isi = ISIcalc(spike{17,SoM}.mgs(i, :), 2400,3000);
    fefisi_trials_mem = cat(2,fefisi_trials_mem, this_isi);
end

%% Q3 - ISI plot
tiledlayout(4,2);
sgtitle('ISI Histograms');

nexttile;
histogram(isi_trials_fix, 'BinWidth', 1,'EdgeColor','none');
xlim([0 600]);
xlabel('mS');
ylabel('Prevalence');
title("V4 Neuron's Fixation ISI Histogram");
nexttile;
histogram(fefisi_trials_fix, 'BinWidth', 1,'EdgeColor','none');
xlim([0 100]);
xlabel('mS');
ylabel('Prevalence');
title("FEF Neuron's Fixation ISI Histogram");

nexttile;
histogram(isi_trials_vis, 'BinWidth', 1,'EdgeColor','none');
xlim([0 600]);
xlabel('mS');
ylabel('Prevalence');
title("V4 Neuron's Visual ISI Histogram");
nexttile;
histogram(fefisi_trials_vis, 'BinWidth', 1,'EdgeColor','none');
xlim([0 100]);
xlabel('mS');
ylabel('Prevalence');
title("FEF Neuron's Visual ISI Histogram");

nexttile;
histogram(isi_trials_mem, 'BinWidth', 1,'EdgeColor','none');
xlim([0 600]);
xlabel('mS');
ylabel('Prevalence');
title("V4 Neuron's Memory ISI Histogram");
nexttile;
histogram(fefisi_trials_mem, 'BinWidth', 1,'EdgeColor','none');
xlim([0 100]);
xlabel('mS');
ylabel('Prevalence');
title("FEF Neuron's Memory ISI Histogram");

nexttile;
histogram(isi_trials_all, 'BinWidth', 1,'EdgeColor','none');
xlim([0 600]);
xlabel('mS');
ylabel('Prevalence');
title("V4 Neuron's Entire ISI Histogram");
nexttile;
histogram(fefisi_trials_all, 'BinWidth', 1,'EdgeColor','none');
xlim([0 100]);
xlabel('mS');
ylabel('Prevalence');
title("FEF Neuron's Entire ISI Histogram");

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Q5 - plotting mean average rates
mean_trials =[mean_trials_led51, mean_trials_led52, mean_trials_led53, ...
    mean_trials_led54, mean_trials_led55, mean_trials_led56, ...
    mean_trials_led57, mean_trials_led58];

% figure;
% scatter(51:58, mean_trials, 100, 'red', '*');
% legend('Mean Firing Rate');
% xlim([50,59]);

%% Q5 - fitting mean average rates
fit_coef = polyfit(51:58, mean_trials, 2);
fit_out = fit_coef(1)*((50:0.1:59).^2) + fit_coef(2)*(50:0.1:59) + fit_coef(3);

figure;
hold on;
plot(51:58,mean_trials, '*r');
plot(50:0.1:59,fit_out, 'b');
legend('Mean Firing Rate', 'Fitted Quadratic Equation');
xlim([50,59]);
xlabel('LED NO.')
hold on;

%% FUNCTIONS
function windowed = window(timeser, wind)
    total = ceil(length(timeser)/wind);
    dil_time = zeros(1, total);
    for i = 1:total-1
        dil_time(i) = sum(timeser((i-1)*wind+1:i*wind));
    end
    dil_time(total) = sum(timeser((total-1)*wind+1:end));
    
    windowed = imresize(dil_time, [1,wind*total], 'nearest');
    windowed = windowed(1:length(timeser));
end

function ISI = ISIcalc(series, t1, t2)
    series = series(t1:t2);
    series(series~=0) = 1;
    series = series .* (1:length(series));
    series(series==0) = [];
    series(2:end) = series(2:end)-series(1:end-1);
    % series(1) = []; % empty vector error 
    
    ISI = series;
    if length(ISI)~=0
        ISI(1) = [];
    end
    
    
    % ISI = series; % empty vector error
end