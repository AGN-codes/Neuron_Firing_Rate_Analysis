%% vodr
clc;
clear all;
close all;

%% loading the data
load("data.mat");

%% neuron selection and Single or Multi
neuron_no = 15;
SoM = 1; % 1 = multi, 2 = single

%% Q4 - getting the led numbers and all trials
led_no = [Event.vodr.codes(:,5), Event.vodr.codes(:,7), ...
    Event.vodr.codes(:,9), Event.vodr.codes(:,11), ...
    Event.vodr.codes(:,13), Event.vodr.codes(:,15), ...
    Event.vodr.codes(:,17), Event.vodr.codes(:,19)]-50;

all_trials = spike{neuron_no,SoM}.vodr;

%% Q4 - extracting useful data for RF mapping task
eff_spk = all_trials(:, 252:end);
eff_idx = [1:100, 201:300, 401:500, 601:700, ...
    801:900, 1001:1100, 1201:1300, 1401:1500];
eff_spk = eff_spk(:,eff_idx);

total_spk = zeros(304, 8);
for i = 1:304
    for j = 1:8
        total_spk(i,j) = sum(eff_spk(i, (((j-1)*100)+1):j*100))/100;
    end
end

led_spk = zeros(1, 49);
for i = 1:304
    for j = 1:8
        led_spk(led_no(i,j)) = ...
            led_spk(led_no(i,j)) + total_spk(i,j);
    end
end

for i = 1:49
    led_spk(i) = led_spk(i) / sum( led_no == i, 'all');
end

led_heatmap = reshape(led_spk,7,7);

%% Q4 - plotting the RF mapping task
figure;
led_heatmap_heatmap = heatmap(led_heatmap);
ylabel('Average Firing Rate (Hz)');
title('Receptive Field Heat Map');
