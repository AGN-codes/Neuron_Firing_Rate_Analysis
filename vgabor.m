%% vgabor
clc;
clear all;
close all;

%% loading the data
load("data.mat");
%% neuron selection and Single or Multi
neuron_no = 15;
SoM = 1; % 1 = multi, 2 = single

%% Q5 - extracting data
trial_ang = Event_vgabor.codes(:,5) - 50; % 51:66
all_trials = spike_vgabor{neuron_no,SoM};

ang_sum_trial = zeros(1,16);
for i = 1:129
    this_ang = trial_ang(i);
    ang_sum_trial(this_ang) = ang_sum_trial(this_ang) + ...
        sum(all_trials(i, :));
end
ang_sum_trial = ang_sum_trial / 3001;
for i = 1:16
    ang_sum_trial(i) = ang_sum_trial(i) / ...
        sum(trial_ang == i);
end

%% Q5 - plotting the  data
% figure;
% scatter(51:66, ang_sum_trial, 100, 'red', '*');
% legend('Mean Firing Rate');
% xlim([51, 66]);

%% Q5 - fitting mean average rates
fit_coef = polyfit(51:66, ang_sum_trial, 2);
fit_out = fit_coef(1)*((51:0.1:66).^2) + fit_coef(2)*(51:0.1:66) + fit_coef(3);

figure;
hold on;
plot(51:66, ang_sum_trial, '*r');
plot(51:0.1:66, fit_out, 'b');
legend('Mean Firing Rate', 'Fitted Quadratic Equation');
xlim([50, 67]);
xlabel('LED NO.')
hold on;