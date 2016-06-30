% Plot file

close all
clear all
toto = load('./cvg_kmf.mat');
error_kmf  = toto.error2;
toto = load('./cvg_kmfo.mat');
error_kmfo  = toto.error;
toto = load('./cvg_sp.mat');
error_sp  = toto.error;
toto = load('./cvg_sppp.mat');
error_sppp  = toto.error;
toto = load('./cvg_spd.mat');
error_spd  = toto.error;

hold on
set(gca, 'fontsize', 15);
plot(log10(error_kmf),'Color','blue')
plot(log10(error_kmfo(2:end)),'Color','black')
plot(log10(error_sp(2:end)),'Color','red')
plot(log10(error_sppp(2:end)),'Color','green')
plot(log10(error_spd(2:end)),'Color','magenta')
legend('kmf','kmf+orthodir','SP','preconditionned SP','DSP')