%% Error plots
% p = 1
p1UnOF = [11 21 41 81 161];
p1UH1 = [0.14251 0.07124 0.03564 0.01785 0.00899];

plot(p1UnOF,p1UH1)
figure
semilogx(p1UnOF,p1UH1)
figure
loglog(p1UnOF,p1UH1)


% p2
p2UnOF = [13 23 43 83 163];
p2UH1 = [0.00597729 0.001902932 0.00128554 0.001237 0.001233907 ];

plot(p2UnOF,p2UH1)
figure
semilogx(p2UnOF,p2UH1)
figure
loglog(p2UnOF,p2UH1)
