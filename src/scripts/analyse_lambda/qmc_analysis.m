clear all

load qmene_av.dat
load mmene_av.dat

Q=qmene_av;
M=mmene_av;

q0

Q(:,2)-=Q0+BSCORR;

QMC=[Q(:,1) Q(:,2)+M(:,2)  Q(:,2) M(:,2)];

save QMC_fin.dat QMC
