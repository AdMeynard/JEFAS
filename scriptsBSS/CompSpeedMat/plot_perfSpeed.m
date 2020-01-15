close all; clear all; clc;
addpath('results/')

load results1
sSIR0(1) = mean(SIR0);
sSIRsobi(1) = mean(SIRsobi);
sSIRpsobi(1) = mean(SIRpsobi);
sSIRqtf(1) = mean(SIRqtf);
sSIR(1) = mean(SIR);
sAMsobi(1) = mean(indSOBI);
sAMpsobi(1) = mean(indPSOBI);
sAMqtf(1) = mean(indQTF);
sAM(1) = mean(indJEFAS);

load results2
sSIR0(2) = mean(SIR0);
sSIRsobi(2) = mean(SIRsobi);
sSIRpsobi(2) = mean(SIRpsobi);
sSIRqtf(2) = mean(SIRqtf);
sSIR(2) = mean(SIR);
sAMsobi(2) = mean(indSOBI);
sAMpsobi(2) = mean(indPSOBI);
sAMqtf(2) = mean(indQTF);
sAM(2) = mean(indJEFAS);

load results3
sSIR0(3) = mean(SIR0);
sSIRsobi(3) = mean(SIRsobi);
sSIRpsobi(3) = mean(SIRpsobi);
sSIRqtf(3) = mean(SIRqtf);
sSIR(3) = mean(SIR);
sAMsobi(3) = mean(indSOBI);
sAMpsobi(3) = mean(indPSOBI);
sAMqtf(3) = mean(indQTF);
sAM(3) = mean(indJEFAS);

load results5
sSIR0(4) = mean(SIR0);
sSIRsobi(4) = mean(SIRsobi);
sSIRpsobi(4) = mean(SIRpsobi);
sSIRqtf(4) = mean(SIRqtf);
sSIR(4) = mean(SIR);
sAMsobi(4) = mean(indSOBI);
sAMpsobi(4) = mean(indPSOBI);
sAMqtf(4) = mean(indQTF);
sAM(4) = mean(indJEFAS);

load results10
sSIR0(5) = mean(SIR0);
sSIRsobi(5) = mean(SIRsobi);
sSIRpsobi(5) = mean(SIRpsobi);
sSIRqtf(5) = mean(SIRqtf);
sSIR(5) = mean(SIR);
sAMsobi(5) = mean(indSOBI);
sAMpsobi(5) = mean(indPSOBI);
sAMqtf(5) = mean(indQTF);
sAM(5) = mean(indJEFAS);

load results20
sSIR0(6) = mean(SIR0);
sSIRsobi(6) = mean(SIRsobi);
sSIRpsobi(6) = mean(SIRpsobi);
sSIRqtf(6) = mean(SIRqtf);
sSIR(6) = mean(SIR);
sAMsobi(6) = mean(indSOBI);
sAMpsobi(6) = mean(indPSOBI);
sAMqtf(6) = mean(indQTF);
sAM(6) = mean(indJEFAS);

load results50
sSIR0(7) = mean(SIR0);
sSIRsobi(7) = mean(SIRsobi);
sSIRpsobi(7) = mean(SIRpsobi);
sSIRqtf(7) = mean(SIRqtf);
sSIR(7) = mean(SIR);
sAMsobi(7) = mean(indSOBI);
sAMpsobi(7) = mean(indPSOBI);
sAMqtf(7) = mean(indQTF);
sAM(7) = mean(indJEFAS);

load results100
sSIR0(8) = mean(SIR0);
sSIRsobi(8) = mean(SIRsobi);
sSIRpsobi(8) = mean(SIRpsobi);
sSIRqtf(8) = mean(SIRqtf);
sSIR(8) = mean(SIR);
sAMsobi(8) = mean(indSOBI);
sAMpsobi(8) = mean(indPSOBI);
sAMqtf(8) = mean(indQTF);
sAM(8) = mean(indJEFAS);

s = [1 2 3 5 10 20 50 100];
figure;
subplot(1,2,1);
semilogx(s,sSIRsobi,'k--+',s,sSIRpsobi,'r:+',s,sSIRqtf,'g-.+',s,sSIR,'b--+','linewidth',2); 
xlabel('Normalized $\|\mathbf{A}''\|_\infty$','interpreter','latex'); ylabel('SIR (dB)'); 
xticks(s);
axis tight; grid on;
legend({'SOBI','p-SOBI','QTF-BSS','JEFAS-BSS'},'Orientation','horizontal'); set(gca,'fontsize',24);

subplot(1,2,2);
semilogx(s,sAMsobi,'k--+',s,sAMpsobi,'r:+',s,sAMqtf,'g-.+',s,sAM,'b--+','linewidth',2); 
xlabel('Normalized $\|\mathbf{A}''\|_\infty$','interpreter','latex'); ylabel('Averaged $\rho$ (dB)','interpreter','latex'); 
xticks(s);
axis tight; grid on;
legend({'SOBI','p-SOBI','QTF-BSS','JEFAS-BSS'},'Orientation','horizontal'); set(gca,'fontsize',24);