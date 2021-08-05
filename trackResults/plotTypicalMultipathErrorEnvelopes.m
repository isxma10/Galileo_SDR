close all;
clear all;

numPoints = 1000;

t=linspace(0, 2*pi, numPoints);

PSK_1 = [zeros(1,length(t)) zeros(1,length(t)) ones(1,length(t)) zeros(1,length(t)) zeros(1,length(t))];
BOC_1_1 = [zeros(1,length(t)) zeros(1,length(t)) sign(sin(t)) zeros(1,length(t)) zeros(1,length(t))];
BOC_6_1 = [zeros(1,length(t)) zeros(1,length(t)) sign(sin(6*t)) zeros(1,length(t)) zeros(1,length(t))];
CBOCplus = (sqrt(1/11)*BOC_6_1 + sqrt(10/11)*BOC_1_1) ;
CBOCminus = (sqrt(1/11)*BOC_6_1 - sqrt(10/11)*BOC_1_1) ;

corr_PSK_1 = zeros(1,length(BOC_1_1));
corr_BOC_1_1 = zeros(1,length(BOC_1_1));
corr_CBOCplus = zeros(1,length(BOC_1_1));
corr_CBOCminus = zeros(1,length(BOC_1_1));

for i=1:length(BOC_1_1)
    corr_PSK_1(i)=sum(PSK_1.*circshift(PSK_1,i));
    corr_BOC_1_1(i)=sum(BOC_1_1.*circshift(BOC_1_1,i));
    corr_CBOCplus(i)=sum(CBOCplus.*circshift(CBOCplus,i));
    corr_CBOCminus(i)=sum(CBOCminus.*circshift(CBOCminus,i));
    
end



% PSK multi path error
delta1 = round(numPoints);
shiftPSK1 = circshift(corr_PSK_1, delta1);
discPSK1 = shiftPSK1 - corr_PSK_1;

delta2 = round(numPoints/2);
shiftPSK2 = circshift(corr_PSK_1, delta2);
discPSK2 = shiftPSK2 - corr_PSK_1;

codedelays = linspace(-2.5, 2.5, 5*numPoints);

multiPath1 = discPSK1 - 0.5*circshift(discPSK1, round(0.75*numPoints));


figure(1)
plot(codedelays, circshift(corr_PSK_1, (length(corr_PSK_1)/2)-numPoints/2)/numPoints,'LineWidth',2)
hold on
plot(codedelays, circshift(shiftPSK1, (length(corr_PSK_1)/2)-numPoints/2)/numPoints,'LineWidth',2)

plot(codedelays, circshift(discPSK1, (length(corr_PSK_1)/2)-numPoints/2)/numPoints,'LineWidth',2)
% plot(codedelays, circshift(discPSK2, (length(corr_PSK_1)/2)-numPoints/2+delta2/2)/numPoints,'LineWidth',2)

grid on
% ylim([-0.2 1.2]);
xlabel('Code delay (chips)');
ylabel('Normalised correlation amplitude');
figure(100)
plot(codedelays, circshift(corr_PSK_1, (length(corr_PSK_1)/2)-numPoints/2+delta2/2)/numPoints,'LineWidth',2)
hold on
plot(codedelays, circshift(shiftPSK2, (length(corr_PSK_1)/2)-numPoints/2+delta2/2)/numPoints,'LineWidth',2)

plot(codedelays, circshift(discPSK2, (length(corr_PSK_1)/2)-numPoints/2+delta2/2)/numPoints,'LineWidth',2)

grid on
ylim([-1 1]);
xlabel('Code delay (chips)');
ylabel('Normalised correlation amplitude');

figure(2)
hold on

plot(codedelays, circshift(discPSK1, (length(corr_PSK_1)/2)-numPoints/2)/numPoints,'LineWidth',2)
plot(codedelays, circshift(-0.5*circshift(discPSK1, round(0.75*numPoints)), (length(corr_PSK_1)/2)-numPoints/2)/numPoints,'LineWidth',2)
plot(codedelays, circshift(multiPath1, (length(corr_PSK_1)/2)-numPoints/2)/numPoints,'LineWidth',2)

grid on
% ylim([-0.2 1.2]);
xlabel('Code delay (chips)');
ylabel('Normalised correlation amplitude');
legend('Direct signal','Multipath 180° phase','Discriminator')

multipathIndexPlus = ones(1,round(numPoints+delta1/2));
multipathIndexMinus = ones(1,round(numPoints+delta1/2));

for n = 0:(length(multipathIndexPlus)-1)
    multiPath1 = discPSK1 + 0.5*circshift(discPSK1, n);
    multiPath2 = discPSK1 - 0.5*circshift(discPSK1, n);
    multipathDiscPlus = circshift(multiPath1, (length(corr_PSK_1)/2)-numPoints/2)/numPoints;
    multipathDiscMinus = circshift(multiPath2, (length(corr_PSK_1)/2)-numPoints/2)/numPoints;
    
    for i = 2:length(multipathDiscPlus)
        
        if (multipathDiscPlus(i) > 0) && (multipathDiscPlus(i-1) <= 0)
            multipathIndexPlus(n+1) = i;
        end
        if (multipathDiscMinus(i) > 0) && (multipathDiscMinus(i-1) <= 0)
            multipathIndexMinus(n+1) = i;
        end
    end
end

multipathDelay = linspace(0, (numPoints+delta1)/numPoints, round(numPoints+delta1/2));

figure(3)
hold on

plot(multipathDelay, codedelays(multipathIndexPlus),'LineWidth',2)
plot(multipathDelay, codedelays(multipathIndexMinus),'LineWidth',2)
grid on
% ylim([-0.2 1.2]);
xlabel('Multipath delay (chips)');
ylabel('Multipath error (chips)');

multipathIndexPlus = ones(1,round(numPoints+delta1/2));
multipathIndexMinus = ones(1,round(numPoints+delta1/2));

for n = 0:(length(multipathIndexPlus)-1)
    multiPath1 = discPSK2 + 0.5*circshift(discPSK2, n);
    multiPath2 = discPSK2 - 0.5*circshift(discPSK2, n);
    multipathDiscPlus = circshift(multiPath1, (length(corr_PSK_1)/2)-round(delta1/2))/numPoints;
    multipathDiscMinus = circshift(multiPath2, (length(corr_PSK_1)/2)-round(delta1/2))/numPoints;
    
    for i = 2*numPoints:3*numPoints
        
        if (multipathDiscPlus(i) > 0) && (multipathDiscPlus(i-1) <= 0)
            multipathIndexPlus(n+1) = i;
        end
        if (multipathDiscMinus(i) > 0) && (multipathDiscMinus(i-1) <= 0)
            multipathIndexMinus(n+1) = i;
        end
    end
end

multipathDelay = linspace(0, (numPoints+delta1)/numPoints, round(numPoints+delta1/2));

figure(7)
hold on

plot(multipathDelay, codedelays(multipathIndexPlus),'LineWidth',2)
plot(multipathDelay, codedelays(multipathIndexMinus),'LineWidth',2)
grid on
% ylim([-0.2 1.2]);
xlabel('Multipath delay (chips)');
ylabel('Multipath error (chips)');

% BOC 11 multi path error
delta1 = round(numPoints/10);
multipathDelay = linspace(0, (numPoints+delta1/2)/numPoints, round(numPoints+delta1/2));

shiftBOC11 = circshift(corr_BOC_1_1, delta1);
discBOC11 = shiftBOC11 - corr_BOC_1_1;

multipathIndexPlus = ones(1,round(numPoints+delta1/2));
multipathIndexMinus = ones(1,round(numPoints+delta1/2));

for n = 0:(length(multipathIndexPlus)-1)
    multiPath1 = discBOC11 + 0.5*circshift(discBOC11, n);
    multiPath2 = discBOC11 - 0.5*circshift(discBOC11, n);
    multipathDiscPlus = circshift(multiPath1, (length(corr_PSK_1)/2)-round(delta1/2))/numPoints;
    multipathDiscMinus = circshift(multiPath2, (length(corr_PSK_1)/2)-round(delta1/2))/numPoints;
    
    for i = 2*numPoints:3*numPoints
        
        if (multipathDiscPlus(i) > 0) && (multipathDiscPlus(i-1) <= 0)
            multipathIndexPlus(n+1) = i;
        end
        if (multipathDiscMinus(i) > 0) && (multipathDiscMinus(i-1) <= 0)
            multipathIndexMinus(n+1) = i;
        end
    end
end

figure(4)
hold on

plot(multipathDelay, codedelays(multipathIndexPlus),'LineWidth',2)
plot(multipathDelay, codedelays(multipathIndexMinus),'LineWidth',2)
grid on
% ylim([-0.2 1.2]);
xlabel('Multipath delay (chips)');
ylabel('Multipath error (chips)');

% CBOC multi path error
delta1 = round(numPoints/10);
multipathDelay = linspace(0, (numPoints+delta1/2)/numPoints, round(numPoints+delta1/2));

shiftCBOCm = circshift(corr_CBOCminus, delta1);
discCBOCm = shiftCBOCm - corr_CBOCminus;
shiftCBOCp = circshift(corr_CBOCplus, delta1);
discCBOCp = shiftCBOCp - corr_CBOCplus;

multipathIndexPlus = ones(1,round(numPoints+delta1/2));
multipathIndexMinus = ones(1,round(numPoints+delta1/2));

for n = 0:(length(multipathIndexPlus)-1)
    multiPath1 = discCBOCm + 0.5*circshift(discCBOCm, n);
    multiPath2 = discCBOCm - 0.5*circshift(discCBOCm, n);
    multipathDiscPlus = circshift(multiPath1, (length(corr_PSK_1)/2)-round(delta1/2))/numPoints;
    multipathDiscMinus = circshift(multiPath2, (length(corr_PSK_1)/2)-round(delta1/2))/numPoints;
    
    for i = 2.45*numPoints:2.55*numPoints
        
        if (multipathDiscPlus(i) > 0) && (multipathDiscPlus(i-1) <= 0)
            multipathIndexPlus(n+1) = i;
        end
        if (multipathDiscMinus(i) > 0) && (multipathDiscMinus(i-1) <= 0)
            multipathIndexMinus(n+1) = i;
        end
    end
end

figure(4)
hold on

plot(multipathDelay, codedelays(multipathIndexPlus),'LineWidth',2)
plot(multipathDelay, codedelays(multipathIndexMinus),'LineWidth',2)
grid on
% ylim([-0.2 1.2]);
xlabel('Multipath delay (chips)');
ylabel('Multipath error (chips)');

% do for CBOC pluse
multipathIndexPlus = ones(1,round(numPoints+delta1/2));
multipathIndexMinus = ones(1,round(numPoints+delta1/2));

for n = 0:(length(multipathIndexPlus)-1)
    multiPath1 = discCBOCp + 0.5*circshift(discCBOCp, n);
    multiPath2 = discCBOCp - 0.5*circshift(discCBOCp, n);
    multipathDiscPlus = circshift(multiPath1, (length(corr_PSK_1)/2)-round(delta1/2))/numPoints;
    multipathDiscMinus = circshift(multiPath2, (length(corr_PSK_1)/2)-round(delta1/2))/numPoints;
    
    for i = 2.45*numPoints:2.55*numPoints
        
        if (multipathDiscPlus(i) > 0) && (multipathDiscPlus(i-1) <= 0)
            multipathIndexPlus(n+1) = i;
        end
        if (multipathDiscMinus(i) > 0) && (multipathDiscMinus(i-1) <= 0)
            multipathIndexMinus(n+1) = i;
        end
    end
end

figure(4)
hold on

plot(multipathDelay, codedelays(multipathIndexPlus),'LineWidth',2)
plot(multipathDelay, codedelays(multipathIndexMinus),'LineWidth',2)
grid on
% ylim([-0.2 1.2]);
xlabel('Multipath delay (chips)');
ylabel('Multipath error (chips)');

legend('BOC(1,1) 0°','BOC(1,1) 180°','CBOC(6,1,1/11,-) 0°','CBOC(6,1,1/11,-) 180°','CBOC(6,1,1/11,+) 0°','CBOC(6,1,1/11,+) 180°')
