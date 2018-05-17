%FUNDnational
%The Climate Framework for Uncertainty, Negotiation and Distribution
%version 4.0-matlab-national
%
%This script runs FUND 4.0 MN
%
%Richard Tol, 2 September 2014
%This code is protected by the MIT License

Initialize
SRES
ReadData
SetParameters
InitModules
Calibration

for t=2:NHistYear
     [MRHbox(t,1,:), CO2conc(t,1), pH(t,1)] = stepCO2(MRHbox(t-1,1,:),historicCO2emit(t-1),historicLUemit(t-1));
     [CH4conc(t,1), N2Oconc(t,1), SF6conc(t,1), CFC11conc(t,1), CFC12conc(t,1)] = stepGHG(CH4conc(t-1,1),historicCH4emit(t-1),N2Oconc(t-1,1),historicN2Oemit(t-1),SF6conc(t-1,1),historicSF6emit(t-1),CFC11conc(t-1,1),histCFC11emit(t-1),CFC12conc(t-1,1),histCFC12emit(t-1));
     CFC12conc(t,1) = (1-CFC12life)*CFC12conc(t-1,1) + histCFC12emit(t-1);
     [RadForc(t,1),atmtemp(t,1), oceantemp(t,1),SLR(t,1)] = stepClimate(CO2conc(t,1),CH4conc(t,1),N2Oconc(t,1),SF6conc(t,1),CFC11conc(t,1),CFC12conc(t,1),Semit(t,1),trO3radforc(t,1),atmtemp(t-1,1),oceantemp(t-1,1));
     [K(t,1),Y(t,1),Energy(t,1),CO2emit(t,1)] = stepEconomy(K(t-1,1),Y(t-1,1),TFP(t,1),histPopulation(t),EnInt(t,1),CO2Int(t,1));
     impact(:,t,1) = aggimpact(atmtemp(t,1),imppar,Y(t,1),Population(t,1),YpC2010);
     for c=1:NCountry
            impctr(:,c,t,1) = aggimpact(atmtemp(t,1),squeeze(ctrimppar(:,c,:)),YCtr(c,t,1),PopCtr(c,t,1),YpC2010Ctr(c));
     end
end

for s=2:NScen
    MRHbox(:,s,:) = MRHbox(:,1,:);
    CO2conc(:,s) = CO2conc(:,1);
    CH4conc(:,s) = CH4conc(:,1);
    N2Oconc(:,s) = N2Oconc(:,1);
    SF6conc(:,s) = SF6conc(:,1);
    CFC11conc(:,s) = CFC11conc(:,1);
    CFC12conc(:,s) = CFC12conc(:,1);
    pH(:,s) = pH(:,1);
    RadForc(:,s) = RadForc(:,1);
    atmtemp(:,s) = atmtemp(:,1);
    oceantemp(:,s) = oceantemp(:,1);
    SLR(:,s) = SLR(:,1);
    K(:,s) = K(:,1);
    Y(:,s) = Y(:,1);
    Energy(:,s) = Energy(:,1);
    CO2emit(:,s) = CO2emit(:,1);
    impact(:,:,s) = impact(:,:,1);
    impctr(:,:,:,s) = impctr(:,:,:,1);
end

impactd = impact;

for t=NHistYear+1:NYear
    for s=1:NScen
        if t == 265
            [MRHbox(t,s,:), CO2conc(t,s), pH(t,s)]= stepCO2(MRHbox(t-1,s,:),CO2emit(t-1,s)+1,LUemit(t-1,s));
        else
            [MRHbox(t,s,:), CO2conc(t,s), pH(t,s)]= stepCO2(MRHbox(t-1,s,:),CO2emit(t-1,s),LUemit(t-1,s));
        end
        [CH4conc(t,s), N2Oconc(t,s), SF6conc(t,s),CFC11conc(t,s),CFC12conc(t,s)] = stepGHG(CH4conc(t-1,s),CH4emit(t-1,s),N2Oconc(t-1,s),N2Oemit(t-1,s),SF6conc(t-1,s),SF6emit(t-1,s),CFC11conc(t-1,s),CFC11emit(t-1,s),CFC12conc(t-1,s),CFC12emit(t-1,s));
        [RadForc(t,s),atmtemp(t,s),oceantemp(t,s),SLR(t,s)] = stepClimate(CO2conc(t,s),CH4conc(t,s),N2Oconc(t,s),SF6conc(t,s),CFC11conc(t,s),CFC12conc(t,s),Semit(t,s),trO3radforc(t,s),atmtemp(t-1,s),oceantemp(t-1,s));
        [K(t,s),Y(t,s),Energy(t,s),CO2emit(t,s)]= stepEconomy(K(t-1,s),Y(t-1,s),TFP(t,s),Population(t,s),EnInt(t,s),CO2Int(t,s));
        for c=1:NCountry
             [KCtr(c,t,s),YCtr(c,t,s),EnergyCtr(c,t,s),CO2Ctr(c,t,s)]= stepEconomy(KCtr(c,t-1,s),YCtr(c,t-1,s),TFPCtr(c,t,s),PopCtr(c,t,s),EnIntCtr(c,t,s),CO2IntCtr(c,t,s));
        end
        impactd(:,t,s) = aggimpact(atmtemp(t,s),imppar,Y(t,s),Population(t,s),YpC2010);
        for c=1:NCountry
            impctrd(:,c,t,s) = aggimpact(atmtemp(t,s),squeeze(ctrimppar(:,c,:)),YCtr(c,t,s),PopCtr(c,t,s),YpC2010Ctr(c));
        end
    end
end

for t=NHistYear+1:NYear
    for s=1:NScen
        [MRHbox(t,s,:), CO2conc(t,s), pH(t,s)]= stepCO2(MRHbox(t-1,s,:),CO2emit(t-1,s),LUemit(t-1,s));
        [CH4conc(t,s), N2Oconc(t,s), SF6conc(t,s),CFC11conc(t,s),CFC12conc(t,s)] = stepGHG(CH4conc(t-1,s),CH4emit(t-1,s),N2Oconc(t-1,s),N2Oemit(t-1,s),SF6conc(t-1,s),SF6emit(t-1,s),CFC11conc(t-1,s),CFC11emit(t-1,s),CFC12conc(t-1,s),CFC12emit(t-1,s));
        [RadForc(t,s),atmtemp(t,s),oceantemp(t,s),SLR(t,s)] = stepClimate(CO2conc(t,s),CH4conc(t,s),N2Oconc(t,s),SF6conc(t,s),CFC11conc(t,s),CFC12conc(t,s),Semit(t,s),trO3radforc(t,s),atmtemp(t-1,s),oceantemp(t-1,s));
        [K(t,s),Y(t,s),Energy(t,s),CO2emit(t,s)]= stepEconomy(K(t-1,s),Y(t-1,s),TFP(t,s),Population(t,s),EnInt(t,s),CO2Int(t,s));
        for c=1:NCountry
             [KCtr(c,t,s),YCtr(c,t,s),EnergyCtr(c,t,s),CO2Ctr(c,t,s)]= stepEconomy(KCtr(c,t-1,s),YCtr(c,t-1,s),TFPCtr(c,t,s),PopCtr(c,t,s),EnIntCtr(c,t,s),CO2IntCtr(c,t,s));
        end
        impact(:,t,s) = aggimpact(atmtemp(t,s),imppar,Y(t,s),Population(t,s),YpC2010);
        for c=1:NCountry
            impctr(:,c,t,s) = aggimpact(atmtemp(t,s),squeeze(ctrimppar(:,c,:)),YCtr(c,t,s),PopCtr(c,t,s),YpC2010Ctr(c));
        end
    end
end

for i=1:NImpact,
    for t=1:NYear,
        for s=1:NScen,
            imptot(i,t,s) = impctr(i,:,t,s)*YCtr(:,t,s) / sum(YCtr(:,t,s));
        end
    end
end
     
YpC(:,:) = Y(:,:)./Population(:,:);
YpCCtr(:,:,:) = YCtr(:,:,:)./PopCtr(:,:,:);

SocialCostofCarbon;

[YpC2010s Index] = sort(YpCCtr(:,NHistYear,1));
YpC2100s = YpCCtr(:,NYear,1);
YpC2100s = YpC2100s(Index);
Imp1s = squeeze(impctr(1,:,NYear,1));
Imp2s = squeeze(impctr(2,:,NYear,1));
Imp3s = squeeze(impctr(3,:,NYear,1));
Imp4s = squeeze(impctr(4,:,NYear,1));
Imp5s = squeeze(impctr(5,:,NYear,1));
Imp6s = squeeze(impctr(6,:,NYear,1));
Imp7s = squeeze(impctr(7,:,NYear,1));
Imp8s = squeeze(impctr(8,:,NYear,1));
Imp1s = Imp1s(Index);
Imp2s = Imp2s(Index);
Imp3s = Imp3s(Index);
Imp4s = Imp4s(Index);
Imp5s = Imp5s(Index);
Imp6s = Imp6s(Index);
Imp7s = Imp7s(Index);
Imp8s = Imp8s(Index);
SCC1 = squeeze(SCCc(1,:,1,4,2));
SCC2 = squeeze(SCCc(2,:,1,4,2));
SCC3 = squeeze(SCCc(3,:,1,4,2));
SCC4 = squeeze(SCCc(4,:,1,4,2));
SCC5 = squeeze(SCCc(5,:,1,4,2));
SCC6 = squeeze(SCCc(6,:,1,4,2));
SCC7 = squeeze(SCCc(7,:,1,4,2));
SCC8 = squeeze(SCCc(8,:,1,4,2));
SCC1 = SCC1(Index);
SCC2 = SCC2(Index);
SCC3 = SCC3(Index);
SCC4 = SCC4(Index);
SCC5 = SCC5(Index);
SCC6 = SCC6(Index);
SCC7 = SCC7(Index);
SCC8 = SCC8(Index);

subplot(4,7,1), plot(Year,Population(:,1),Year,Population(:,2),Year,Population(:,3),Year,Population(:,4),Year,Population(:,5),Year,Population(:,6),Year,Population(:,7),Year,Population(:,8),Year,Population(:,9)), xlabel('year'), ylabel('number of people'), title('Population')
subplot(4,7,2), plot(Year,YpC(:,1),Year,YpC(:,2),Year,YpC(:,3),Year,YpC(:,4),Year,YpC(:,5),Year,YpC(:,6),Year,YpC(:,7),Year,YpC(:,8),Year,YpC(:,9)), xlabel('year'), ylabel('dollar per person per year'), title('Average income')
subplot(4,7,3), plot(Year,EnInt(:,1)*10^12,Year,EnInt(:,2)*10^12,Year,EnInt(:,3)*10^12,Year,EnInt(:,4)*10^12,Year,EnInt(:,5)*10^12,Year,EnInt(:,6)*10^12,Year,EnInt(:,7)*10^12,Year,EnInt(:,8)*10^12,Year,EnInt(:,9)*10^12), xlabel('year'), ylabel('gram oil equivalent per dollar'), title('Energy intensity')
subplot(4,7,4), plot(Year,CO2Int(:,1),Year,CO2Int(:,2),Year,CO2Int(:,3),Year,CO2Int(:,4),Year,CO2Int(:,5),Year,CO2Int(:,6),Year,CO2Int(:,7),Year,CO2Int(:,8),Year,CO2Int(:,9)), xlabel('year'), ylabel('gram carbon per gram oil equivalent'), title('Carbon intensity')
subplot(4,7,5), plot(Year,CO2emit(:,1)/1000,Year,CO2emit(:,2)/1000,Year,CO2emit(:,3)/1000,Year,CO2emit(:,4)/1000,Year,CO2emit(:,5)/1000,Year,CO2emit(:,6)/1000,Year,CO2emit(:,7)/1000,Year,CO2emit(:,8)/1000,Year,CO2emit(:,9)/1000), xlabel('year'), ylabel('billion tonnes of carbon'), title('Carbon dioxide emissions')
subplot(4,7,6), plot(Year,CO2conc(:,1),Year,CO2conc(:,2),Year,CO2conc(:,3),Year,CO2conc(:,4),Year,CO2conc(:,5),Year,CO2conc(:,6),Year,CO2conc(:,7),Year,CO2conc(:,8),Year,CO2conc(:,9)), xlabel('year'), ylabel('parts per million by volume'), title('Carbon dioxide concentration')
subplot(4,7,7), plot(Year,atmtemp(:,1),Year,atmtemp(:,2),Year,atmtemp(:,3),Year,atmtemp(:,4),Year,atmtemp(:,5),Year,atmtemp(:,6),Year,atmtemp(:,7),Year,atmtemp(:,8),Year,atmtemp(:,9)), xlabel('year'), ylabel('degree Celsius'), title('Temperature')
subplot(4,7,8), plot(Year,imptot(1,:,1),Year,imptot(1,:,2),Year,imptot(1,:,3),Year,imptot(1,:,4),Year,imptot(1,:,5),Year,imptot(1,:,6),Year,imptot(1,:,7),Year,imptot(1,:,8),Year,imptot(1,:,9)), xlabel('year'), ylabel('percent income'), title('Tol (parabola)')
subplot(4,7,9), plot(Year,imptot(2,:,1),Year,imptot(2,:,2),Year,imptot(2,:,3),Year,imptot(2,:,4),Year,imptot(2,:,5),Year,imptot(2,:,6),Year,imptot(2,:,7),Year,imptot(2,:,8),Year,imptot(2,:,9)), xlabel('year'), ylabel('percent income'), title('Weitzman (6)')
subplot(4,7,10), plot(Year,imptot(3,:,1),Year,imptot(3,:,2),Year,imptot(3,:,3),Year,imptot(3,:,4),Year,imptot(3,:,5),Year,imptot(3,:,6),Year,imptot(3,:,7),Year,imptot(3,:,8),Year,imptot(3,:,9)), xlabel('year'), ylabel('percent income'), title('Weitzman (7)')
subplot(4,7,11), plot(Year,imptot(4,:,1),Year,imptot(4,:,2),Year,imptot(4,:,3),Year,imptot(4,:,4),Year,imptot(4,:,5),Year,imptot(4,:,6),Year,imptot(4,:,7),Year,imptot(4,:,8),Year,imptot(4,:,9)), xlabel('year'), ylabel('percent income'), title('Nordhaus')
subplot(4,7,12), plot(Year,imptot(5,:,1),Year,imptot(5,:,2),Year,imptot(5,:,3),Year,imptot(5,:,4),Year,imptot(5,:,5),Year,imptot(5,:,6),Year,imptot(5,:,7),Year,imptot(5,:,8),Year,imptot(5,:,9)), xlabel('year'), ylabel('percent income'), title('Hope')
subplot(4,7,13), plot(Year,imptot(6,:,1),Year,imptot(6,:,2),Year,imptot(6,:,3),Year,imptot(6,:,4),Year,imptot(6,:,5),Year,imptot(6,:,6),Year,imptot(6,:,7),Year,imptot(6,:,8),Year,imptot(6,:,9)), xlabel('year'), ylabel('percent income'), title('van der Ploeg')
subplot(4,7,14), plot(Year,imptot(8,:,1),Year,imptot(8,:,2),Year,imptot(8,:,3),Year,imptot(8,:,4),Year,imptot(8,:,5),Year,imptot(8,:,6),Year,imptot(8,:,7),Year,imptot(8,:,8),Year,imptot(8,:,9)), xlabel('year'), ylabel('percent income'), title('Tol (bilinear)')
subplot(4,7,15),bar(Imp1s)
subplot(4,7,16),bar(Imp2s)
subplot(4,7,17),bar(Imp3s)
subplot(4,7,18),bar(Imp4s)
subplot(4,7,19),bar(Imp5s)
subplot(4,7,20),bar(Imp6s)
subplot(4,7,21),bar(Imp8s)
subplot(4,7,22),bar(SCC1)
subplot(4,7,23),bar(SCC2)
subplot(4,7,24),bar(SCC3)
subplot(4,7,25),bar(SCC4)
subplot(4,7,26),bar(SCC5)
subplot(4,7,27),bar(SCC6)
subplot(4,7,28),bar(SCC8)