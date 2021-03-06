%initKaya
%The Climate Framework for Uncertainty, Negotiation and Distribution,
%version 4.1-matlab-national
%
%This script is part of FUND 4.1 MN
%It initializes variables and sets parameters
%
%Richard Tol, 14 May 2018
%This code is protected by the MIT License

%% global
TFP = zeros(NYear,NScen);
K = zeros(NYear,NScen);
Y = zeros(NYear,NScen);

TFP(1,1) = (histGDP(1)/histPopulation(1))^LabourElast * (Depreciation/SavingsRate)^(1-LabourElast);
K(1,1) = SavingsRate*histGDP(1)/Depreciation;
Y(1,1) = TFP(1,1)*histPopulation(1)^LabourElast*K(1,1)^(1-LabourElast);

gTFP = 0.0075;

for t=2:NHistYear
    TFP(t,1) = (1+gTFP)*TFP(t-1,1);
end

CalibTFP

for s=2:NScen
    TFP(:,s) = TFP(:,1);
end

Population= zeros(NYear,NScen);
EnInt = zeros(NYear,NScen);
Energy = zeros(NYear,NScen);
CO2Int = zeros(NYear,NScen);
for s=1:NScen
    Population(1:NHistYear,s)=histPopulation;
    EnInt(1:NHistYear,s)=histEnInt;
    CO2Int(1:NHistYear,s)=histCO2Int;
end

Energy(1,1) = EnInt(1,1)*Y(1,1);
CO2emit(1,1) = CO2Int(1,1)*Energy(1,1);

for t=NHistYear+1:NYear
    ts = t-NHistYear;
    for s=1:NSRES
        Population(t,s) = (1+SRESdPop(s,ts))*Population(t-1,s);
        TFP(t,s) = (1+SRESdInc(s,ts))*TFP(t-1,s);
        EnInt(t,s)= (1+SRESdEnInt(s,ts))*EnInt(t-1,s);
        CO2Int(t,s)= (1+SRESdCO2Int(s,ts))*CO2Int(t-1,s);
    end
    for s=NSRES+1:NScen
        Population(t,s) = (1+SSPdPop(s-NSRES,ts))*Population(t-1,s);
        TFP(t,s) = (1+SSPdInc(s-NSRES,ts))*TFP(t-1,s);
        EnInt(t,s)= (1+SSPdEnInt(s-NSRES,ts))*EnInt(t-1,s);
        CO2Int(t,s)= (1+SSPdCO2Int(s-NSRES,ts))*CO2Int(t-1,s);
    end
end

%% national

%population
xdPop = histPopCtr(:,NHistCtrYr)./histPopCtr(:,1);
xdPop = xdPop.^(1/NHistCtrYr);
xdPopave = sum(xdPop)/NCountry;
PopCtr = zeros(NCountry,NYear,NScen);

for s=1:NScen
    PopCtr(:,NHistYear-NHistCtrYr+1:NHistYear,s) = histPopCtr;
end

for t=NHistYear+1:NYear
    xweight = ((NYear-t)/(NYear-NHistYear))^2;
    for s=1:NScen
        PopCtr(:,t,s) = PopCtr(:,t-1,s).*(xdPop*xweight + xdPopave*(1-xweight));
        xnorm = sum(PopCtr(:,t,s));
        PopCtr(:,t,s) = PopCtr(:,t,s)*Population(t,s)/xnorm;
    end
end

%economic output
TFPCtr = zeros(NCountry,NYear,NScen);
KCtr = zeros(NCountry,NYear,NScen);
YCtr = zeros(NCountry,NYear,NScen);
YpCCtr = zeros(NCountry,NYear,NScen);

%historical values
for t = NHistYear-NHistCtrYr+1:NHistYear,
    ts = t - NHistYear + NHistCtrYr;
    for c = 1:NCountry,
        TFPCtr(c,t,1) = (histGDPCtr(c,ts)/histPopCtr(c,ts))^LabourElast * (Depreciation/SavingsRate)^(1-LabourElast);
        KCtr(c,t,1) = SavingsRate*histGDPCtr(c,ts)/Depreciation;
        YCtr(c,t,1) = TFPCtr(c,t,1)*histPopCtr(c,ts)^LabourElast*KCtr(c,t,1)^(1-LabourElast);
    end
end

%impute missing observations
for t = NHistYear-NHistCtrYr+2:NHistYear,
    ts = t - NHistYear + NHistCtrYr;
    for c = 1:NCountry,
        if YCtr(c,t,1) == 0 & YCtr(c,t-1,1) > 0,
           TFPCtr(c,t,1) = (1+gTFP)*TFPCtr(c,t-1,1); 
           KCtr(c,t,1) = (1-Depreciation)*KCtr(c,t-1,1) + SavingsRate*YCtr(c,t-1,1);
           YCtr(c,t,1) = TFPCtr(c,t,1)*histPopCtr(c,ts)^LabourElast*KCtr(c,t,1)^(1-LabourElast);
        end
    end
end

YpCCtr = YCtr./PopCtr;

YpC2010Ctr = squeeze(YpCCtr(:,NHistYear,1));

%create scenarios
for s=1:NScen,
    TFPCtr(:,:,s) = TFPCtr(:,:,1);
    KCtr(:,:,s) = KCtr(:,:,1);
    YCtr(:,:,s) = YCtr(:,:,1);
end

for t=NHistYear+1:NYear
    ts = t-NHistYear;
    for s=1:NSRES
        for c = 1:NCountry
            TFPCtr(c,t,s) =  (1+SRESdInc(s,ts))*TFPCtr(c,t-1,s);
        end
    end
    for s=NSRES+1:NScen
        for c = 1:NCountry
            TFPCtr(c,t,s) =  (1+SSPdInc(s-NSRES,ts))*TFPCtr(c,t-1,s);
        end
    end
end

%energy
EnergyCtr = zeros(NCountry,NYear,NScen);

EnergyCtr(:,NHistYear-NHistCtrYr+1:NHistYear,1) = histEnergyCtr;

for i = 1:5,
    ti = 5-i;
    xYpC = YpCCtr(:,NHistYear-ti,1);
    xEnInt = EnergyCtr(:,NHistYear-ti,1)./YCtr(:,NHistYear-ti,1);
    xxYpC = xYpC(xEnInt > 0);
    xxEnInt = xEnInt(xEnInt > 0);
    xxYpC = log(xxYpC);
    xxEnInt = log(xxEnInt);
    xX = [ones(size(xxYpC,1),1) xxYpC];
    xb = inv(xX'*xX)*xX'*xxEnInt;
    for c=1:NCountry
        if xEnInt(c) ==0
            if EnergyCtr(c,NHistYear-ti-1,1) == 0
                xEnInt(c) = exp(xb(1) + xb(2)*log(xYpC(c)));
            else
                xEnInt(c) = EnergyCtr(c,NHistYear-ti-1,1)./YCtr(c,NHistYear-ti-1,1);
                xEnInt(c) = xEnInt(c)*(1+xb(2)*(YpCCtr(c,NHistYear-ti,1)/YpCCtr(c,NHistYear-ti-1,1)-1));
            end
            EnergyCtr(c,NHistYear-ti,1) = xEnInt(c)*YCtr(c,NHistYear-ti,1);    
        end
    end
end

%carbon dioxide
CO2Ctr = zeros(NCountry,NYear,NScen);

CO2Ctr(:,NHistYear-NHistCtrYr+1:NHistYear,1) = histCO2Ctr;

for i = 1:4,
    ti = 4-i;
    xYpC = YpCCtr(:,NHistYear-ti,1);
    xCO2Int = CO2Ctr(:,NHistYear-ti,1)./EnergyCtr(:,NHistYear-ti,1);
    xxYpC = xYpC(xCO2Int > 0);
    xxCO2Int = xEnInt(xCO2Int > 0);
    xxYpC = log(xxYpC);
    xxCO2Int = log(xxCO2Int);
    xX = [ones(size(xxYpC,1),1) xxYpC];
    xb = inv(xX'*xX)*xX'*xxCO2Int;
    for c=1:NCountry
        if xCO2Int(c) ==0
            xCO2Int(c) = exp(xb(1) + xb(2)*log(xYpC(c)));
            CO2Ctr(c,NHistYear-ti,1) = xCO2Int(c)*EnergyCtr(c,NHistYear-ti,1);    
        end
    end
end

clear x*

%scenarios
for s=1:NScen,
    EnergyCtr(:,:,s) = EnergyCtr(:,:,1);
    CO2Ctr(:,:,s) = CO2Ctr(:,:,1);
end

EnIntCtr = EnergyCtr./YCtr;
CO2IntCtr = CO2Ctr./EnergyCtr;

for t=NHistYear+1:NYear
    ts = t-NHistYear;
    for s=1:NSRES
        EnIntCtr(:,t,s)= (1+SRESdEnInt(s,ts))*EnIntCtr(:,t-1,s);
        CO2IntCtr(:,t,s)= (1+SRESdCO2Int(s,ts))*CO2IntCtr(:,t-1,s);
    end
    for s=NSRES+1:NScen
        EnIntCtr(:,t,s)= (1+SSPdEnInt(s-NSRES,ts))*EnIntCtr(:,t-1,s);
        CO2IntCtr(:,t,s)= (1+SSPdCO2Int(s-NSRES,ts))*CO2IntCtr(:,t-1,s);
    end
end