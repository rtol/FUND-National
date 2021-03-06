%readdata
%The Climate Framework for Uncertainty, Negotiation and Distribution,
%version 4.0-matlab-national
%
%This script is part of FUND 4.0 MN
%It reads historic data
%
%Richard Tol, 2 September 2014
%This code is protected by the MIT License

historicCO2emit = csvread('histCO2emit.csv');
historicLUemit = csvread('histLUemit.csv');
histCO2conc = csvread('histCO2conc.csv');

historicCH4emit = csvread('histCH4emit.csv');
histCH4conc = csvread('histCH4conc.csv');
historicN2Oemit = csvread('histN2Oemit.csv');
histN2Oconc = csvread('histN2Oconc.csv');
historicSF6emit = csvread('histSF6emit.csv');
histSF6conc = csvread('histSF6conc.csv');
histCFC11emit = csvread('histCFC11emit.csv');
histCFC11conc = csvread('histCFC11conc.csv');
histCFC12emit = csvread('histCFC12emit.csv');
histCFC12conc = csvread('histCFC12conc.csv');
histO3radforc = csvread('histO3radforc.csv');

historicSemit = csvread('histSemit.csv');

historicTemp = csvread('histTemp.csv');
NatTemp = csvread('NationalTemperature.csv');

observedImpact = csvread('TotalImpact.csv');
observedNatImp = csvread('NationalImpact.csv');

histPopulation = csvread('histPopulation.csv');   %global population, 1750-2010
histGDP = csvread('histGDP.csv');
histYpC = histGDP./histPopulation;
histEnergy = csvread('histEnergy.csv');
histEnInt = histEnergy./histGDP;
histCO2Int = historicCO2emit./histEnergy;

histPopCtr = csvread('histPopCountry.csv');       %national population, 1960-2010
ctr = size(histPopCtr);
NHistCtrYr = ctr(2);
NCountry = ctr(1);
StartCtrYr = EndHistYear-NHistCtrYr+1;
clear ctr

histGDPCtr = csvread('histGDPCountry.csv');       %national GDP, 1960-2010, note many missing observations
histEnergyCtr = csvread('histEnergyCountry.csv');%national primary energy use, 1960-2010, note many missing observations
histEnergyCtr = histEnergyCtr/1000;
histCO2Ctr = csvread('histCO2Country.csv');        %national carbon dioxide emissions from fossil fuel combustion, 1960-2010, note many missing observations
histCO2Ctr = histCO2Ctr*12/44/1000;
