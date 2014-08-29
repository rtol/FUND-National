%Initialize
%The Climate Framework for Uncertainty, Negotiation and Distribution,
%version 4.0-matlab-national
%
%This script is part of FUND 4.0 MN
%It initializes variables and sets parameters
%
%Richard Tol, 29 August 2014
%This code is protected by the MIT License

clear all

StartYear = 1750;
StartCtrYr = 1960;
EndHistYear = 2010;
NHistYear = EndHistYear - StartYear + 1;
NHistCtrYr = EndHistYear - StartCtrYr + 1;
EndYear = 2100;
NYear = EndYear - StartYear + 1;
Year = zeros(NYear,1);
Year(1) = StartYear;
for t=2:NYear,
    Year(t) = Year(t-1) + 1;
end

NScen = 4;
NCountry = 209;