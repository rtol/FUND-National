%CalibImpact
%The Climate Framework for Uncertainty, Negotiation and Distribution,
%version 4.1-matlab-global
%
%This script is part of FUND 4.1 MG
%It calibrates five alternative aggregrate damage functions
%
%Richard Tol, 14 May 2018
%This code is protected by the MIT License

vtemp = observedImpact(:,1);
vimp = observedImpact(:,2);
vN = length(vtemp);
vtemp = vtemp(1:vN-1);
vimp = vimp(1:vN-1);

imppar = zeros(NImpact,4);

[imppar(1,1) imppar(1,2) imppar(1,4)] = fittol(vtemp,vimp);
[imppar(2,1) imppar(2,2) imppar(2,4)] = fitweitzman(vtemp,vimp);
[imppar(3,1) imppar(3,2) imppar(3,4)] = fitweitzman(vtemp,vimp);
[imppar(4,1) imppar(4,4)] = fitnordhaus(vtemp,vimp);
[imppar(5,1) imppar(5,4)] = fithope(vtemp,vimp);
[imppar(6,1) imppar(6,4)] = fitploeg(vtemp,vimp);
[imppar(7,1) imppar(7,4)] = fitgolosov(vtemp,vimp);
[imppar(8,1) imppar(8,2) imppar(8,3) imppar(8,4)] = fittol2(vtemp,vimp);

%national impact
ctrimppar = zeros(NImpact,NCountry,2);
xnatimp = zeros(NImpact,NCountry);

for c = 1:NCountry,
    ctrimppar(1,c,1) = imppar(1,1)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    ctrimppar(1,c,2) = imppar(1,2)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    xnatimp(1,c) = ctrimppar(1,c,1)*2.5 + ctrimppar(1,c,2)*2.5^2;
    ctrimppar(2,c,1) = imppar(2,1)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    ctrimppar(2,c,2) = imppar(2,2)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    xnatimp(2,c) = ctrimppar(2,c,1)*2.5^2 + ctrimppar(2,c,2)*2.5^6;
    ctrimppar(3,c,1) = imppar(3,1)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    ctrimppar(3,c,2) = imppar(3,2)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    xnatimp(3,c) = ctrimppar(3,c,1)*2.5^2 + ctrimppar(3,c,2)*2.5^7;
    ctrimppar(4,c,1) = imppar(4,1)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    xnatimp(4,c) = ctrimppar(4,c,1)*2.5^2;
    ctrimppar(5,c,1) = imppar(5,1)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    xnatimp(5,c) = ctrimppar(5,c,1)*2.5;
    ctrimppar(6,c,1) = imppar(6,1)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    xnatimp(6,c) = ctrimppar(6,c,1)*(exp(2.5)-1);
    ctrimppar(7,c,1) = imppar(7,1)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    xnatimp(7,c) = ctrimppar(7,c,1)*(exp(exp(2.5))-exp(1));
    ctrimppar(8,c,1) = imppar(8,1)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    ctrimppar(8,c,2) = imppar(8,2)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    ctrimppar(8,c,3) = imppar(8,3)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    xnatimp(8,c) = evaltol2(ctrimppar(8,c,1:3),2.5);
end

%rescale to match global impact for 2010 income and 2.5K warming
xtotimp = sum(xnatimp(1,:)'.*YCtr(:,NHistYear,1));
xglbimp = (imppar(1,1)*2.5 + imppar(1,2)*(2.5^2))*Y(NHistYear,1);
ctrimppar(1,:,1) = ctrimppar(1,:,1)*xtotimp/xglbimp;
ctrimppar(1,:,2) = ctrimppar(1,:,2)*xtotimp/xglbimp;

xtotimp = sum(xnatimp(2,:)'.*YCtr(:,NHistYear,1));
xglbimp = (imppar(2,1)*(2.5^2) + imppar(2,2)*(2.5^6))*Y(NHistYear,1);
ctrimppar(2,:,1) = ctrimppar(2,:,1)*xtotimp/xglbimp;
ctrimppar(2,:,2) = ctrimppar(2,:,2)*xtotimp/xglbimp;

xtotimp = sum(xnatimp(3,:)'.*YCtr(:,NHistYear,1));
xglbimp = (imppar(3,1)*(2.5^2) + imppar(2,2)*(2.5^7))*Y(NHistYear,1);
ctrimppar(3,:,1) = ctrimppar(3,:,1)*xtotimp/xglbimp;
ctrimppar(3,:,2) = ctrimppar(3,:,2)*xtotimp/xglbimp;

xtotimp = sum(xnatimp(4,:)'.*YCtr(:,NHistYear,1));
xglbimp =imppar(4,1)*(2.5^2)*Y(NHistYear,1);
ctrimppar(4,:,1) = ctrimppar(4,:,1)*xtotimp/xglbimp;

xtotimp = sum(xnatimp(5,:)'.*YCtr(:,NHistYear,1));
xglbimp =imppar(5,1)*2.5*Y(NHistYear,1);
ctrimppar(5,:,1) = ctrimppar(5,:,1)*xtotimp/xglbimp;

xtotimp = sum(xnatimp(6,:)'.*YCtr(:,NHistYear,1));
xglbimp =imppar(6,1)*(exp(2.5)-1)*Y(NHistYear,1);
ctrimppar(6,:,1) = ctrimppar(6,:,1)*xtotimp/xglbimp;

xtotimp = sum(xnatimp(7,:)'.*YCtr(:,NHistYear,1));
xglbimp =imppar(7,1)*(exp(exp(2.5))-exp(1))*Y(NHistYear,1);
ctrimppar(7,:,1) = ctrimppar(7,:,1)*xtotimp/xglbimp;

xtotimp = sum(xnatimp(8,:)'.*YCtr(:,NHistYear,1));
xglbimp =evaltol2(imppar(8,1:3),2.5)*Y(NHistYear,1);
ctrimppar(8,:,1) = ctrimppar(8,:,1)*xtotimp/xglbimp;
ctrimppar(8,:,2) = ctrimppar(8,:,2)*xtotimp/xglbimp;
ctrimppar(8,:,3) = ctrimppar(8,:,3)*xtotimp/xglbimp;

clear x*