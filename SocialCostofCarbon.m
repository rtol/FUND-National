%SocialCostofCarbon
%The Climate Framework for Uncertainty, Negotiation and Distribution,
%version 4.0-matlab-national
%
%This script is part of FUND 4.0 MN
%It computes the social cost of carbon
%
%Richard Tol, 24 September 2014
%This code is protected by the MIT License

dimpact = impactd - impact;
dimpctr = impctrd - impctr;

for i=1:NImpact,
    dimpabs(i,:,:) = squeeze(dimpact(i,:,:)).*Y;
    dimpabsctr(i,:,:,:) = squeeze(dimpctr(i,:,:,:)).*YCtr;
end

dimpabsctr(isnan(dimpabsctr)) = 0;

gy = zeros(NYear,NScen);
gyc = zeros(NCountry,NYear,NScen);
for t=2:NYear,
    for s=1:NScen,
        gy(t,s) = (YpC(t,s)-YpC(t-1,s))/YpC(t-1,s);
        for c=1:NCountry,
            gyc(c,t,s) = (YpCCtr(c,t,s)-YpCCtr(c,t-1,s))/YpCCtr(c,t-1,s);
        end
    end
end

gp = zeros(NYear,NScen);
gpc = zeros(NCountry,NYear,NScen);
for t=2:NYear,
    for s=1:NScen,
        gp(t,s) = (Population(t,s)-Population(t-1,s))/Population(t-1,s);
        for c=1:NCountry,
            gpc(c,t,s) = (PopCtr(c,t,s)-PopCtr(c,t-1,s))/PopCtr(c,t-1,s);
        end
    end
end

NDR = 6;
PRTP = [0.001 0.010 0.020 0.030 0.040 0.050];
NRA = 5;
RA = [0.5 1.0 1.5 2.0 2.5];
NCDR = 5;
DR = [0.0700 0.0487 0.05688 0.0300 0.0250];

df = zeros(NYear,NScen,NDR,NRA);
df2 = zeros(NYear,NCDR);
dfc = zeros(NCountry,NYear,NScen,NDR,NRA);
for s=1:NScen,
    for p=1:NDR,
        for r=1:NRA,
            df(SCCYear,s,p,r) = 1;
            for c=1:NCountry
                dfc(c,SCCYear,s,p,r) = 1;
            end
        end
    end
end
for r=1:NCDR,
    df2(SCCYear,r) = 1;
end

for t=SCCYear+1:NYear,
    for s=1:NScen,
        for p=1:NDR,
            for r=1:NRA
                df(t,s,p,r) = df(t-1,s,p,r)/(1+PRTP(p)+gp(t-1,s)+RA(r)*gy(t-1,s));
                for c=1:NCountry
                    dfc(c,t,s,p,r) = dfc(c,t-1,s,p,r)/(1+PRTP(p)+gpc(c,t-1,s)+RA(r)*gyc(c,t-1,s));
                end
            end
        end
    end
    for r=1:NCDR
        df2(t,r) = df2(t-1,r)/(1+DR(r));
    end
end

SCC = zeros(NImpact, NScen, NDR, NRA);
SCC2 = zeros(NImpact, NScen, NCDR);
SCCc = zeros(NImpact, NCountry, NScen, NDR, NRA);
SCCc2 = zeros(NImpact, NCountry, NScen, NCDR);

for i=1:NImpact,
    for s=1:NScen,
        for p=1:NDR,
            for r=1:NRA,
                SCC(i,s,p,r) = squeeze(df(:,s,p,r))'*dimpabs(i,:,s)';
                for c= 1:NCountry
                    SCCc(i,c,s,p,r) = dfc(c,:,s,p,r)*squeeze(dimpabsctr(i,c,:,s));
                end
            end
        end
        for p=1:NCDR,
            SCC2(i,s,p) = df2(:,p)'*dimpabs(i,:,s)';
            for c=1:NCountry
                SCCc2(i,c,s,p) = df2(:,p)'*squeeze(dimpabsctr(i,c,:,s));
            end
        end
    end
end

SCC = -0.01*SCC/1000000;
SCCc = -0.01*SCCc/1000000;
SCC2 = -0.01*SCC2/1000000;
SCCc2 = -0.01*SCCc2/1000000;

s = PrintTable;
s.addRow('Time pref \ Risk aversion', RA(1), RA(2), RA(3), RA(4), RA(5));
s.addRow(num2str(PRTP(1),5),num2str(SCC(1,1,1,1),7),num2str(SCC(1,1,1,2),7),num2str(SCC(1,1,1,3),7),num2str(SCC(1,1,1,4),7),num2str(SCC(1,1,1,5),7));
s.addRow(num2str(PRTP(2),5),num2str(SCC(1,1,2,1),7),num2str(SCC(1,1,2,2),7),num2str(SCC(1,1,2,3),7),num2str(SCC(1,1,2,4),7),num2str(SCC(1,1,2,5),7));
s.addRow(num2str(PRTP(3),5),num2str(SCC(1,1,3,1),7),num2str(SCC(1,1,3,2),7),num2str(SCC(1,1,3,3),7),num2str(SCC(1,1,3,4),7),num2str(SCC(1,1,3,5),7));
s.addRow(num2str(PRTP(4),5),num2str(SCC(1,1,4,1),7),num2str(SCC(1,1,4,2),7),num2str(SCC(1,1,4,3),7),num2str(SCC(1,1,4,4),7),num2str(SCC(1,1,4,5),7));
s.addRow(num2str(PRTP(5),5),num2str(SCC(1,1,5,1),7),num2str(SCC(1,1,5,2),7),num2str(SCC(1,1,5,3),7),num2str(SCC(1,1,5,4),7),num2str(SCC(1,1,5,5),7));
s.addRow(num2str(PRTP(6),5),num2str(SCC(1,1,6,1),7),num2str(SCC(1,1,6,2),7),num2str(SCC(1,1,6,3),7),num2str(SCC(1,1,6,4),7),num2str(SCC(1,1,6,5),7));
disp('Social cost of carbon ($/tC)')
disp('alternative rates of time preference (rows) and risk aversion (columns)')
str = ['model = Tol, scenario = SRES A1 ' num2str(SCCYear+StartYear)];
disp(str)
s.display

line = sprintf('\n');
disp(line)

t = PrintTable;
t.addRow('Model \ Scenario', 'SRES A1', 'SRES A2', 'SRES B1', 'SRES B2', 'SSP1', 'SSP2', 'SSP3', 'SSP4', 'SSP5');
t.addRow('Tol (parabola)',num2str(SCC(1,1,4,2),7),num2str(SCC(1,2,4,2),7),num2str(SCC(1,3,4,2),7),num2str(SCC(1,4,4,2),7),num2str(SCC(1,5,4,2),7),num2str(SCC(1,6,4,2),7),num2str(SCC(1,7,4,2),7),num2str(SCC(1,8,4,2),7),num2str(SCC(1,9,4,2),7));
t.addRow('Weitzman (6)',num2str(SCC(2,1,4,2),7),num2str(SCC(2,2,4,2),7),num2str(SCC(2,3,4,2),7),num2str(SCC(2,4,4,2),7),num2str(SCC(2,5,4,2),7),num2str(SCC(2,6,4,2),7),num2str(SCC(2,7,4,2),7),num2str(SCC(2,8,4,2),7),num2str(SCC(2,9,4,2),7));
t.addRow('Weitzman (7)',num2str(SCC(3,1,4,2),7),num2str(SCC(3,2,4,2),7),num2str(SCC(3,3,4,2),7),num2str(SCC(3,4,4,2),7),num2str(SCC(3,5,4,2),7),num2str(SCC(3,6,4,2),7),num2str(SCC(3,7,4,2),7),num2str(SCC(3,8,4,2),7),num2str(SCC(3,9,4,2),7));
t.addRow('Nordhaus',num2str(SCC(4,1,4,2),7),num2str(SCC(4,2,4,2),7),num2str(SCC(4,3,4,2),7),num2str(SCC(4,4,4,2),7),num2str(SCC(4,5,4,2),7),num2str(SCC(4,6,4,2),7),num2str(SCC(4,7,4,2),7),num2str(SCC(4,8,4,2),7),num2str(SCC(4,9,4,2),7));
t.addRow('Hope',num2str(SCC(5,1,4,2),7),num2str(SCC(5,2,4,2),7),num2str(SCC(5,3,4,2),7),num2str(SCC(5,4,4,2),7),num2str(SCC(5,5,4,2),7),num2str(SCC(5,6,4,2),7),num2str(SCC(5,7,4,2),7),num2str(SCC(5,8,4,2),7),num2str(SCC(5,9,4,2),7));
t.addRow('Ploeg',num2str(SCC(6,1,4,2),7),num2str(SCC(6,2,4,2),7),num2str(SCC(6,3,4,2),7),num2str(SCC(6,4,4,2),7),num2str(SCC(6,5,4,2),7),num2str(SCC(6,6,4,2),7),num2str(SCC(6,7,4,2),7),num2str(SCC(6,8,4,2),7),num2str(SCC(6,9,4,2),7));
t.addRow('Golosov',num2str(SCC(7,1,4,2),7),num2str(SCC(7,2,4,2),7),num2str(SCC(7,3,4,2),7),num2str(SCC(7,4,4,2),7),num2str(SCC(7,5,4,2),7),num2str(SCC(7,6,4,2),7),num2str(SCC(7,7,4,2),7),num2str(SCC(7,8,4,2),7),num2str(SCC(7,9,4,2),7));
t.addRow('Tol2 (piecewise linear)',num2str(SCC(8,1,4,2),7),num2str(SCC(8,2,4,2),7),num2str(SCC(8,3,4,2),7),num2str(SCC(8,4,4,2),7),num2str(SCC(8,5,4,2),7),num2str(SCC(8,6,4,2),7),num2str(SCC(8,7,4,2),7),num2str(SCC(8,8,4,2),7),num2str(SCC(8,9,4,2),7));
disp('Social cost of carbon ($/tC)')
disp('alternative impact models (rows) and scenarios (columns)')
str = ['pure rate of time preference = 0.03; rate of risk aversion = 1 ' num2str(SCCYear+StartYear)];
disp(str)
t.display