%% Loss Map vs Pwr Req for aging, Temp,
clear;
clc
cd('D:\Fuel_cell\Fuel_cell_model')
aging={'BOL', 'EOL', 'BOL_eTC', 'EOL_eTC'};
modelname='FCS.slx';uiopen(modelname,1)

% Define your DutyCycle and Pwrdes arrays
PDesMin=5;
PDesMax=150;
PDesSpcing=5;
DCMin=10;
DCSpacing=10;
DCMax=100;
Pwrdes=(PDesMin:PDesSpcing:PDesMax);
DutyCycle=(DCMin:DCSpacing:DCMax);

stabtime=40; %Time to run simulation for each condition for acheviong steady state
twdw=100;
% timearray=(0:stabtime:numel(Pwrdes)*stabtime')';
UseFCPSPwrReq=1;
Temp=25;%[48, 0, -30];%DegC
FCS_Age=aging{2};
FCSData=zeros(numel(Pwrdes),5);
% Initialize an empty table without headers
tableData = zeros(length(DutyCycle), length(Pwrdes));

% Initialize cell arrays for row and column names
rowNames = cell(length(DutyCycle), 1);
colNames = cell(length(Pwrdes), 1);

% Iterate through each combination of DutyCycle and Pwrdes
for j =1:numel(Pwrdes)
    Pdes = Pwrdes(j);
    FCSimPdes=[Pdes Pdes];
    FCSimPdesTime=[ 0   stabtime+twdw];
    for i =1:numel(DutyCycle)
        DC =DutyCycle(i);
        Preq = Pdes * 100 / DC;
        % Check if Preq is within the desired range (PDesMin to PDesMax)
        if PDesMin <= Preq && Preq <= PDesMax
            run('fuelCellOrg_v1_1_FcsSySp150d2_v2_standalone.m');
            run('Map4d_DCDC_BoostFCS_ConstEff.m');
            PwrReq     = [Preq Preq 0];
            PwrReqTime = [ 0    stabtime   stabtime+(DC/100)*twdw];
            FCPSTemp       = [Temp      Temp]+273;
            simstrtt=0;
            simendt=max(PwrReqTime);
            stpsz=0.2;
%             time=(simstrtt+.0001:stpsz:simendt)';
            SOC_overrideValue     = [50  50      50                50       ];%[50 50 ];            
            SOC_overrideEnable    = [1   0       1                 0        ];
            SOCInitTime           = [0   1   stabtime-1            stabtime ];
            Energy_overrideEnable = [1   0       1                 0        ];
            Enery_Time            = [0   1   stabtime-1            stabtime ];
            sim(modelname);
            % Store Pdes value in the table at the corresponding position
            FCResPReq(i, j) = Preq;
            FCResPDes(i,j)=  Pdes;
            FCResH2Cons(i,j)   =   FlwRate(end); %cumsum(mflow_H2.Data(stabtime:stabtime+twdw)*1000);%Mass of Hé consumed
            FCResChElLoss(i,j)= FCResH2Cons(i,j)*120  - Energy_FcsNetAct(end);
            FCSResEff = 100*Energy_FcsNetAct(end)/FCResH2Cons(i,j)*120;
            tempvar=getsampleusingtime(SOC,stabtime);tempvar=tempvar.Data;
            FCSResSOCStrt=tempvar(end);
           % FCSResSOCEnd
            
        else
            FCResPReq(i,j)=0;
        end
        % Assign DutyCycle values as row names
        rowNames{i} = num2str(DutyCycle(i));
        
    end
   
    
end

% Assign Pwrdes values as column names
for j = 1:length(Pwrdes)
    colNames{j} = num2str(Pwrdes(j));
end

% Create a table with string row and column names
%   dataTable = array2table(FCResPReq, 'RowNames', rowNames, 'VariableNames', colNames);
%dataTable1 = array2table(FCResH2Cons, 'RowNames', rowNames, 'VariableNames', colNames);
%dataTable2 = array2table(FCResChElLoss, 'RowNames', rowNames, 'VariableNames', colNames);
%dataTable3 = array2table(FCSResEff, 'RowNames', rowNames, 'VariableNames', colNames);
%dataTable4 = array2table(FCSResSOCStrt, 'RowNames', rowNames, 'VariableNames', colNames);


% Display the resulting table
%disp(dataTable);

