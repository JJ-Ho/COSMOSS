%% TCO
GUI_Inputs = struct;
GUI_Inputs.Theta_D1 = -30;
GUI_Inputs.Theta_D2 = 30;
GUI_Inputs.FreqRange = 1650:1800;
GUI_Inputs.LocFreqType = 'Symbolic';
GUI_Inputs.CouplingType = 'Symbolic';
SD_TCO = ConstructTCO(GUI_Inputs);

SD_TCO_H2ex = H_handler(SD_TCO,GUI_Inputs,'TwoEx');

%% PDB_AmideI
GUI_Inputs = struct;
GUI_Inputs.Preprocessed = false;
GUI_Inputs.LocFreqType = 'Symbolic';
GUI_Inputs.CouplingType = 'Symbolic';

app = struct;
app.Parent = [];

S_PDB = Load_PDB(app,GUI_Inputs);
S_PDB_AmideI = ConstructAmideIPDB(S_PDB,GUI_Inputs);