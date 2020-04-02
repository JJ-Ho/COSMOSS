classdef testCOSMOSS < matlab.uitest.TestCase & matlab.mock.TestCase
    methods(TestClassSetup)
        function addPath(tc)
            % get the upper layer path before the tests begin
            [COSMOSSPath,~,~] = fileparts(pwd);
            addpath(COSMOSSPath)
        end
    end
    
    methods (Test)
        function test_TCO(tc)
            app = COSMOSS;
            tc.addTeardown(@closeCOSMOSS,app);
            
            tc.choose(app.ListBox_Model,1);
            tc.press(app.Button_SelectModel);
            tc.press(app.hModel.Button_Generate);
            tc.press(app.Button_FTIR);
            tc.press(app.Button_SFG);
            tc.press(app.Button_2DIR);
            tc.press(app.Button_2DSFG);
            
            calculated.Structure  = app.Structure;
            calculated.Data_FTIR  = app.Data_FTIR;
            calculated.Data_SFG   = app.Data_SFG;
            calculated.Data_2DIR  = app.Data_2DIR;
            calculated.Data_2DSFG = app.Data_2DSFG;
            
            load('expectedResults_TCO','expectedResults')
            tc.verifyEqual(calculated,expectedResults)
        end
        
        function test_PDB(tc)
            import matlab.mock.actions.AssignOutputs
            FilesName = '5-glycine-helix-z.pdb';
            PathName = '/Users/jjho/Documents/MATLAB/COSMOSS/StructureFiles/PDB/Helix/';
            
            [mockChooser,behavior] = tc.createMock(?FileChooser);
            when(behavior.chooseFile('*.*'),AssignOutputs(FilesName,PathName,1))
            
            app = COSMOSS(mockChooser);
            tc.addTeardown(@closeCOSMOSS,app);
            
            tc.choose(app.ListBox_Model,2);
            tc.press(app.Button_SelectModel);
            tc.press(app.hModel.Button_LoadPDB);
            tc.press(app.hModel.Button_Generate);
            tc.press(app.Button_FTIR);
            tc.press(app.Button_SFG);
            tc.press(app.Button_2DIR);
            tc.press(app.Button_2DSFG);
            
            calculated.Structure  = app.Structure;
            calculated.Data_FTIR  = app.Data_FTIR;
            calculated.Data_SFG   = app.Data_SFG;
            calculated.Data_2DIR  = app.Data_2DIR;
            calculated.Data_2DSFG = app.Data_2DSFG;
            
            load('expectedResults_PDB_5GlyHelix','expectedResults')
            tc.verifyEqual(calculated,expectedResults)
        end

    end
    
end