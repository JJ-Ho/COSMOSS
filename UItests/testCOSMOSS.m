classdef testCOSMOSS < matlab.uitest.TestCase & matlab.mock.TestCase
    properties
        hApp
    end
    
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
            tc.hApp = app;
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
            tc.hApp = app;
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

        function test2DGrid(tc)
            import matlab.mock.actions.AssignOutputs
            FilesName = '131029_MBA_Reverse_TDV.txt';
            PathName = '/Users/jjho/Documents/MATLAB/COSMOSS/StructureFiles/G09/';
            
            [mockChooser,behavior] = tc.createMock(?FileChooser);
            when(behavior.chooseFile('*.*'),AssignOutputs(FilesName,PathName,1))
            
            
            app = COSMOSS(mockChooser);
            tc.hApp = app;
            tc.choose(app.ListBox_Model,3);
            tc.press(app.Button_SelectModel);
            tc.press(app.hModel.Button_LoadInput);
            tc.press(app.hModel.Button_Generate);
            tc.press(app.hModel.Button_Draw);
        end
        
        function test_BetaSheet(tc)
            app = COSMOSS;
            tc.hApp = app;
            
            tc.choose(app.ListBox_Model,4);
            tc.press(app.Button_SelectModel);
            tc.type(app.hModel.N_Residue,3);
            tc.type(app.hModel.N_Strand,2);
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
            
            load('expectedResults_BetaSheet','expectedResults')
            tc.verifyEqual(calculated,expectedResults)
        end
        
        function test_Comb2(tc)
            app = COSMOSS;
            tc.hApp = app;
            
            tc.choose(app.ListBox_Model,6);
            tc.press(app.Button_SelectModel);
            tc.press(app.hModel.Button_Model1);
            tc.press(app.hModel.Button_Model2);
            tc.press(app.hModel.hModel1.Button_Generate);
            tc.press(app.hModel.hModel2.Button_Generate);
            tc.type(app.hModel.Trans_Z,5);
            tc.press(app.hModel.Button_Combine);
            tc.press(app.Button_FTIR);
            tc.press(app.Button_SFG);
            tc.press(app.Button_2DIR);
            tc.press(app.Button_2DSFG);
            
            calculated.Structure  = app.Structure;
            calculated.Data_FTIR  = app.Data_FTIR;
            calculated.Data_SFG   = app.Data_SFG;
            calculated.Data_2DIR  = app.Data_2DIR;
            calculated.Data_2DSFG = app.Data_2DSFG;
            
            load('expectedResults_Comb2','expectedResults')
            tc.verifyEqual(calculated,expectedResults)
        end
        
    end
    
    methods (TestMethodTeardown)
        function closeCOSMOSS(tc)
            delete(tc.hApp.hModel)
            delete(tc.hApp)
            close all
        end
    end
end