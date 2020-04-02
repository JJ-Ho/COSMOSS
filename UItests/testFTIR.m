classdef testFTIR < matlab.uitest.TestCase & matlab.mock.TestCase
    methods(TestClassSetup)
        function addPath(tc)
            % get the upper layer path before the tests begin
            [COSMOSSPath,~,~] = fileparts(pwd);
            addpath(COSMOSSPath)
        end
    end
    
    methods (Test)
        function test_FTIR_TCO(tc)
            app = COSMOSS;
            tc.addTeardown(@closeCOSMOSS,app);
            
            tc.choose(app.ListBox_Model,1);
            tc.press(app.Button_SelectModel);
            tc.press(app.hModel.Button_Generate);
            tc.press(app.Button_FTIR);
            
            Response1D = app.Data_FTIR.Response1D;
            
            load('expectedResults','FTIR')
            tc.verifyEqual(Response1D,FTIR.TCO);
        end
        
        function test_FTIR_PDB(tc)
            import matlab.mock.actions.AssignOutputs
            FilesName = '2RRI.pdb';
            PathName = '/Users/jjho/Documents/MATLAB/COSMOSS/StructureFiles/PDB/';
            
            [mockChooser,behavior] = tc.createMock(?FileChooser);
            when(behavior.chooseFile('*.*'),AssignOutputs(FilesName,PathName,1))
            
            app = COSMOSS(mockChooser);
            tc.addTeardown(@closeCOSMOSS,app);
            
            tc.choose(app.ListBox_Model,2);
            tc.press(app.Button_SelectModel);
            tc.press(app.hModel.Button_LoadPDB);
            tc.press(app.hModel.Button_Generate);
            tc.press(app.Button_FTIR);
            
            Response1D = app.Data_FTIR.Response1D;
            
            load('expectedResults','FTIR')
            tc.verifyEqual(Response1D,FTIR.PDB_2RRI);
        end

    end
    
end