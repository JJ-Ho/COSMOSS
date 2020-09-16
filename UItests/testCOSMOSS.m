classdef testCOSMOSS < matlab.uitest.TestCase & matlab.mock.TestCase
    properties
        hApp
        RelTol = 1e-12
        AbsTol = 1e-8
        COSMOSSPath
        expectedPath
        saveResult = 0
    end
    
    methods(TestClassSetup)
        function addPath(tc)
            % get the upper layer path before the tests begin
            [Path,~,~] = fileparts(pwd);
            addpath(Path)
            
            tc.COSMOSSPath = Path;
            tc.expectedPath = './expected/';
        end
    end

    methods (Test)
        
        function test_NCO(tc) 
            caseName = 'NCO';
            app = COSMOSS;
            tc.hApp = app;
            tc.choose(app.ListBox_Model,1);
            tc.press(app.Button_SelectModel);
            tc.press(app.hModel.Button_Generate);
            tc.choose(app.DropDown_RotAvg,'Isotropic')
            tc.press(app.Button_FTIR);
            tc.choose(app.DropDown_RotAvg,'C1')
            tc.press(app.Button_SFG);
            tc.choose(app.DropDown_RotAvg,'Isotropic')
            tc.press(app.Button_2DIR);
            tc.choose(app.DropDown_RotAvg,'C1')
            tc.press(app.Button_2DSFG);
            
            calculated.Structure  = app.Structure;
            calculated.Data_FTIR  = app.Data_FTIR;
            calculated.Data_SFG   = app.Data_SFG;
            calculated.Data_2DIR  = app.Data_2DIR;
            calculated.Data_2DSFG = app.Data_2DSFG;
            
            if tc.saveResult
                disp('Saving result...')
                expectedResults = calculated;
                save(['expectedResults_',caseName],'expectedResults')
            else
                load([tc.expectedPath,'expectedResults_',caseName],'expectedResults')
                tc.verifyEqual(calculated,expectedResults,'AbsTol',tc.AbsTol)
            end
        end
        
        function test_PDB(tc)
            caseName = 'PDB_5GlyHelix';
            import matlab.mock.actions.AssignOutputs
            FilesName = '5-glycine-helix-z.pdb';
            PathName = [tc.COSMOSSPath,'/StructureFiles/PDB/Helix/'];
            
            [mockChooser,behavior] = tc.createMock(?FileChooser);
            when(behavior.chooseFile('*.*'),AssignOutputs(FilesName,PathName,1))
            
            app = COSMOSS(mockChooser);
            tc.hApp = app;
            tc.choose(app.ListBox_Model,2);
            tc.press(app.Button_SelectModel);
            tc.press(app.hModel.Button_LoadPDB);
            tc.press(app.hModel.Button_Generate);
            tc.choose(app.DropDown_RotAvg,'Isotropic')
            tc.press(app.Button_FTIR);
            tc.choose(app.DropDown_RotAvg,'C1')
            tc.press(app.Button_SFG);
            tc.choose(app.DropDown_RotAvg,'Isotropic')
            tc.press(app.Button_2DIR);
            tc.choose(app.DropDown_RotAvg,'C1')
            tc.press(app.Button_2DSFG);
            
            calculated.Structure  = app.Structure;
            calculated.Data_FTIR  = app.Data_FTIR;
            calculated.Data_SFG   = app.Data_SFG;
            calculated.Data_2DIR  = app.Data_2DIR;
            calculated.Data_2DSFG = app.Data_2DSFG;
            
            if tc.saveResult
                disp('Saving result...')
                expectedResults = calculated;
                save(['expectedResults_',caseName],'expectedResults')
            else
                load([tc.expectedPath,'expectedResults_',caseName],'expectedResults')
                tc.verifyEqual(calculated,expectedResults,'AbsTol',tc.AbsTol)
            end
        end

        function test2DGrid(tc)
            caseName = '2DGrid';
            import matlab.mock.actions.AssignOutputs
            FilesName = '131029_MBA_Reverse_TDV.txt';
            PathName = [tc.COSMOSSPath,'/StructureFiles/G09/'];
            
            [mockChooser,behavior] = tc.createMock(?FileChooser);
            when(behavior.chooseFile('*.*'),AssignOutputs(FilesName,PathName,1))
            
            
            app = COSMOSS(mockChooser);
            tc.hApp = app;
            tc.choose(app.ListBox_Model,3);
            tc.press(app.Button_SelectModel);
            tc.press(app.hModel.Button_LoadInput);
            tc.press(app.hModel.Button_Generate);
            tc.choose(app.DropDown_RotAvg,'Isotropic')
            tc.press(app.Button_FTIR);
            tc.choose(app.DropDown_RotAvg,'C1')
            tc.press(app.Button_SFG);
            tc.choose(app.DropDown_RotAvg,'Isotropic')
            tc.press(app.Button_2DIR);
            tc.choose(app.DropDown_RotAvg,'C1')
            tc.press(app.Button_2DSFG);

            calculated.Structure  = app.Structure;
            calculated.Data_FTIR  = app.Data_FTIR;
            calculated.Data_SFG   = app.Data_SFG;
            calculated.Data_2DIR  = app.Data_2DIR;
            calculated.Data_2DSFG = app.Data_2DSFG;
            
            if tc.saveResult
                disp('Saving result...')
                expectedResults = calculated;
                save(['expectedResults_',caseName],'expectedResults')
            else
                load([tc.expectedPath,'expectedResults_',caseName],'expectedResults')
                tc.verifyEqual(calculated,expectedResults,'AbsTol',tc.AbsTol)
            end
        end
        
        function test_BetaSheet(tc)
            caseName = 'BetaSheet';
            app = COSMOSS;
            tc.hApp = app;
            
            tc.choose(app.ListBox_Model,4);
            tc.press(app.Button_SelectModel);
            tc.type(app.hModel.N_Residue,2);
            tc.type(app.hModel.N_Strand,2);
            tc.press(app.hModel.Button_Generate);
            tc.choose(app.DropDown_RotAvg,'Isotropic')
            tc.press(app.Button_FTIR);
            tc.choose(app.DropDown_RotAvg,'C1')
            tc.press(app.Button_SFG);
            tc.choose(app.DropDown_RotAvg,'Isotropic')
            tc.press(app.Button_2DIR);
            tc.choose(app.DropDown_RotAvg,'C1')
            tc.press(app.Button_2DSFG);
            
            calculated.Structure  = app.Structure;
            calculated.Data_FTIR  = app.Data_FTIR;
            calculated.Data_SFG   = app.Data_SFG;
            calculated.Data_2DIR  = app.Data_2DIR;
            calculated.Data_2DSFG = app.Data_2DSFG;
            
            if tc.saveResult
                disp('Saving result...')
                expectedResults = calculated;
                save(['expectedResults_',caseName],'expectedResults')
            else
                load([tc.expectedPath,'expectedResults_',caseName],'expectedResults')
                tc.verifyEqual(calculated,expectedResults,'AbsTol',tc.AbsTol)
            end
        end
        
        function test_Comb2(tc)
            caseName = 'Comb2';
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
            tc.choose(app.DropDown_RotAvg,'Isotropic')
            tc.press(app.Button_FTIR);
            tc.choose(app.DropDown_RotAvg,'C1')
            tc.press(app.Button_SFG);
            tc.choose(app.DropDown_RotAvg,'Isotropic')
            tc.press(app.Button_2DIR);
            tc.choose(app.DropDown_RotAvg,'C1')
            tc.press(app.Button_2DSFG);
            
            calculated.Structure  = app.Structure;
            calculated.Data_FTIR  = app.Data_FTIR;
            calculated.Data_SFG   = app.Data_SFG;
            calculated.Data_2DIR  = app.Data_2DIR;
            calculated.Data_2DSFG = app.Data_2DSFG;
            
            if tc.saveResult
                disp('Saving result...')
                expectedResults = calculated;
                save(['expectedResults_',caseName],'expectedResults')
            else
                load([tc.expectedPath,'expectedResults_',caseName],'expectedResults')
                tc.verifyEqual(calculated,expectedResults,'AbsTol',tc.AbsTol)
            end
        end
        
    end
    
    methods (TestMethodTeardown)
        function closeCOSMOSS(tc)
            UUID = tc.hApp.UUID;
            delete(tc.hApp)
            
            % find childrens
            hUIFigures = findall(0, 'HandleVisibility', 'off');
            for j = 1:length(hUIFigures)
                if isequal(hUIFigures(j).Tag,UUID)
                    delete(hUIFigures(j))
                end
            end
            close all
        end
    end
end