classdef testFTIR < matlab.uitest.TestCase
    
    methods (Test)
        function test_FTIR_TCO(tc)
            app = COSMOSS;
            tc.addTeardown(@closeCOSMOSS,app);
            
            tc.choose(app.ListBox_Model,1);
            tc.press(app.Button_SelectModel);
            tc.press(app.hModel.Button_Generate);
            tc.press(app.Button_FTIR);
            
            Response1D = app.Data_FTIR.Response1D;
            
            load('expectedResults','expectedResponse1D')
            tc.verifyEqual(Response1D,expectedResponse1D);
        end
        
%         function test_FTIR_PDB(tc)
%             app = COSMOSS;
%             tc.addTeardown(@closeCOSMOSS,app);
%             
%             tc.choose(app.ListBox_Model,2);
%             tc.press(app.Button_SelectModel);
%             tc.press(app.hModel.Button_LoadPDB);
%             tc.press(app.hModel.Button_Generate);
%             tc.press(app.Button_FTIR);
%             
%             Response1D = app.Data_FTIR.Response1D;
%             
%             load('expectedResults','expectedResponse1D')
%             tc.verifyEqual(Response1D,expectedResponse1D);
%         end

    end
    
end