function EJLR_T = EJLR_Tensor(COSMOSS_Input,SpecType,N_Grid)
%% Debug
% % clear all
% SpecType = 'SFG';
% N_Grid   = 10;
% 
% Exp1D = zeros(2,3);
% Exp2D = zeros(2,5);
% Exp1D(2,:) = 90;
% Exp2D(2,:) = 90;
% COSMOSS_Input.Exp1D = Exp1D;
% COSMOSS_Input.Exp2D = Exp2D;

%% Decide Inputs

switch SpecType
    case 'SFG'
        Angle = COSMOSS_Input.Exp1D(1,:)./180.*pi;
        Polar = COSMOSS_Input.Exp1D(2,:)./180.*pi;
        
        [Jc,Jb,Ja] = ndgrid(1:3,1:3,1:3);
        Ja = Ja(:);
        Jb = Jb(:);
        Jc = Jc(:);  
        
        EJLR1 = EJLR_Grid(Polar(1),Angle(1),N_Grid,'Transimisive'); 
        EJLR2 = EJLR_Grid(Polar(2),Angle(2),N_Grid,'Transimisive');
        EJLR3 = EJLR_Grid(Polar(3),Angle(3),N_Grid,'Reflective');
        
        EJLR_T = EJLR3(:,Ja).*EJLR2(:,Jb).*EJLR1(:,Jc); % [NG * R2^3]
        
    case '2DSFG'
        Angle = COSMOSS_Input.Exp2D(1,:)./180.*pi;
        Polar = COSMOSS_Input.Exp2D(2,:)./180.*pi;
        
        [Je,Jd,Jc,Jb,Ja] = ndgrid(1:3,1:3,1:3,1:3,1:3);
        Ja = Ja(:);
        Jb = Jb(:);
        Jc = Jc(:);
        Jd = Jd(:);
        Je = Je(:);
        
        EJLR1 = EJLR_Grid(Polar(1),Angle(1),N_Grid,'Transimisive');
        EJLR2 = EJLR_Grid(Polar(2),Angle(2),N_Grid,'Transimisive');
        EJLR3 = EJLR_Grid(Polar(3),Angle(3),N_Grid,'Transimisive');
        EJLR4 = EJLR_Grid(Polar(4),Angle(4),N_Grid,'Transimisive');
        EJLR5 = EJLR_Grid(Polar(5),Angle(5),N_Grid,'Reflective');
        
        EJLR_T = EJLR5(:,Ja).*EJLR4(:,Jb).*EJLR3(:,Jc).*EJLR2(:,Jd).*EJLR1(:,Je); % [NG * R2^5]  
end

