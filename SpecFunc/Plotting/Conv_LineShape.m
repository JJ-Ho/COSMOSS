function [ConvL,AvalibleL] = Conv_LineShape(Dimension,LineShape,FreqRange,LineWidth)
AvalibleL = {'Lorentzian','Gaussian','No Conv','Spy'};

N_Grid = length(FreqRange);
spec_array = (1:N_Grid) - ceil(N_Grid/2);

switch Dimension
    case 'List'
        ConvL = '';
        % this case is use to export "AvalibleL"
    case 1
        switch LineShape 
            case 'Lorentzian' 
                LineWidth = LineWidth/2;
                ConvL = spec_array./((spec_array.^2)+(LineWidth^2)) + 1i*LineWidth./(spec_array.^2+LineWidth^2);
            case 'Gaussian'
                ConvL = 1i*exp(-(spec_array.^2)./(LineWidth^2));
            case 'KK' % K-K use K-K relation to generate Re part
                disp('not support KK lineshape in 2D yet...')
%                 LineWidth = LineWidth/2;
%                 Im = (1/pi)*(LineWidth)./(spec_array.^2+(LineWidth)^2);
%                 Re = kkimbook2(FreqRange,Im,1);
%                 ConvL = Im+Re;    
            case {'No Conv','Spy'}
                ConvL = ones(1,N_Grid);
        end

    case 2
        center  = ceil(N_Grid/2);
        [p1,p2] = meshgrid(1:N_Grid,1:N_Grid);
        
        switch LineShape
            case 'Lorentzian' 
                lnshpf_R =((-1./(-(p2-center)+1i*LineWidth)).*(1./((p1-center)+1i*LineWidth)));
                lnshpf_N =((-1./( (p2-center)+1i*LineWidth)).*(1./((p1-center)+1i*LineWidth)));

            case 'Gaussian'
                lnshpf_R = ngaussval(sqrt((p1-center).^2+(p2-center).^2),LineWidth);
                lnshpf_N = ngaussval(sqrt((p1-center).^2+(p2-center).^2),LineWidth);
            case 'KK'
                disp('not support KK lineshape in 2D yet...')
            case {'No Conv','Spy'}
                lnshpf_R = 1;
                lnshpf_N = 1;
                disp('Plotting stick spectrum')    
            otherwise
                disp(['Lineshape: ',LineShape,' is not supported...'])
        end
        ConvL.lnshpf_R = lnshpf_R;
        ConvL.lnshpf_N = lnshpf_N;
        
    otherwise
        disp(['Dimension: ',Dimension,' is not supported...'])
end
