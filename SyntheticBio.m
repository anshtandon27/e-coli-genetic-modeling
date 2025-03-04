
%Identify steady state GFP expression
function normGFP = SyntheticBio(AHL, params)
    %Input parameterse and parameters to modify
    LuxR = params(1);
    rho_R = 0.5;
    delta_R = 2.31e-1;
    KR = params(4);
    alpha_TXGFP = 0.05;
    delta_TXGFP = 0.2;
    alpha_GFP = 2;
    delta_GFP = 4e-4;
    n1 = 1;
    
    % Solve for steady-state dimeric R concentration
    R = rho_R * LuxR^2 .* AHL.^2 ./ delta_R;
    % Solve for steady-state TXGFP mRNA concentration
    TXGFP = alpha_TXGFP .* (R.^n1 ./ (KR^n1 + R.^n1)) ./ delta_TXGFP;
    % Solve for steady-state GFP protein concentration
    GFP = alpha_GFP .* TXGFP ./ delta_GFP;
    normGFP = GFP/max(GFP);
end