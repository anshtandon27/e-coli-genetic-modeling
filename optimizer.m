function paramOptim = optimizer(conc, expData, params)
    param0 = params;

    %Set boundaries for biologically feasible values
    lb = 0.1*param0;
    ub = 10*param0;
    objFun = @(x) compute_rmse(x, conc, expData);
    options = optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'sqp');

    %call optimizer function, fmincon
    param_opt = fmincon(objFun, param0, [], [], [], [], lb, ub, [], options);
    paramOptim = param_opt;
end
