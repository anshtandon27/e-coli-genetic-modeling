function rmse = compute_rmse(params, conc, expData)
    %Compute RMSE given experimental and computational data
    model_output = SyntheticBio(conc, params);
    rmse = sqrt(mean((model_output-expData).^2));
end