function [rootMeanSquaredError] = rmse(truth_angles,estimated_angles)

n = length(estimated_angles);

if numel(truth_angles) ~= numel(estimated_angles)
    disp('User, number of truth_angles ~= estimated_angles');
    return;
end

rootMeanSquaredError = sqrt( (1./n) .* sum( (truth_angles - estimated_angles).^2 ) );

end
