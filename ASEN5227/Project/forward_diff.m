%% forward differentiation of f(t).
% Remember that there will be one fewer points in the deriviative
function deriv = forward_diff(f,t)
% All the delta t's
t_diffs = (t(2:end) - t(1:end-1));

% Take the delta of the function, divide by the appropriate dt.
deriv = f(:,2:end) - f(:,1:end-1);
for ii = 1:length(deriv)
    deriv(:,ii) = deriv(:,ii)/t_diffs(ii);
end