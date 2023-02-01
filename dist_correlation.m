function dcor = dist_correlation(xx, yy)
%% Program starts from here
% Check if the sizes of the inputs match
if size(xx,1) ~= size(yy,1)
    error('Inputs must have the same number of rows')
end
nn = length(xx);
X = repmat(xx, 1, nn);
Y = repmat(yy, 1, nn);
a_jk = abs(X-X');
b_jk = abs(Y-Y');
a_bar_j = mean(a_jk, 2);
a_bar_k = mean(a_jk, 1);
a_bar_j_mtrx = repmat(a_bar_j, 1, nn);
a_bar_k_mtrx = repmat(a_bar_k, nn, 1);
a_mean = mean(mean(a_jk));
b_bar_j = mean(b_jk, 2);
b_bar_k = mean(b_jk, 1);
b_bar_j_mtrx = repmat(b_bar_j, 1, nn);
b_bar_k_mtrx = repmat(b_bar_k, nn, 1);
b_mean = mean(mean(b_jk));
A_jk = a_jk - a_bar_j_mtrx - a_bar_k_mtrx + a_mean;
B_jk = b_jk - b_bar_j_mtrx - b_bar_k_mtrx + b_mean;
% Calculate squared sample distance covariance and variances
dcov = sum(sum(A_jk.*B_jk));
dvarx = sum(sum(A_jk.*A_jk));
dvary = sum(sum(B_jk.*B_jk));
% Calculate the distance correlation
dcor = sqrt(dcov/sqrt(dvarx*dvary));
end
