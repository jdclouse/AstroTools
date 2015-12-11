%% Problem 3b)
Ab = [8     2     6     7     2     3     7
    8     2     3     9     3    10     5
    3     5   8    10     9     4     4
    7    10     3     6     3     2     9
    7     4     6     2     9     3     6];

[Q,R] = qr(Ab)
rank_A = 5;
range_A_basis = Q
[Q_AT, R_AT] = qr(Ab')
null_A_basis = Q_AT(:,rank_A+1:end)