function matrixCompare( test_str, mat, tol )
%matrixCompare compare the matrix, give a pass/fail message
error = sum(sum(mat));
result = 'FAIL';
if error < tol
    result = 'PASS';
end

fprintf('%s: %s\n', test_str, result);

end

