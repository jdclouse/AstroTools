function result = rss(vec)
[nr, nc] = size(vec);
result = zeros(1,nc);
for idx = 1:nc
  sum = 0;
  for elem = 1:nr
    sum += vec(elem,idx)*vec(elem,idx);
  end
  result(idx) = sqrt(sum);
end
%result = sum;
end
  