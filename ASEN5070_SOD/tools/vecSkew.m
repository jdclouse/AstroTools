function mat = vecSkew( vec )
fcnPrintQueue(mfilename('fullpath'))
mat = [0, -vec(3), vec(2); vec(3), 0, -vec(1); -vec(2), vec(1), 0];

end

