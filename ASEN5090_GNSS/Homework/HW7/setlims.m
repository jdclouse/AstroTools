function setlims(minimum, maximum)

for ii = 1:3
subplot(1,3,ii)
grid on
ylim([minimum, maximum])
end

end