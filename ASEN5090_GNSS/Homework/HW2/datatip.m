function datatip( hl, loc )
%datatip make a datatip

% create a line and get its handle
% hl = plot(1:10);
% put the data cursor mode on so
% that a data cursor manager is created
datacursormode on
% get the handle to the data cursor manager
wbdf = get(gcf,'WindowButtonDownFcn');
dcm = wbdf{3}{2};
% Specify the required position as a structure
props.Position = [loc 0]
% create the data tip
dcm.createDatatip(hl,props);
end

