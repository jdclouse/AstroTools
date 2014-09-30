classdef PRNCode < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        G1;
        G2;
        G1_out = [];
        G2_out = [];
        S1;
        S2;
        CA_code = [];
        PRN_Number = -1;
    end
    
    methods
        function PRN = PRNCode(PRN_number)
            fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 
            PRN.G1 = BitShiftRegister(10, [3,10]);
            PRN.G2 = BitShiftRegister(10, [2, 3, 6, 8, 9, 10]);
            
            %phase selection data obtained from
            %IS-GPS-200D
            phase_selections = [2 6;3 7;4 8;5 9;1 9;2 10;1 8;2 9;3 10;...
                2 3;3 4;5 6;6 7;7 8;8 9;9 10;1 4;2 5;3 6;4 7;5 8;6 9;...
                1 3;4 6;5 7;6 8;7 9;8 10;1 6;2 7;3 8;4 9;5 10;4 10;...
                1 7;2 8;4 10];
            PRN.S1 = phase_selections(PRN_number, 1);
            PRN.S2 = phase_selections(PRN_number, 2);
            PRN.PRN_Number = PRN_number;
        end
        function update(PRN)
            G1_out_scalar = PRN.G1.update();
            
            %TODO: consider reversing this?
            PRN.CA_code = [ PRN.CA_code recurse_xor([PRN.G2.bits(PRN.S1), ...
                PRN.G2.bits(PRN.S2), G1_out_scalar])];
            
            G2_out_scalar = PRN.G2.update();
            PRN.G1_out = [G1_out_scalar PRN.G1_out];
            PRN.G2_out = [G2_out_scalar PRN.G2_out];
        end
    end
    
end

