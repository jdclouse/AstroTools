%==============================================================================================%
%==============================================================================================%
%                                                                                              %
% README file for Sidera GPS Visibility GUI                                                    %
%                                                                                              %
% Author: Ben K. Bradley                                                                       %
% Date: 02/16/2013                                                                             %
%                                                                                              %
% ASTER Labs, Inc.                                                                             %
%                                                                                              %
%==============================================================================================% 
%==============================================================================================%

* The driver program is called Sidera
	- To RUN, type   >> Sidera    into the MATLAB
          command window and press enter/return.	
	- This will run the GUI (graphical user interface) that
          allows you to enter in the scenario settings and run
          the visibility computations.
	


* Most of the files are .p files instead of the usual .m files.
	- This simply indicates that these files are encrypted
          using MATLAB's encryption algorithm. 



* Each time Sidera is run the computed visibility data is is saved
	as GPSvis_Data#.mat where # is incremented as necessary.


* Variables Saved:

  time_wntow    - week number and time of week (seconds) for all computations  [WN TOW] (nx2) (WN wrt 06jan80 (no rollover))
  GPSdata       - structure of GPS results. structure is below

             		.Xecef: ECEF X-coordinates of healthy satellites (m)   (nx32)
             		.Yecef: ECEF Y-coordinates of healthy satellites (m)   (nx32)
             		.Zecef: ECEF Z-coordinates of healthy satellites (m)   (nx32)
         	    .Xecef_vis: ECEF X-coordinates of visible satellites (m)   (nx32)(visible to antenna)
         	    .Yecef_vis: ECEF Y-coordinates of visible satellites (m)   (nx32)(visible to antenna)
         	    .Zecef_vis: ECEF Z-coordinates of visible satellites (m)   (nx32)(visible to antenna)
      		 .topo_numsats: number of satellites above topomask             (nx1)
                      .topo_az: topocentric azimuth of satellites, deg         (nx32)
           	      .topo_el: topocentric elevation of satellites, deg       (nx32)
       		  .ant_numsats: number of satellites above antenna mask         (nx1)
            	       .ant_el: satellite elevation wrt antenna, deg           (nx32)
     		.ant_dopplerL1: doppler shift on L1 frequency, Hz              (nx32)
     		.ant_dopplerL2: doppler shift on L2 frequency, Hz              (nx32)
     		.ant_dopplerL5: doppler shift on L5 frequency, Hz              (nx32)
               		  .DOP: structure containing DOPs:  .GDOP  .HDOP   each (nx1)  (Dilution of Precision)
                                                            .VDOP  .PDOP
                                                            .TDOP


