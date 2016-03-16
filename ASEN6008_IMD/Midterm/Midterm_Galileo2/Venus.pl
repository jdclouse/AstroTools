stk.v.10.0
WrittenBy    STK_v10.1.3

BEGIN Planet

Name        Venus

BEGIN PathDescription

		CentralBody                    Venus
		UseCbEphemeris                 No

BEGIN               EphemerisData

	EphemerisSource             JplDE

	JplIndex                    1

	JplSpiceId                  299

	ApplyTDTtoTDBCorrectionForDE     Yes

      OrbitEpoch                  2451545.0000000
      OrbitMeanDist               108208925006.86
      OrbitEcc                    0.0067732300000002
      OrbitInc                    24.432965702573
      OrbitRAAN                   8.0077488974455
      OrbitPerLong                132.21942759358
      OrbitMeanLong               182.66617759357
      OrbitMeanDistDot            3.7681051444216
      OrbitEccDot                 -1.3519507186858e-009
      OrbitIncDot                 4.1244766351354e-007
      OrbitRAANDot                -4.3180134699804e-007
      OrbitPerLongDot             -8.5280571764923e-007
      OrbitMeanLongDot            1.6021304488902

END     EphemerisData

END PathDescription

	BEGIN PhysicalData

		GM                   3.248585920790e+014
		Radius               6.051800000000e+006
		Magnitude            0.000000000000e+000
		ReferenceDistance    0.000000000000e+000

	END PhysicalData

	BEGIN AutoRename

		AutoRename           Yes

	END AutoRename

BEGIN Extensions
    
    BEGIN Graphics

			BEGIN Attributes

				MarkerColor             #ff00ff
				LabelColor              #ff00ff
				LineColor               #ff00ff
				LineStyle               0
				LineWidth               1.0
				MarkerStyle             2
				FontStyle               0

			END Attributes

			BEGIN Graphics

				Show                     On
				Inherit                  On
				ShowLabel                On
				ShowPlanetPoint          On
				ShowSubPlanetPoint       Off
				ShowSubPlanetLabel       Off
				ShowOrbit                On
				NumOrbitPoints           360
				OrbitTime                0.000000000000e+000
				OrbitDisplay             OneOrbit
				TransformTrajectory      On

			END Graphics
    END Graphics
    
    BEGIN ExternData
    END ExternData
    
    BEGIN ADFFileData
    END ADFFileData
    
    BEGIN AccessConstraints
		LineOfSight   IncludeIntervals 
    END AccessConstraints
    
    BEGIN Desc
        ShortText    0

        LongText    0

    END Desc
    
    BEGIN Crdn
    END Crdn
    
    BEGIN VO
    END VO

END Extensions

END Planet

