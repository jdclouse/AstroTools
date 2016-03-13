stk.v.10.0
WrittenBy    STK_v10.1.3

BEGIN Planet

Name        Jupiter

BEGIN PathDescription

		CentralBody                    Jupiter
		UseCbEphemeris                 No

BEGIN               EphemerisData

	EphemerisSource             JplDE

	JplIndex                    4

	JplSpiceId                  599

	ApplyTDTtoTDBCorrectionForDE     Yes

      OrbitEpoch                  2451545.0000000
      OrbitMeanDist               778412023132.79
      OrbitEcc                    0.048392660000000
      OrbitInc                    23.233623726315
      OrbitRAAN                   3.2543575060730
      OrbitPerLong                15.019945665615
      OrbitMeanLong               34.670475665615
      OrbitMeanDistDot            2487.6456756167
      OrbitEccDot                 -3.5263518138261e-009
      OrbitIncDot                 -2.0487360627962e-007
      OrbitRAANDot                -1.4970999423279e-007
      OrbitPerLongDot             6.3732431943191e-006
      OrbitMeanLongDot            0.083086747568238

END     EphemerisData

END PathDescription

	BEGIN PhysicalData

		GM                   1.266865350000e+017
		Radius               7.149200000000e+007
		Magnitude            0.000000000000e+000
		ReferenceDistance    0.000000000000e+000

	END PhysicalData

	BEGIN AutoRename

		AutoRename           Yes

	END AutoRename

BEGIN Extensions
    
    BEGIN Graphics

			BEGIN Attributes

				MarkerColor             #ffff00
				LabelColor              #ffff00
				LineColor               #ffff00
				LineStyle               0
				LineWidth               1.0
				MarkerStyle             2
				FontStyle               0

			END Attributes

			BEGIN Graphics

				Show                     On
				Inherit                  Off
				ShowLabel                On
				ShowPlanetPoint          On
				ShowSubPlanetPoint       Off
				ShowSubPlanetLabel       Off
				ShowOrbit                On
				NumOrbitPoints           360
				OrbitTime                3.749769292471e+008
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
    END Desc
    
    BEGIN Crdn
    END Crdn
    
    BEGIN VO
    END VO

END Extensions

END Planet

