stk.v.10.0
WrittenBy    STK_v10.1.3

BEGIN Planet

Name        Mars

BEGIN PathDescription

		CentralBody                    Mars
		UseCbEphemeris                 Yes

BEGIN               EphemerisData

	EphemerisSource             None

	JplIndex                    3

	JplSpiceId                  499

	ApplyTDTtoTDBCorrectionForDE     Yes

      OrbitEpoch                  2451545.0000000
      OrbitMeanDist               227936636175.28
      OrbitEcc                    0.093412330000000
      OrbitInc                    24.677208276709
      OrbitRAAN                   3.3758259392290
      OrbitPerLong                336.33376729917
      OrbitMeanLong               355.74624729917
      OrbitMeanDistDot            -295.75529617248
      OrbitEccDot                 3.2585900068446e-009
      OrbitIncDot                 4.8386904489709e-008
      OrbitRAANDot                -7.4964810511602e-007
      OrbitPerLongDot             1.1805536709966e-005
      OrbitMeanLongDot            0.52403297064432

END     EphemerisData

END PathDescription

	BEGIN PhysicalData

		GM                   4.282837190128e+013
		Radius               3.396190000000e+006
		Magnitude            0.000000000000e+000
		ReferenceDistance    0.000000000000e+000

	END PhysicalData

	BEGIN AutoRename

		AutoRename           Yes

	END AutoRename

BEGIN Extensions
    
    BEGIN Graphics

			BEGIN Attributes

				MarkerColor             #ff0000
				LabelColor              #ff0000
				LineColor               #ff0000
				LineStyle               0
				LineWidth               1.0
				MarkerStyle             2
				FontStyle               0

			END Attributes

			BEGIN Graphics

				Show                     On
				Inherit                  On
				ShowLabel                Off
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

