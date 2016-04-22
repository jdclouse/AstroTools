stk.v.10.0
WrittenBy    STK_v10.1.3

BEGIN Planet

Name        Callisto

BEGIN PathDescription

		CentralBody                    Callisto
		UseCbEphemeris                 Yes

BEGIN               EphemerisData

	EphemerisSource             None

	JplIndex                    -1

	JplSpiceId                  504

	ApplyTDTtoTDBCorrectionForDE     Yes

END     EphemerisData

END PathDescription

	BEGIN PhysicalData

		GM                   7.179292000000e+012
		Radius               2.410300000000e+006
		Magnitude            0.000000000000e+000
		ReferenceDistance    0.000000000000e+000

	END PhysicalData

	BEGIN AutoRename

		AutoRename           Yes

	END AutoRename

BEGIN Extensions
    
    BEGIN Graphics

			BEGIN Attributes

				MarkerColor             #0000ff
				LabelColor              #0000ff
				LineColor               #0000ff
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

