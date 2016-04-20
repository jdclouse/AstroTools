stk.v.10.0
WrittenBy    STK_v10.1.3

BEGIN Planet

Name        Europa

BEGIN PathDescription

		CentralBody                    Europa
		UseCbEphemeris                 No

BEGIN               EphemerisData

	EphemerisSource             JplSpice

	JplIndex                    -1

	JplSpiceId                  502

	ApplyTDTtoTDBCorrectionForDE     Yes

END     EphemerisData

END PathDescription

	BEGIN PhysicalData

		GM                   3.202720000000e+012
		Radius               1.562600000000e+006
		Magnitude            0.000000000000e+000
		ReferenceDistance    0.000000000000e+000

	END PhysicalData

	BEGIN AutoRename

		AutoRename           Yes

	END AutoRename

BEGIN Extensions
    
    BEGIN Graphics

			BEGIN Attributes

				MarkerColor             #00ffff
				LabelColor              #00ffff
				LineColor               #00ffff
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

