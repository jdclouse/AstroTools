function DCM_frame_to_T = triadDCM( vec_primary, vec_secondary )
fcnPrintQueue(mfilename('fullpath'))
%Assumes input vectors are unitized!
%input vectors in the same frame
t1 = vec_primary;
t2 = cross(t1, vec_secondary)/norm(cross(t1, vec_secondary));
t3 = cross(t1, t2);

% Put the T-frame vectors, expressed in the given frame, into the DCM:
DCM_frame_to_T = [t1 t2 t3];

end

