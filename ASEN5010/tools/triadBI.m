function DCM_BI = triadBI( body_vec_primary, body_vec_secondary,...
    inrtl_vec_primary, inrtl_vec_secondary )
fcnPrintQueue(mfilename('fullpath'))
%Assumes input vectors are unitized!

DCM_BI = triadDCM(body_vec_primary, body_vec_secondary)...
    *(triadDCM(inrtl_vec_primary, inrtl_vec_secondary))';
end

