void bundleErrUrban(x, calib_data, M, robust)
{
    //global weights     % only used if robust is enabled

    double a=x[0];
    double b=x[1];
    double c=x[2];
    double d=x[3];
    double e=x[4];
    unsigned offset = calib_data.taylor_order+5;

    ssc=x(6:offset);

    //Each matrix is one checkerboard
    Xpp=[];
    Ypp=[];
    unsigned lauf = 0;
    for(unsigned board_num=0; board_num<M.size(); ++board_num)
    {
        std::vector<Eigen::Array3d> Mc;
        Eigen::Matrix3d R = rodrigues([x[offset+1+lauf],x[offset+2+lauf],x[offset+3+lauf]]);
        Eigen::Array3d T;
        T<< x[offset+4+lauf],x[offset+5+lauf],x[offset+6+lauf];
        num_points = M[board_num].size();

        Mc.push_back(R*M[im_num] + T*ones(1,num_points));
        Xpp=[Xpp; calib_data.Xp_abs{i}];
        Ypp=[Ypp; calib_data.Yp_abs{i}];
        lauf = lauf+6;
    }

    [xp1,yp1] = omni3d2pixel(calib_data.ocam_model.ss.*ssc, Mc, calib_data.ocam_model.width, calib_data.ocam_model.height);
    xp = xp1*c + yp1*d + calib_data.ocam_model.xc*a;
    yp = xp1*e + yp1   + calib_data.ocam_model.yc*b;

    errx = Xpp-xp;
    erry = Ypp-yp;

    if (!robust)
    {
        errW = [errx
                erry];
    }
    else
    {
        errn = [errx
                erry];
        w = huberWeight(errn);
        weights = w;
        errW = sqrt(w).*errn;
    }

}
double huberWeight(v)
{
    k = 1;
    a = abs(v) <= k;
    b = abs(v) > k;
    weight(find(a)) = v(a)./v(a);
    weight(find(b)) = k./abs(v(b));
    weight = weight';
}
