void omni3d2pixel(ss, xx_vec, width, height)
{
    //convert 3D coordinates vector into 2D pixel coordinates
    //These three lines overcome problem when xx = [0,0,+-1]

    std::vector<double> poly_coef = std::reverse(ss.begin(),ss.end());
    std::vector<double> poly_coef_tmp = poly_coef;
    std::vector<double> x;
    std::vector<double> y;
    for(unsigned j = 0;j<xx.size();++j)
    {
        std::vector<double> res;
        Eigen::Array3d xx = xx_vec[j];
        if (xx[0] == 0 && x[1] == 0)
        {
            xx[0] = eps;
            xx[1] = eps;
        }
        m = xx[3]/sqrt(pow(xx[0], 2) + pow(xx[1], 2));
        poly_coef_tmp[poly_coef_tmp.size() - 2] = poly_coef[poly_coef_tmp.size() - 2] - m;
        if(poly_coef_tmp.size() == 5)
        {
            res = root_fast(poly_coef_tmp[0], poly_coef_tmp[1], poly_coef_tmp[2],poly_coef_tmp[3], poly_coef_tmp[4]);
        }
        else
            std::cout<<"Need to adapt the code for taylor order more than 4."<<std::endl;
    }
    if(res.size()=0)
    {
        rho = nan;
    }
    else if (res.size() >1) {
       rho = min(res);
    }
    else
        rho = res;
    x.push_back(xx[0]/sqrt(pow(xx[0],2)+ pow(xx[1], 2)) * rho);
    y.push_back(xx[1]/sqrt(pow(xx[0],2)+ pow(xx[1], 2)) * rho);
}
