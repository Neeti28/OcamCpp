struct CostFunctor
{
    CostFunctor (std::vector<double> X, std::vector<double> Y, std::vector<double>Z,
                 std::vector<double>detected_pixel_X, std::vector<double>detected_pixel_Y)
        : X_(X)
        ,Y_(Y)
        ,Z_(Z)
        ,detected_pixel_X_(detected_pixel_X)
        ,detected_pixel_Y_(detected_pixel_Y)
    {
        ss_ = std::reverse(ss.begin(), ss.end());
    }

    bool operator()(std::vector<double> extrinsic_params, std::vector<double>affine_parameters, std::vector<double> ss, double* residuals) const
    {

        Eigen::Matrix3d R = rodrigues([extrinsic_params[offset+1+lauf],extrinsic_params[offset+2+lauf],extrinsic_params[offset+3+lauf]]);
        Eigen::Array3d T;
        T<< extrinsic_params[offset+4+lauf],extrinsic_params[offset+5+lauf],extrinsic_params[offset+6+lauf];
        std::vector<double> roots;
        for(unsigned i =0; i<X.size(); ++i )
        {
            Eigen::MatrixXd corner;
            corner << X[i], Y[i], Z[i];
            Eigen::MatrixXd corner_ccd = R*corner + T;
            if(corner_ccd[0][0] == 0 && corner_ccd[1][0] == 0 )
            {
                corner_ccd[0][0] = eps;
                corner_ccd[1][0] = eps;
            }
            double m = corner_ccd[2][0]/sqrt(pow(corner_ccd[0][0], 2) + pow(corner_ccd[1][0], 2));
            std::vector<double> poly_coef_tmp_ = ss;
            poly_coef_tmp_[poly_coef_tmp_.size() -2] = poly_coef_tmp_[poly_coef_tmp_.size() -2] - m;
            if(poly_coef_tmp_.size() == 5)
            {
                roots = getQuatricRoots(poly_coef_tmp[0], poly_coef_tmp[1], poly_coef_tmp[2],poly_coef_tmp[3], poly_coef_tmp[4]);
            }
            else
                std::cout<<"Need to adapt the code for taylor order more than 4."<<std::endl;
            if(res.size()=0)
            {
                rho = nan;
            }
            else if (res.size() >1) {
                rho = min(res);
            }
            else
                rho = res;
            double tmp_x = (corner_ccd[0][0]/sqrt(pow(corner_ccd[0][0], 2) + pow(corner_ccd[1][0], 2)) * rho);
            double tmp_y = (corner_ccd[1][0]/sqrt(pow(corner_ccd[0][0], 2) + pow(corner_ccd[1][0], 2)) * rho);
            double a = affine_parameters[0];
            double b = affine_parameters[1];
            double c = affine_parameters[2];
            double d = affine_parameters[3];
            double e = affine_parameters[4];
            double xc = affine_parameters[5];
            double yc = affine_parameters[6];

            double pixel_x = tmp_x*c + tmp_y*d + xc*a;
            double pixel_y = tmp_x*e + tmp_y + yc*b;
        }



        return true;
    }
    std::vector<double>X_;
    std::vector<double>Y_;
    std::vector<double>Z_;
    std::vector<double>detected_pixel_X_;
    std::vector<double>detected_pixel_Y_;
};

void optimizer(std::vector<Eigen::MatrixXd> board_corners_x, std::vector<Eigen::MatrixXd> board_corners_y, std::vector<Eigen::MatrixXd> board_corners_z
               ,std::vector<Eigen::MatrixXd> detected_corners_x, std::vector<Eigen::MatrixXd> detected_corners_y, ocam_model)
{
    ceres::Problem problem;
    for(unsigned board_num=0; board_num<num_checkerbords; ++board_num)
    {
        std::vector<double> board_corner_x = board_corners_x[board_num];
        std::vector<double> board_corner_y = board_corners_y[board_num];
        std::vector<double> board_corner_z = board_corners_z[board_num];
        std::vector<double> detected_corner_x = detected_corners_x[board_num];
        std::vector<double> detected_corner_y = detected_corners_x[board_num];

        DynamicNumericDiffCostFunction<CostFunctor>* cost_function =
                new DynamicNumericDiffCostFunction<CostFunctor> (new CostFunctor (board_corner_x, board_corner_y, board_corner_z,
                                                                                  detected_corner_x, detected_corner_y));

        problem.AddResidualBlock(cost_function, NULL, anchor, marker);

    }
}
