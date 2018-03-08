#include <vector>
#include <math.h>
#include<eigen2/Eigen/Array>
#include<eigen2/Eigen/Dense>
#include<eigen2/Eigen/SVD>

struct complex { double r,i; };
struct pair<T> { T p1, p2; };

pair<complex> getRoots(double a, double b, double c)
{
    pair<complex> result={0};

    if(a<0.000001)    // ==0
    {
        if(b>0.000001)  // !=0
            result.p1.r=result.p2.r=-c/b;
        else
            if(c>0.00001) throw exception("no solutions");
        return result;
    }

    double delta=b*b-4*a*c;
    if(delta>=0)
    {
        result.p1.r=(-b-sqrt(delta))/2/a;
        result.p2.r=(-b+sqrt(delta))/2/a;
    }
    else
    {
        result.p1.r=result.p2.r=-b/2/a;
        result.p1.i=sqrt(-delta)/2/a;
        result.p2.i=-sqrt(-delta)/2/a;
    }

    return result;
}

template <class MatT>
Eigen::Matrix<typename MatT::Scalar, MatT::ColsAtCompileTime, MatT::RowsAtCompileTime>
pseudoinverse(const MatT &mat, typename MatT::Scalar tolerance = typename MatT::Scalar{1e-4}) // choose appropriately
{
    typedef typename MatT::Scalar Scalar;
    auto svd = mat.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
    const auto &singularValues = svd.singularValues();
    Eigen::Matrix<Scalar, MatT::ColsAtCompileTime, MatT::RowsAtCompileTime> singularValuesInv(mat.cols(), mat.rows());
    singularValuesInv.setZero();
    for (unsigned int i = 0; i < singularValues.size(); ++i) {
        if (singularValues(i) > tolerance)
        {
            singularValuesInv(i, i) = Scalar{1} / singularValues(i);
        }
        else
        {
            singularValuesInv(i, i) = Scalar{0};
        }
    }
    return svd.matrixV() * singularValuesInv * svd.matrixU().adjoint();
}

void findExtrinsicParameters(std::vector<Eigen::ArrayXd >board_corners_x, std::vector<std::ArrayXd >board_corners_y,
                             std::vector<Eigen::ArrayXd >detected_corners_x, std::vector<Eigen::ArrayXd >detected_corners_y,
                             double xc, double yc, unsigned taylor_order)
{
    std::vector<Eigen::Matrix3d> RRfin;
    unsigned num_checkerboard = detected_corners_x.size();
    for(unsigned im_num=0; im_num<num_checkerboard; ++im_num)
    {
        Eigen::ArrayXd detected_corner_x = detected_corners_x[im_num] - xc;
        Eigen::ArrayXd detected_corner_y = detected_corners_y[im_num] - yc;

        Eigen::ArrayXd board_corner_x = board_corners_x[im_num];
        Eigen::ArrayYd board_corner_y = board_corners_y[im_num];

        //building matrix
        Eigen::MatrixXd A(detected_corner_x.rows(), 6);
        Eigen::ArrayXd xy = board_corner_x * detected_corner_y;
        Eigen::ArrayXd y = board_corner_y * detected_corner_y;
        Eigen::ArrayXd x = -1 * board_corner_x * detected_corner_x;
        Eigen::ArrayXd yx = -1 * board_corner_y * detected_corner_x;
        Eigen::ArrayXd Y = detected_corner_y;
        Eigen::ArrayXd X = -1 * detected_corner_x;

        for(unsigned i=0; i< xy.size(); ++i)
        {
            A[i][0] = xy[i];
            A[i][1] = y[i];
            A[i][2] = x[i];
            A[i][3] = yx[i];
            A[i][4] = Y[i];
            A[i][5] = X[i];
        }

        //SVD
        Eigen::JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);
        Eigen::MatrixXd V = svd.matrixV();

        double R11 = V[0][V.col() -1 ];
        double R12 = V[1][V.col() -1 ];
        double R21 = V[2][V.col() -1 ];
        double R22 = V[3][V.col() -1 ];
        double T1 = V[4][V.col() -1 ];
        double T2 = V[5][V.col() -1 ];

        double AA = pow(R11*R12 + R21*R22, 2);
        double BB = pow(R11, 2) + pow(R21, 2);
        double CC = pow(R12, 2) + pow(R22, 2);

        //Finding R31 and R32
        pair<complex> R32_2_temp = getRoots(1, CC-BB, -AA);

        std::vector<double>R32_2;
        if (R32_2_temp.p1.i == 0 & R32_2_temp.p1.r >=0)
            R32_2.push_back(R32_2_temp.p1.r);
        if(R32_2_temp.p2.i == 0 & R32_2_temp.p2.r >=0)
            R32_2.push_back(R32_2_temp.p2.r);

        std::vector<double>R31;
        std::vector<double>R32;
        std::vector<int>sg(1, -1);
        for(i =0;i<R32_2.size();++i)
        {
            for(j=0;j<sg.size();++j)
            {
                double sqrtR32_2 = sg[j] * sqrt(R32_2[i]);
                R32.push_back(sqrtR32_2);
                if (R32_2 == 0)
                {
                    R31.push_back(sqrt(CC-BB));
                    R31.push_back(-sqrt(CC-BB));
                    R32.push_back(sqrtR32_2);
                }
                else
                {
                    R31.push_back(-sqrtR32_2/)(R11*R12+R21*R22);
                }
            }
        }

        std::vector<Eigen::Matrix3d> RR;
        for(unsigned i1=0;i1<R32.size();++i1)
        {
            for(unsigned i2=0;i2<sg.size(); ++i2)
            {
                double Lb=sqrt(R11^2+R21^2+R31(i1)^2);
                Eigen::Matrix3d RR_tmp;
                double factor = sg(i2) * Lb ;
                RR_tmp<<R11, R12, T1,
                        R21, R22, T2,
                        R31[i1], R32[i1], 0;
                RR_tmp = factor * RR_tmp;
                RR.push_back(RR_tmp);
            }
        }

        std::vector<Eigen::Matrix3d>RR1;
        double minRR = std::numeric_limits::infinity();
        int minRR_ind = -1;
        for(unsigned min_count=0; min_count<RR.size();++min_count)
        {
            Eigen::Matrix3d RR_matrix = RR[min_count];
            double norm = sqrt(pow(RR_matrix[1][3] - detected_corner_x[1], 2) + pow(RR_matrix[2][3] - detected_corner_y[1], 2));
            if(norm < minRR)
            {
                minRR = norm;
                minRR_ind = min_count;
            }

        }

        if(minRR_ind != -1)
        {
            unsigned count2 = 0;
            for(unsigned count = 1; count <RR.size(); ++count)
            {
                Eigen::Matrix3d RR_matrix = RR[count];
                Eigen::Matrix3d RR_matrix_min = RR[minRR_ind];
                if(signbit(RR_matrix[1][3]) == signbit(RR_matrix_min[1][3]) & signbit(RR_matrix[2][3])==signbit(RR_matrix_min[2][3]))
                {
                    count2 = count2 + 1;
                    RR1.push_back(RR_matrix);
                }
            }
        }

        if(RR1.size()==0)
        {
            double ss = 0;

        }

        RRfin.push_back(RR1);

    }
}

void omni_find_parameters_fun(Xt, Yt, detected_corner_x, detected_corner_y, xc, yc, RRfin, taylor_order, num_checkerboard)
{
    unsigned count = 0;
    for(unsigned im_num=0; im_num<num_checkerboard; ++im_num)
    {
        Eigen::ArrayXd detected_corner_x = detected_corners_x[im_num] - xc;
        Eigen::ArrayXd detected_corner_y = detected_corners_y[im_num] - yc;

        count=count+1;
        Eigen::Matrix3d RRdef=RRfin[im_num];

        double R11 = RRdef[1][1];
        double R21 = RRdef[2][1];
        double R31 = RRdef[3][1];
        double R12 = RRdef[1][2];
        double R22 = RRdef[2][2];
        double R32 = RRdef[3][2];
        double T1 = RRdef[1][3];
        double T2 = RRdef[2][3];

        Eigen::ArrayXd MA= R21*board_corner_x + R22*board_corner_y + T2;
        Eigen::ArrayXd MB= Ypt*( R31*board_corner_x + R32*board_corner_y );
        Eigen::ArrayXd MC= R11*board_corner_x + R12*board_corner_y + T1;
        Eigen::ArrayXd MD= Xpt*( R31*board_corner_x + R32*board_corner_y );

        std::vector<Eigen::ArrayXd>rho;
        for(unsigned j=2; j<taylor_order; ++j)
        {
            Eigen::ArrayXd tmp = detected_corner_x * detected_corner_x + detected_corner_y * detected_corner_y;
            tmp = tmp.sqrt();
            tmp = tmp.pow(j);
            rho.push_back(tmp);
        }

        std::vector<Eigen::ArrayXd> PP1;
        PP1.push_back(MA);
        PP1.push_back(MC);

        for(unsigned j=2; j<taylor_order; ++j)
        {
            PP1.push_back(MA * rho[j]);
            PP1.push_back(MC * rho[j]);
        }


        PP=[PP   zeros(size(PP,1),1);
            PP1, zeros(size(PP1,1),count-1) [-Ypt;-Xpt]];
        QQ=[QQ;
            MB; MD];

    }

}
unsigned plot_RR(RR, board_corner_x, board_corner_y, detected_corner_x, detected_corner_y, figure_number)
{

    for(unsigned i=1; i<RR.size(); ++i)
    {
        Eigen::Matrix3d RRdef = RR[i];
        double R11 = RRdef[1][1];
        double R21 = RRdef[2][1];
        double R31 = RRdef[3][1];
        double R12 = RRdef[1][2];
        double R22 = RRdef[2][2];
        double R32 = RRdef[3][2];
        double T1 = RRdef[1][3];
        double T2 = RRdef[2][3];

        Eigen::ArrayXd MA= R21*board_corner_x + R22*board_corner_y + T2;
        Eigen::ArrayXd MB= Ypt*( R31*board_corner_x + R32*board_corner_y );
        Eigen::ArrayXd MC= R11*board_corner_x + R12*board_corner_y + T1;
        Eigen::ArrayXd MD= Xpt*( R31*board_corner_x + R32*board_corner_y );
        double rho= sqrt(pow(detected_corner_x,2) + pow(detected_corner_y,2));
        double rho2= pow(detected_corner_x,2) + pow(detected_corner_y,2);


    }
}
