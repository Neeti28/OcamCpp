class calib_data
{
public:
    calib_data() {}
    float xc;
    float yc;
    float dx;
    float dy;
    unsigned n_sq_x;
    unsigned n_sq_y;
    unsigned taylor_order;
};


struct ocam_model
{
    std::vector<double>ss;
    float xc;
    float yc;
    unsigned width;
    unsigned height;
    float c;
    float d;
    float e;
};
