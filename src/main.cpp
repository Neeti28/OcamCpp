bool readCorners(const string& filename, std::vector<cv::Point2f>& image_points)
{
    image_points.clear ();
    try
    {
        ifstream inputFile;

        //inputFile.exceptions (!ios_base::goodbit);
        inputFile.open(filename.c_str ());
        unsigned point_count;
        inputFile >> point_count;
        image_points.resize(point_count);
        for(unsigned pointIdx = 0; pointIdx < point_count; ++pointIdx)
        {
            inputFile >> image_points [pointIdx].x;
            inputFile >> image_points [pointIdx].y;
            //cout << pointIdx << ". " << image_points [pointIdx] << endl;
        }
        inputFile.close();
    }
    catch (std::ifstream::failure e)
    {
        return false;
    }
    return true;
}

int main (int argc, char** argv)
{
    po::options_description desc("Options");
    desc.add_options()
            ("help,h", "Print help messages")
            //why do i need this
            ("imfolder, i", po::value<string>(), "folder containing images" )
            ("corners, c", po::value<string>(), "folder containing corners")
            ("board, b", po::value<string>(), "checkerboard xm file")
            ("taylor, t", po::value<unsigned>(), "taylor order")
            ("verbodśe,v", po::value<unsigned>(),  "verbose")
            ("frames", po::value<vector<string> > (), "frames to be compared");

    po::variables_map var_map;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), var_map);

    if (var_map.count("help"))
    {
        cout << desc << "\n";
        return 1;
    }
    try
    {
        po::notify(var_map);
    } catch (const std::exception& ex)
    {
        cerr << ex.what ()<< endl;
        cerr << desc << endl;
        exit (1);
    }

    string corner_dir = var_map["corners"].as<string> ();
    string image_dir = var_map["imfolder"].as<string> ();
    string board_filename = var_map["board"].as<string> ();
    unsigned verbose_level = var_map["verbose"].as<unsigned> ();
    unsigned taylor_order = var_map["taylor"].as<unsigned> ();

    CheckerBoard board (Size(1, 1), Size2f (0, 0));
    if (!CheckerBoardXMLReader::readBoard (board_filename, board))
    {
        cerr << "Could not read calibration pattern description file \"" << board_filename << "\"\n";
        exit (-1);
    }

    if (verbose_level)
    {
        cout << "Calibration Camera with the following parameters:\n";
        cout << "board size      : " << board.getPatternSize ().width << " x " << board.getPatternSize ().height << endl;
        cout << "pattern size    : " << board.getRectSize ().width << " x " << board.getRectSize ().height << endl;
        cout << "output file name: " << output_filename << endl;
        cout << "verbose level   : " << verbose_level << endl;
    }

    //read corners
    std::vector<Eigen::MatrixXd> detected_corners_x;
    std::vector<Eigen::MatrixXd> detected_corners_y;
    std::vector<Eigen::MatrixXd> board_corners_x;
    std::vector<Eigen::MatrixXd> board_corners_y;
    std::vector<Eigen::MatrixXd> board_corners_z;
    std::vector<std::string> corner_names;
    unsigned num_corners;
    getFiles (corner_dir, "*.corners", corner_names);
    for(unsigned i=0; i<corner_names.size();++i)
    {
        std::vector<cv::Point2f> images_points;
        readCorners(corner_names[i], image_points);
        for(unsigned col=0;col<board.getPatternSize ().width; ++col)
        {
            for(unsigned row=0;row<board.getPatternSize().height; ++row)
            {
                detected_corners_x[row][col] = images_points[col*row + board.getPatternSize().height].x;
                detected_corners_x[row][col] = images_points[col*row + board.getPatternSize().height].y;
                board_corners_x[row][col] = board.getRectSize().height * row;
                board_corners_y[row][col] = board.getRectSize().width * col;
                board_corners_z[row][col] = 0;
            }
        }

    }
    findExtrinsicParameters(board_corners_x, board_corners_y, detected_corners_x, detected_corners_y, ocam_model.xc, ocam_model.yc, calib_data.taylor_order);
    optimizer(board_corners_x, board_corners_y, board_corners_z, detected_corners_x, detected_corners_x, detected_corners_y, ocam_model)
    //get all images
    return 0;
}
