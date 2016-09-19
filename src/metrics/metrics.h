/**
Metrics (scaled to [0,1] unless otherwise noted)

"mean_l"        // l for lightness
"stddev_l"
"mean_h"      
"stddev_h"
"mean_s"
"stddev_s"
"entropy_grey"
"mean_edge"
"stddev_edge"
"mean_corner"
"n_corners"     // not scaled
"n_patches"     // up to 5, not scaled
"patch_1_mean_l"
"patch_2_mean_l"
"patch_3_mean_l"
"patch_4_mean_l"
"patch_5_mean_l"
...             // as above for h and s
"patch_1_size"
"patch_2_size"
"patch_3_size"
"patch_4_size"
"patch_5_size"
"convexity"
**/

//const int num_metrics = 19;


#include <vector>

// this is a mess, clean it up


using namespace cv;
using namespace std;

typedef vector<double> metrics_t;

metrics_t getMetrics(Mat image, bool colour = false);
