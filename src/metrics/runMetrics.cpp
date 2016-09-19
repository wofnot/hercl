#include <stdio.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <assert.h>

#include "metrics.h"

using namespace std;

int main(int argc, char** argv )
{
  assert(argc>=2);
  
  int i,j;
  
  bool colour = false;
  
  if (argv[1][0] == '-' && argv[1][1] == 'c') {
    colour = true;
  }
  
  for (j = 1; j < argc; j++) {
    char *fn = argv[j];
    
    if (fn[0] == '-') continue;
    
    cerr << "reading " << fn << "\n";
    
    Mat image;
    image = imread(fn, colour ? 1 : 0);

    if ( !image.data )
    {
        cerr << "No image data\n";
        return -1;
    }
    
    Mat image2;
    if(colour) {
   //   cout << "loaded bgr\n";
      cvtColor(image, image2, CV_BGR2HLS );
    } else {
      image2=image;
    } 
    
    //namedWindow("Display Image", WINDOW_AUTOSIZE );
    //imshow("Display Image", image);
    
    metrics_t metrics = getMetrics(image2, colour);
    
    for (i = 0; i < (int) metrics.size(); i++) {
      cout << "," << metrics[i];
    }
    cout << "\n";
  }
  

  return 0;
}
