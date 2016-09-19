#include <stdio.h>
#include <iostream>
#include <opencv2/opencv.hpp>

#include <assert.h>
#include <unordered_map>
#include <utility>

#include "metrics.h"

using namespace std;
using namespace cv;

typedef pair<double, Scalar> patch;

bool compareContours(const patch & a, const patch & b) { 
  return a.first > b.first; // descending sort 
}

const int num_metrics_colour = 33;
string metric_names_colour[num_metrics_colour] = {
  "mean_l",
  "stddev_l",
  "mean_s",
  "stddev_s",
  "mean_h",
  "stddev_h",
  "entropy_grey",
  "mean_edge",
  "stddev_edge",
  "n_patches",
  "patch_1_mean_l",
  "patch_2_mean_l",
  "patch_3_mean_l",
  "patch_4_mean_l",
  "patch_5_mean_l",
  "patch_1_mean_s",
  "patch_2_mean_s",
  "patch_3_mean_s",
  "patch_4_mean_s",
  "patch_5_mean_s",
  "patch_1_mean_h",
  "patch_2_mean_h",
  "patch_3_mean_h",
  "patch_4_mean_h",
  "patch_5_mean_h",
  "patch_1_size",
  "patch_2_size",
  "patch_3_size",
  "patch_4_size",
  "patch_5_size",
  "convexity",
  "mean_corner",
  "n_corners"
};

const int num_metrics_grey = 19;
string metric_names_grey[num_metrics_grey] = {
  "mean_l",
  "stddev_l",
  "entropy_grey",
  "mean_edge",
  "stddev_edge",
  "n_patches",
  "patch_1_mean_l",
  "patch_2_mean_l",
  "patch_3_mean_l",
  "patch_4_mean_l",
  "patch_5_mean_l",
  "patch_1_size",
  "patch_2_size",
  "patch_3_size",
  "patch_4_size",
  "patch_5_size",
  "convexity",
  "mean_corner",
  "n_corners"
};

metrics_t getMetrics(Mat img_in, bool colour) {

  unordered_map <string, double > metrics;

  Mat image;
  
  Mat grey;
  
  if (img_in.type() != CV_8UC1) {
    // colour input
    // assume hls
    assert(img_in.type() == CV_8UC3);
    
    // extract lightness channel
    Mat hls_channels[3];
    split( img_in, hls_channels );
    grey = hls_channels[1];
    
    if (colour) {
      image = img_in;
    } else {
      image = grey;
    }
  } else {
    // greyscale input
    grey = img_in;
    
    if (colour) {
      // create hls image with hue and saturation set to 0
      Mat hls_channels[3];
      hls_channels[1]=img_in;
      hls_channels[0]=Mat(img_in.rows,img_in.cols, CV_8UC1, Scalar(0));
      hls_channels[2]=Mat(img_in.rows,img_in.cols, CV_8UC1, Scalar(0));
      
      merge(hls_channels, 3, image);
    } else {
      // use image as is
      image = img_in;
    }
  }
  
    
  //imshow("grey",grey);
  //waitKey(0);

  int i,j;
  // cout <<"0\n";
  Mat out;
  // ** MEAN, STDDEV
  Scalar mean;
  Scalar stddev;
  
  meanStdDev(image,mean,stddev);

  if (colour) {  
    metrics[string("mean_h")] = mean[0]/255;
    metrics[string("mean_l")] = mean[1]/255;
    metrics[string("mean_s")] = mean[2]/255;
    
    metrics[string("stddev_h")] = stddev[0]/(255*0.5); // (scale according to popoviciu inequality)
    metrics[string("stddev_l")] = stddev[1]/(255*0.5);
    metrics[string("stddev_s")] = stddev[2]/(255*0.5);
  } else {
    metrics[string("mean_l")] = mean[0]/255;
    metrics[string("stddev_l")] = stddev[0]/255;
  }
  
  // ** ENTROPY
  // cout <<"1\n";
  
  Mat hist;
  /// Establish the number of bins
  int histSize = 256;
  
  int n_pix = image.total();
  
  /// Set the ranges ( for B,G,R) )
  float range[] = { 0, 256 } ;
  const float* histRange = { range };
  int channels[] = {0};
  /// Compute the histograms:
  calcHist( &grey, 1, channels, Mat(), hist, 1, &histSize, &histRange);
  hist /= n_pix;

  Mat logP;
  cv::log(hist,logP);
  
  metrics["entropy_grey"] = (-1*sum(hist.mul(logP)).val[0])/log(histSize);
  
  // ** EDGE
  // cout <<"2\n";
  int scale = 1;
  int delta = 0;
  int ddepth = CV_32F;
  /// Generate grad_x and grad_y
  Mat grad_x, grad_y;


  /// Gradient X
  Sobel( image, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT );
  /// Gradient Y
  Sobel( image, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT );

  cv::Mat grad_mag(grad_x.size(), grad_x.type());

  magnitude(grad_x, grad_y, grad_mag);

  // scale matrix elems to [0,1]:
  grad_mag *= 1.0/(sqrt(2)*1020);

  meanStdDev(grad_mag,mean,stddev);

  metrics["mean_edge"] = mean[0];
  metrics["stddev_edge"] = stddev[0]/(0.5); // (scale according to popoviciu inequality)

  // ** CORNERS

  Mat dst, dst_norm, dst_norm_scaled;
  dst = Mat::zeros( image.size(), CV_32FC1 );

  /// Detector parameters
  int blockSize = 2;
  int apertureSize = 3;
  double k = 0.04;
  

  /// Detecting corners
  cornerHarris( grey, dst, blockSize, apertureSize, k, BORDER_DEFAULT );

  /// Normalizing
  normalize( dst, dst_norm, 0, 255, NORM_MINMAX, CV_32FC1, Mat() );
  convertScaleAbs( dst_norm, dst_norm_scaled );

  int thresh = 90;

  threshold(dst_norm_scaled, dst_norm_scaled, thresh, 0, THRESH_TOZERO);

  meanStdDev(dst_norm_scaled,mean,stddev);

  metrics["mean_corner"] = mean[0]/255.0;

  vector<vector<Point> > contours;

  findContours(dst_norm_scaled, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);

  metrics["n_corners"] = contours.size();

  // ** SEGMENTATION
  // cout <<"3\n";
  Mat blur;
//  GaussianBlur(image, blur, Size(15,15),1);

  // TODO: implement colour clustering

  resize(grey,blur,Size(),0.5,0.5);//downsample image by half for kmeans

  Mat pixels;
  blur.reshape(0, blur.total()).convertTo(pixels, CV_32F);
  Mat mini_labels;
  int n_clusters = 4;

  vector<float > centers;
  kmeans(pixels,n_clusters,mini_labels, TermCriteria(TermCriteria::EPS, 30, 0.1), 4, 0, centers);
  
  mini_labels = mini_labels.reshape(0, blur.rows);
  mini_labels.convertTo(mini_labels, CV_8U);
   
  // now apply cluster centers to full-size image
   
  Mat labels(grey.rows,grey.cols,grey.type());
  for(j = 0; j < grey.rows; j++){
    for(i = 0; i < grey.cols; i++){
      int p = grey.at<unsigned char>(j,i);
      
      int min_dist = 255;
      int q;
      
      // find cluster center nearest to this pixel
      for (q = 0; q < n_clusters; q++) {
        int dist = abs(centers[q] - p);
        if (dist<min_dist) {
          min_dist = dist;
          labels.at<unsigned char>(j,i) = centers[q];
        }
      }
    }
  }
  
  //imshow("labels",labels);
  //waitKey(0);

   
  // find connected components
  
  vector<patch > bigContours;
  // each item has format (patch_size/image_size, patch_mean). only contours with patch_size/image_size > 0.1 are added
  
//  double maxLabel;
//  minMaxLoc(labels, NULL, &maxLabel);
  
  Mat singleLabel;
  Mat mask(labels.rows,labels.cols,labels.type());
  
  double convexity = 0;
  
  // cout <<"4\n";
  vector<Vec4i> hierarchy;
  for (j = 0; j < n_clusters; j++) {
    inRange(labels,centers[j]-1,centers[j]+1,singleLabel);
    findContours(singleLabel, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_NONE);
    
    
    for (i = 0; i < (int) contours.size(); i++) {
      if (hierarchy[i][3] < 0) { // no parents i.e. not a hole
        double area = contourArea(contours[i]);
        
        int hole_idx = hierarchy[i][2];
        for( ; hole_idx >= 0; hole_idx = hierarchy[hole_idx][0] ) {
          // subtract hole area
          area -= contourArea(contours[hole_idx]);
        }
        
        double proportion = area/n_pix;
        if (proportion > 0.005) {
          // do convex hull stuff
          vector<Point> hull;
          convexHull(contours[i], hull);
          double hullArea = contourArea(hull);
          if (hullArea/area >= 0.8) {
            convexity += proportion;
          }
          
          if (proportion > 0.01) {
            mask = Scalar(0,0,0);
            drawContours(mask, contours, i, Scalar(255,255,255), -1, 4, hierarchy, 10);
            
            Scalar mean = cv::mean(image, mask);
            bigContours.push_back(patch(proportion, mean/255.0));
          }
        }
      }
    }
  }
  // cout <<"5\n";
  
  sort(bigContours.begin(), bigContours.end(), compareContours);
  
  metrics["n_patches"] = MIN(5, bigContours.size());
  
  // put stats for 5 biggest patches into metrics
  for (i = 0; i < 5; i++) {
    string mean_keystring = "patch_"+to_string(i+1)+"_mean";
    string size_keystring = "patch_"+to_string(i+1)+"_size";
    if (i < (int) bigContours.size()) {
      if (colour) {
        metrics[mean_keystring+"_h"] = bigContours[i].second[0];
        metrics[mean_keystring+"_l"] = bigContours[i].second[1];
        metrics[mean_keystring+"_s"] = bigContours[i].second[2];
      } else {
        metrics[mean_keystring+"_l"] = bigContours[i].second[0];
      }
      
      metrics[size_keystring] = bigContours[i].first;
      
      
    } else {
      if (colour) {
        metrics[mean_keystring+"_h"] = 0;
        metrics[mean_keystring+"_l"] = 0;
        metrics[mean_keystring+"_s"] = 0;
      } else {
        metrics[mean_keystring+"_l"] = 0;    
      }
      metrics[size_keystring] = 0; 
    }
  }
  
  metrics["convexity"] = convexity;
  
  int num_metrics = colour ? num_metrics_colour : num_metrics_grey;
  
  metrics_t metrics_vec;
  
  for (i = 0; i < num_metrics; i++) {
    if (colour) {
      metrics_vec.push_back(metrics[metric_names_colour[i]]);
    } else {
      metrics_vec.push_back(metrics[metric_names_grey[i]]);
    }
  }
//  for( auto n : metrics ) {
//    std::cout << "Key:[" << n.first << "] Value:[" << n.second << "]\n";
//  }
//  waitKey(0);
  return metrics_vec;
}


