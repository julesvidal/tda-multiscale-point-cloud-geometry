#include "PCAnalysisUtils.h"
#include <PointCloudNoiser.h>
#include <PCAnalysis.h>
#include <boost/mpl/size.hpp>
#include <string>

using namespace ttk;

ttk::PointCloudNoiser::PointCloudNoiser() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("PointCloudNoiser");
}

int PointCloudNoiser::addGaussianNoise(const double standard_deviation,
                                      const double mean,
                                      float* coordinates,
                                      const SimplexId n) const {
  int dim=3;
    // std::random_device rd{};
    std::mt19937 gen{0};
    std::normal_distribution d{mean, standard_deviation};

    for(SimplexId i=0; i<n; i++){
      for(int iDim=0; iDim<dim; iDim++){
      coordinates[dim*i+iDim] += (float)d(gen);
      }
    }
    return 1;
}

int PointCloudNoiser::SmoothNormalsWithFit(double _scale, float* _coordinates, float* _normals, int _nbPoints, float* _smoothed_normals){
  
  pointCloud::PCAnalysis pca{};
  pca.setThreadNumber(threadNumber_);
  pca.setDebugLevel(debugLevel_);

  printMsg("Smoothing normals with scale "+std::to_string(_scale));
      
  pca.fitNormals(_scale, _coordinates, _normals, _nbPoints, _smoothed_normals);

  // std::cout<<" DOne \n\n"<<std::endl;
  // for(int i=0; i<_nbPoints; i++){
  //       std::cout<<"Point "<<i<<" , coords "<<_coordinates[3*i]<<"  "<<_coordinates[3*i+1]<<"  "<<_coordinates[3*i+2]<<"  "<<", normals "<<_smoothed_normals[3*i]<<"  "<<_smoothed_normals[3*i+1]<<"  "<<_smoothed_normals[3*i+2]<<"  "<<std::endl;
  // }
  return 1;
}

