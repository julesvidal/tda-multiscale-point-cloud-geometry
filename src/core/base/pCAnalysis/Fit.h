/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::PCAnalysis
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %PCAnalysis class that computes for each vertex of a
/// triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'PCAnalysis'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
// #include "Ponca/src/SpatialPartitioning/KdTree/kdTree.h"
#include <Sampling.h>
#include <PCAnalysisUtils.h>
#include <RIMLS.h>
#include <Debug.h>
#include <Triangulation.h>
#include <boost/mpl/size.hpp>
#include <cctype>
#include <cstddef>
#include <immintrin.h>
#include <iostream>
#include <limits>
#include <cmath>
#include <sstream>
#include <string>
#include <utility>

#include <cmath>

using namespace Ponca;
using namespace ttk;
using namespace pointCloud;

namespace ttk {

namespace pointCloud {


  class PCAnalysis : virtual public Debug {

  enum {Dim = 3};
  using Scalar = float;
  using VectorType = Eigen::Matrix<Scalar, Dim, 1>;
  using MatrixType = Eigen::Matrix<Scalar, Dim, Dim>;

  protected:
    int nbIterationsMLS_{5};
    ScaleSpaceSampling<MyPoint::Scalar, MyPoint> scaleSpaceSampler_{};

  public:
    PCAnalysis();


    inline ScaleSpaceSampling<MyPoint::Scalar, MyPoint>* sampler(){
      return &scaleSpaceSampler_;
    }


    template<typename Fit, typename AccelerationStruct>
    int test_fit_kdtree(double tmax, Fit& _fit, Scalar* coordinates, Scalar* normals, const VectorType& _p, AccelerationStruct& _tree)
    {

      // Set a weighting function instance
      _fit.setWeightFunc(WeightFunc(tmax));

      // Set the evaluation position
      _fit.init(_p);

      auto tree_query = _tree.range_neighbors(_p,tmax);

      for(auto it=tree_query.begin(); it != tree_query.end(); it++){
        _fit.addNeighbor(MyPoint(coordinates,normals,*it));
      }

      _fit.finalize();

      return _fit.isStable();
    }




    template<typename Fit, typename CoordsType, typename AccelerationStruct>
    int computeMLSKnnGraph(double _tmax, Fit& _fit, CoordsType* _coordinates, CoordsType* _normals, const VectorType _p0, const SimplexId _p0_index, AccelerationStruct& _tree, VectorType& _p_proj){

      // initial fit
      _fit.setWeightFunc(WeightFunc(_tmax));
      _fit.init(_p0);
      auto tree_query = _tree.range_neighbors(_p0_index, _tmax);

      for(auto it : tree_query){
        _fit.addNeighbor(MyPoint(_coordinates,_normals, it));
      }
      _fit.finalize();
      int stable = _fit.isStable();

      int i=1;
      int go_on = stable && (i < nbIterationsMLS_);

      VectorType p = _p0;


      // mls loop
      while(go_on){

        p = _fit.project(p);

        // fit point p
        _fit.setWeightFunc(WeightFunc(_tmax));
        _fit.init(_p0);

        for(auto it : tree_query){
          _fit.addNeighbor(MyPoint(_coordinates,_normals, it));
        }
        _fit.finalize();
        stable = _fit.isStable();
        i++;
        go_on = stable && (i<nbIterationsMLS_);
      } // end mls loop

      _p_proj=p;

      return stable;
    } 

    template<typename Fit, typename CoordsType, typename AccelerationStruct>
    int computeMLS(double _tmax, Fit& _fit, CoordsType* _coordinates, CoordsType* _normals, const VectorType _p0, AccelerationStruct& _tree, VectorType& _p_proj){

      // initial fit
      test_fit_kdtree(_tmax, _fit, _coordinates, _normals, _p0, _tree);
      int stable = _fit.isStable();

      int i=1;
      int go_on = stable && (i < nbIterationsMLS_);

      VectorType p = _p0;


      // mls loop
      while(go_on){

        p = _fit.project(p);

        // fit point p
        test_fit_kdtree(_tmax, _fit, _coordinates, _normals, p, _tree);
        stable = _fit.isStable();
        i++;
        go_on = stable && (i<nbIterationsMLS_);
      } // end mls loop

      _p_proj=p;

      return stable;
    } 



template<class Fit, int FitDim, class SamplerType>
bool computeFitCoordinates(float* coordinates, float* normals, SimplexId
    nbPoints, double* fitCoordinates, int* isStablePtr, float* mls3DCoordinates, float* mlsNormals, SamplerType&
    sampler, int iScale, std::vector<float>& resultBarycenter,
    std::vector<SimplexId>& sampling_indexes,
    bool useNormalization=true){

    printMsg("compute fit coordinates");


    Timer tm{};

    using VectorType = pointCloud::MyPoint::VectorType;


    bool useKnnGraph = sampler.isUseKnnGraph();
    auto& tree = sampler.kdTree(iScale);
    auto& knnGraph = sampler.knnGraph();
    auto scale = sampler.scale(iScale);

    VectorType barycenter;
    std::cout<<"Basis set to Full barycenter"<<std::endl;

    // find barycenter
    double tm_bary=tm.getElapsedTime();
    float bx=0, by=0, bz=0;
    for(SimplexId I = 0; I < sampling_indexes.size(); I++){
      SimplexId iPoint = sampling_indexes[I];
      bx+=coordinates[3*iPoint];
      by+=coordinates[3*iPoint+1];
      bz+=coordinates[3*iPoint+2];
    }
    barycenter << bx/nbPoints, by/nbPoints, bz/nbPoints;
    printMsg("computed barycenter", 1, tm.getElapsedTime()-tm_bary);
    std::cout<<" bary "<<barycenter.transpose()<<std::endl;
    resultBarycenter.resize(3);
    resultBarycenter[0]=bx/nbPoints;
    resultBarycenter[1]=by/nbPoints;
    resultBarycenter[2]=bz/nbPoints;


    double tm_fitCoords=tm.getElapsedTime();
    unsigned char unstable = 0;
    std::srand(std::time(0));
    // for each point, compute the fit in the barycenter basis
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId I = 0; I < sampling_indexes.size(); I++){
      
      SimplexId iPoint = sampling_indexes[I];


      VectorType p, p_proj;
      p << coordinates[3*iPoint], coordinates[3*iPoint+1], coordinates[3*iPoint+2];
      Fit pointFit;

      bool isStable{false}; 

        if(useKnnGraph){
          isStable = computeMLSKnnGraph(scale, pointFit, coordinates, normals, p, iPoint, knnGraph, p_proj);
        }else{
          isStable = computeMLS(scale, pointFit, coordinates, normals, p, tree, p_proj);
        }

      for(int iFitDim=0; iFitDim<FitDim; iFitDim++){
        fitCoordinates[FitDim*I+iFitDim]=0;
      }

      if(!isStable){
        unstable = 255;
        isStablePtr[I] = 0;
      }else{

        isStablePtr[I]=1;


        VectorType mlsNormal;
        mlsNormal = pointFit.primitiveGradient();

        for(int iDim=0; iDim<3; iDim++){
          mls3DCoordinates[3*I+iDim]=p_proj.coeff(iDim);
          mlsNormals[3*I+iDim]=mlsNormal.coeff(iDim);
        }


        pointFit.changeBasis(barycenter);


        MyPoint::VectorType ul;
        MyPoint::Scalar uc;
        MyPoint::Scalar uq;

        uc=pointFit.m_uc;
        ul=pointFit.m_ul;
        uq = pointFit.m_uq;

        // uc
        fitCoordinates[FitDim*I]=uc;

        // ul
        fitCoordinates[FitDim*I+1]=ul.coeff(0);
        fitCoordinates[FitDim*I+2]=ul.coeff(1);
        fitCoordinates[FitDim*I+3]=ul.coeff(2);

        // uq
        fitCoordinates[FitDim*I+4]=uq;


        // normalize;
        if(useNormalization){
          MyPoint::Scalar norm=0;
          for(int k=0; k<FitDim; k++){
            norm += fitCoordinates[FitDim*I+k]*fitCoordinates[FitDim*I+k];
          }
          norm = sqrt(norm);
          for(int k=0; k<FitDim; k++){
            fitCoordinates[FitDim*I+k]/=norm;
          }
        }
      }
    }
    if(unstable){
        printWrn("At least one fit was not stable, expect segfault");
    }

    printMsg("Computed fit "+std::to_string(FitDim)+"D coordinates with scale "+std::to_string(scale), 1, tm.getElapsedTime()-tm_fitCoords);


    return true;
}


  }; // PCAnalysis class
  

} // namespace pointCloud
} // namespace ttk
