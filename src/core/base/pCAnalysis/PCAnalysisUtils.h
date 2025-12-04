#pragma once

#include "Ponca/src/Fitting/covariancePlaneFit.h"
#include "Ponca/src/Fitting/curvature.h"
#include "Ponca/src/Fitting/meanPlaneFit.h"
#include <RIMLS.h>
#include <Eigen/Eigen>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/Map.h>

#include <Eigen/src/Core/Stride.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Ponca/SpatialPartitioning>
#include <Ponca/Fitting>
#include <iostream>

using namespace Ponca;

namespace  ttk::pointCloud{

double computeNDDistanceSquared(const float* coordsi, const float* coordsj, int dim);
double computeNDDistance(const float* coordsi, const float* coordsj, int dim);


// convenient helper to pass either computeNDDistanceLocalBasis or this
// function as argument
double computeNDDistanceHelper(const double* fitCoordsi,
                         const double* fitCoordsj,
                         const float* mlsCoordsi,
                         const float* mlsCoordsj,
                         const float* globalBasisCenter,
                         int dim);

double computeNDDistanceLocalBasis(const double* fitCoordsi,
                                   const double* fitCoordsj,
                                   const float* mlsCoordsi,
                                   const float* mlsCoordsj,
                                   const float* globalBasisCenter,
                                   int dim);

double compute5DDistanceSquared(const float* coordsi, const float* coordsj);
double compute5DDistance(const float* coordsi, const float* coordsj);
double compute5DDistanceHelper(const double* fitCoordsi,
                         const double* fitCoordsj,
                         const float* mlsCoordsi,
                         const float* mlsCoordsj,
                         const float* globalBasisCenter);

double compute5DDistanceLocalBasis(const double* fitCoordsi,
                                   const double* fitCoordsj,
                                   const float* mlsCoordsi,
                                   const float* mlsCoordsj,
                                   const float* globalBasisCenter);

void changeBasis(const double* inputFitCoords,
                 const float* diffCenter,
                 double* newFitCoords);

  struct SamplingGrid{
    int dimX; int dimY; int dimZ;
    double dX; double dY; double dZ;
    double originX; double originY; double originZ;

    inline size_t size() const {
      return this->dimX * this->dimY * this->dimZ;
    }
  };

template <int _Dim>
  class MyPointND {

public:
  enum {Dim = _Dim};
  using Scalar = double;
  using VectorType = Eigen::Matrix<Scalar, Dim, 1>;
  using MatrixType = Eigen::Matrix<Scalar, Dim, Dim>;


    PONCA_MULTIARCH inline MyPointND(Scalar* _coordinatesArray,  int _pId)
      : m_pos   (Eigen::Map< const VectorType >(_coordinatesArray + Dim*_pId  )),
      m_normal   (nullptr)
    {}

  PONCA_MULTIARCH inline const Eigen::Map< const VectorType >& pos()    const { return m_pos; }
  PONCA_MULTIARCH inline const Eigen::Map< const VectorType >& normal() const { return m_normal; }

  void print(){
    std::cout<<" p "<<m_pos<<std::endl;
  }

private:
    Eigen::Map< const VectorType > m_pos, m_normal;

};


  class MyPoint5D {

public:
  enum {Dim = 5};
  using Scalar = double;
using VectorType = Eigen::Matrix<Scalar, Dim, 1>;
  using MatrixType = Eigen::Matrix<Scalar, Dim, Dim>;


    PONCA_MULTIARCH inline MyPoint5D(Scalar* _coordinatesArray,  int _pId)
      : m_pos   (Eigen::Map< const VectorType >(_coordinatesArray + Dim*_pId  )),
      m_normal   (nullptr)
    {}

  PONCA_MULTIARCH inline const Eigen::Map< const VectorType >& pos()    const { return m_pos; }
  PONCA_MULTIARCH inline const Eigen::Map< const VectorType >& normal() const { return m_normal; }

  void print(){
    std::cout<<" p "<<m_pos<<std::endl;
  }

private:
    Eigen::Map< const VectorType > m_pos, m_normal;

};

  class MyPoint6D {

public:
  enum {Dim = 6};
  using Scalar = float;
  using VectorType = Eigen::Matrix<Scalar, Dim, 1>;
  using MatrixType = Eigen::Matrix<Scalar, Dim/2, 2>;


    PONCA_MULTIARCH inline MyPoint6D(Scalar* _coordinatesArray, Scalar* _normalArray, int _pId, float r)
      : m_pos  (Eigen::Map< const MatrixType, 0, Eigen::OuterStride<>>(_coordinatesArray + 3 * _pId, 3, 2,  Eigen::OuterStride<>((long unsigned int)(_normalArray-_coordinatesArray))).reshaped().cwiseProduct(Eigen::Matrix<Scalar, Dim, 1>{1/(1+r),1/(1+r),1/(1+r),r/(1+r),r/(1+r),r/(1+r)})),
        m_normal (nullptr)
    {}

  PONCA_MULTIARCH inline const VectorType& pos()    const { return m_pos; }
  PONCA_MULTIARCH inline const Eigen::Map< const VectorType >& normal() const { return m_normal; }

  void print(){
    std::cout<<" p "<<m_pos<<std::endl;
  }

private:
    const VectorType m_pos;
    Eigen::Map< const VectorType> m_normal;

};

  class MyPoint {

public:
  enum {Dim = 3};
  using Scalar = float;
  using VectorType = Eigen::Matrix<Scalar, Dim, 1>;
  using MatrixType = Eigen::Matrix<Scalar, Dim, Dim>;


    PONCA_MULTIARCH inline MyPoint(Scalar* _coordinatesArray, Scalar* _normalArray, int _pId)
      : m_pos   (Eigen::Map< const VectorType >(_coordinatesArray + 3 * _pId  )),
      m_normal   (Eigen::Map< const VectorType >(_normalArray + 3 * _pId  ))
  {}

    PONCA_MULTIARCH inline MyPoint(Scalar* _coordinatesArray, Scalar* _normalArray, int _pId, float r)
      : m_pos   (Eigen::Map< const VectorType >(_coordinatesArray + 3 * _pId  )),
      m_normal   (Eigen::Map< const VectorType >(_normalArray + 3 * _pId  ))
  {}
    PONCA_MULTIARCH inline MyPoint(Scalar* _coordinatesArray, int _pId)
      : m_pos   (Eigen::Map< const VectorType >(_coordinatesArray + 3*_pId  )),
      m_normal   (nullptr)
  {}

  PONCA_MULTIARCH inline const Eigen::Map< const VectorType >& pos()    const { return m_pos; }
  PONCA_MULTIARCH inline const Eigen::Map< const VectorType >& normal() const { return m_normal; }

  void print(){
    std::cout<<" p "<<m_pos<<std::endl;
  }


private:
    Eigen::Map< const VectorType > m_pos, m_normal;

};

using Scalar = MyPoint::Scalar;
using WeightFunc = DistWeightFunc<MyPoint, Ponca::SmoothWeightKernel<Scalar>>;
using WeightFuncRIMLS = RIMLSWeightFunc<MyPointRIMLS<3>, Ponca::SmoothWeightKernel<Scalar>>;
using PointRIMLS = MyPointRIMLS<3>;

using CovPlaneFit = Basket<MyPoint, WeightFunc, Ponca::CovariancePlaneFit>;
using mySphereFit = Basket<MyPoint, WeightFunc, Ponca::SphereFit, Ponca::GLSParam>;
using myOrientedSphereFitGLS = Basket<MyPoint, WeightFunc, Ponca::OrientedSphereFit, Ponca::GLSParam>;
using geomVarFit = BasketDiff<myOrientedSphereFitGLS,
                            Ponca::FitScaleSpaceDer, Ponca::OrientedSphereDer, Ponca::GLSDer /*, Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator*/>;

using MeanPlaneFit = Basket<MyPoint, WeightFunc, Ponca::MeanPlaneFit>;
using MeanPlaneFitRIMLS = Basket<PointRIMLS, WeightFuncRIMLS, Ponca::MeanPlaneFit>;

using myOrientedSphereFit = Basket<MyPoint, WeightFunc, Ponca::OrientedSphereFit>;
using FitASODiff = Ponca::BasketDiff<
        myOrientedSphereFit,
        Ponca::DiffType::FitSpaceDer,
        Ponca::OrientedSphereDer, Ponca::MlsSphereFitDer,
        Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;

using myOrientedSphereFitRIMLS = Basket<PointRIMLS, WeightFuncRIMLS, Ponca::OrientedSphereFit>;
using FitASODiffRIMLS = Ponca::BasketDiff<
        myOrientedSphereFitRIMLS,
        Ponca::DiffType::FitSpaceDer,
        Ponca::OrientedSphereDer, Ponca::MlsSphereFitDer,
        Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;


using myOrientedEllipsoidFit = Basket<MyPoint, WeightFunc, Ponca::OrientedEllipsoidFit>;
using FitEllipsoid = Ponca::BasketDiff<
        myOrientedEllipsoidFit,
        Ponca::DiffType::FitSpaceDer,
        Ponca::CurvatureEstimatorBase,
        Ponca::NormalDerivativesCurvatureEstimator>;

using myOrientedEllipsoidFitRIMLS = Basket<PointRIMLS, WeightFuncRIMLS, Ponca::OrientedEllipsoidFit>;
using FitEllipsoidRIMLS = Ponca::BasketDiff<
        myOrientedEllipsoidFitRIMLS,
        Ponca::DiffType::FitSpaceDer,
        Ponca::CurvatureEstimatorBase,
        Ponca::NormalDerivativesCurvatureEstimator>;



}
