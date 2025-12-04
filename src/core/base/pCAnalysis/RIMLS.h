#pragma once

// ttk common includes
// #include "Ponca/src/SpatialPartitioning/KdTree/kdTree.h"
// #include <PCAnalysisUtils.h>
#include <RIMLSWeightFunc.h>
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

    template <int _Dim>
      class MyPointRIMLS {

        public:
          enum {Dim = 3};
          using Scalar = float;
          using VectorType = Eigen::Matrix<Scalar, Dim, 1>;
          using MatrixType = Eigen::Matrix<Scalar, Dim, Dim>;


          PONCA_MULTIARCH inline MyPointRIMLS(const Scalar* _coordinatesArray, const Scalar* _normalArray, int _pId)
            : m_pos   (Eigen::Map< const VectorType >(_coordinatesArray + 3 * _pId  )),
            m_normal   (Eigen::Map< const VectorType >(_normalArray + 3 * _pId  )),
            m_gradientDiff(0),
            m_potentialDiff(0)
        {}

          PONCA_MULTIARCH inline MyPointRIMLS(const Scalar* _coordinatesArray, const Scalar* _normalArray, int _pId, float r)
            : m_pos   (Eigen::Map< const VectorType >(_coordinatesArray + 3 * _pId  )),
            m_normal   (Eigen::Map< const VectorType >(_normalArray + 3 * _pId  )),
            m_gradientDiff(0),
            m_potentialDiff(0)
        {}
          PONCA_MULTIARCH inline MyPointRIMLS(const Scalar* _coordinatesArray, int _pId)
            : m_pos   (Eigen::Map< const VectorType >(_coordinatesArray + 3*_pId  )),
            m_normal   (nullptr),
            m_gradientDiff(0),
            m_potentialDiff(0)
        {}

          PONCA_MULTIARCH inline const Eigen::Map< const VectorType >& pos()    const { return m_pos; }
          PONCA_MULTIARCH inline const Eigen::Map< const VectorType >& normal() const { return m_normal; }
          inline Scalar gradientDiff()  const {return m_gradientDiff;}
          inline Scalar potentialDiff() const {return m_potentialDiff;}


          inline void setGradientDiff(Scalar data){
            m_gradientDiff = data;
          }
          inline void setPotentialDiff(Scalar data){
            m_potentialDiff = data;
          }

          void print(){
            std::cout<<" p "<<m_pos<<std::endl;
          }


        private:
          Eigen::Map< const VectorType > m_pos, m_normal;
          Scalar m_gradientDiff{0};
          Scalar m_potentialDiff{0};

      };


    template <class FitStep, class FitFinal>
      class RIMLS : virtual public Debug {

        public:
          enum {Dim = 3};

          using Scalar = float;
          using VectorType = Eigen::Matrix<Scalar, Dim, 1>;
          using MatrixType = Eigen::Matrix<Scalar, Dim, Dim>;
          using Vector3 = VectorType;

          using Point = MyPointRIMLS<Dim>;
          using WeightFuncRIMLS = RIMLSWeightFunc<MyPointRIMLS<Dim>, Ponca::SmoothWeightKernel<Scalar>>;

          // using OrientedSphereFitRIMLS = Basket<Point, WeightFuncRIMLS, FitFinal>;
          // using FitASODiffRIMLS = Ponca::BasketDiff<
          //   OrientedSphereFitRIMLS,
          //   Ponca::DiffType::FitSpaceDer,
          //   Ponca::OrientedSphereDer, Ponca::MlsSphereFitDer,
          //   Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;

          // using FitStep = FitASODiffRIMLS;
          // using FitFinal = FitASODiffRIMLS;

        public:
          RIMLS();
          RIMLS(const RIMLS& other);
          ~RIMLS();


        private: 

          Scalar m_scale;

          int         m_reweighting_step;

          int         m_step_max;
          Scalar      m_convergence_ratio_min;

          int         m_step;
          bool        m_stable;
          int         m_neighbor_count;

          WeightFuncRIMLS  m_weight_func;
          FitStep         m_fit_step;
          FitFinal         m_fit_final;

        public:

          template<typename CoordsType, typename AccelerationStruct>
            void compute(const CoordsType* coordinates, const CoordsType* normals, const Vector3& point_init, AccelerationStruct& tree, VectorType& point );

        public:
          void set_scale(Scalar scale);
          void set_step_max(int step_max);
          void set_convergence_ratio_min(Scalar convergence_ratio_min);

          void set_reweighting_step(int step);

          // Accessors ---------------------------------------------------------------
        public:
          bool  stable() const;
          const FitFinal& fit() const{
            return m_fit_final;
          }
          FitFinal& fit(){
            return m_fit_final;
          }

      };


    template<typename FitStep, typename FitFinal>
RIMLS<FitStep, FitFinal>::RIMLS() :
    m_scale(1.0),
    m_reweighting_step(1),
    m_step_max(100),
    m_convergence_ratio_min(0.01),
    m_step(0),
    m_stable(false),
    m_neighbor_count(0),
    m_weight_func(1.0),
    m_fit_step(),
    m_fit_final()
{
}

    template<typename FitStep, typename FitFinal>
RIMLS<FitStep, FitFinal>::~RIMLS()
{
}

    template<typename FitStep, typename FitFinal>
RIMLS<FitStep, FitFinal>::RIMLS(const RIMLS& other) :
    m_scale(other.m_scale),
    m_reweighting_step(other.m_reweighting_step),
    m_step_max(other.m_step_max),
    m_convergence_ratio_min(other.m_convergence_ratio_min),
    m_step(other.m_step),
    m_stable(other.m_stable),
    m_neighbor_count(other.m_neighbor_count),
    m_weight_func(other.m_weight_func),
    m_fit_step(other.m_fit_step),
    m_fit_final(other.m_fit_final)
{
}



// Accessors -------------------------------------------------------------------


    template<typename FitStep, typename FitFinal>
bool RIMLS<FitStep, FitFinal>::stable() const
{
    return m_stable;
}



// Parameters ------------------------------------------------------------------


  template<typename FitStep, typename FitFinal>
void RIMLS<FitStep, FitFinal>::set_scale(Scalar scale)
{
    m_weight_func = WeightFuncRIMLS(scale);
    m_scale = scale;
}


    template<typename FitStep, typename FitFinal>
void RIMLS<FitStep, FitFinal>::set_step_max(int step_max)
{
    m_step_max = step_max;
}


    template<typename FitStep, typename FitFinal>
void RIMLS<FitStep, FitFinal>::set_convergence_ratio_min(Scalar convergence_ratio_min)
{
    m_convergence_ratio_min = convergence_ratio_min;
}



    template<typename FitStep, typename FitFinal>
void RIMLS<FitStep, FitFinal>::set_reweighting_step(int step)
{
    m_reweighting_step= step;
}


    template<typename FitStep, typename FitFinal>
      template<typename CoordsType, typename AccelerationStruct>
      void RIMLS<FitStep, FitFinal>::compute(const CoordsType* coordinates, const CoordsType* normals, const Vector3& point_init, AccelerationStruct& tree, VectorType& point )
      {

        // First steps -------------------------------------------------------------
        m_fit_step.setWeightFunc(m_weight_func);
        Scalar dist_min = m_convergence_ratio_min * m_scale;
        Scalar dist     = 10 * dist_min; // arbitrary initialization > dist_min
        m_step          = 0;

        // Scalar  uc       = Scalar(0);
        // Vector3 ul       = Vector3::Zero();
        // Scalar  uq       = Scalar(0);

        Scalar  potential = Scalar(0);
        Vector3 gradient  = Vector3::Zero();

        m_stable        = true;
        bool converge   = dist < dist_min;
        bool reach_max  = m_step >= m_step_max-1;

        point = point_init;

        while(m_stable && !reach_max && !converge)
        {
          // Reweighting
          for(int n=0; n<m_reweighting_step; ++n)
          {
            // init
            m_fit_step.init(point);

            // add neighbors
            auto tree_query = tree.range_neighbors(point, m_scale);
            for(auto idx_nei : tree_query)
            {
              Scalar diffN = Scalar(0.);
              Scalar diffP = Scalar(0.);
              Point p = Point(coordinates, normals, idx_nei);
              if(n>0)
              {
                Vector3 q = p.pos() - point;
                Scalar  s = m_fit_step.potential(q);
                diffN = (gradient - p.normal()).norm();
                diffP = potential-s;
              }
              p.setGradientDiff(diffN);
              p.setPotentialDiff(diffP);
              m_fit_step.addNeighbor(p);
            }

            // finalize
            m_stable = (m_fit_step.finalize() == Ponca::STABLE);
            if(m_stable)
            {
              // uc = m_fit_step.m_uc;
              // ul = m_fit_step.m_ul;
              // uq = m_fit_step.m_uq;
              potential = m_fit_step.potential();
              gradient  = m_fit_step.primitiveGradient();
            }
          }

          // project
          if(m_stable)
          {
            Vector3 proj = m_fit_step.project(point);
            dist         = (proj-point).norm();
            converge     = dist < dist_min;
            point        = proj;
          }

          // increment
          ++m_step;
          reach_max = m_step >= m_step_max-1;
        }

        // Last step ---------------------------------------------------------------
        if(m_stable)
        {
          // init
          m_fit_final.setWeightFunc(m_weight_func);
          m_fit_final.init(point);
          m_neighbor_count = 0;

          // add neighbors
          auto tree_query = tree.range_neighbors(point, m_scale);
          for(auto idx_nei : tree_query)
          {
            Scalar diffN = Scalar(0.);
            Scalar diffP = Scalar(0.);
            Point p = Point(coordinates, normals, idx_nei);
            if(m_step > 0 && m_reweighting_step > 1)
            {
              Vector3 q = p.pos()- point;
              Scalar  s = m_fit_step.potential(q);
              diffN = (gradient - p.normal()).norm();
              diffP = potential-s;
            }
            p.setPotentialDiff(diffP);
            p.setGradientDiff(diffN);

            m_fit_final.addNeighbor(p);
            ++m_neighbor_count;
          }

          // finalize
          m_stable = (m_fit_final.finalize() == Ponca::STABLE);
          if(m_stable)
          {
            point = m_fit_final.project(point);
            // m_fit_final.computeCurvature(true); // useNormal = true
          }
          ++m_step;
        }
      }

  } //end pointcloud
} // end ttk
