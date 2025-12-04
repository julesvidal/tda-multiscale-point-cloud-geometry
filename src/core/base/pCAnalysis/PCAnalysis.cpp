// #include <Eigen/Eigen>
// #include <Eigen/src/Core/Matrix.h>
#include "PCAnalysisUtils.h"
#include <Eigen/Core>
#include <PCAnalysis.h>
#include <Ponca/Fitting>
#include <limits>

using namespace ttk;
using namespace Ponca;
using namespace pointCloud;

PCAnalysis::PCAnalysis() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("PCAnalysis");
}

// typedef MyPoint::Scalar Scalar;
// typedef DistWeightFunc<MyPoint,SmoothWeightKernel<Scalar> > WeightFunc;
// typedef Basket<MyPoint, WeightFunc, CovariancePlaneFit> CovPlaneFit;
// typedef Basket<MyPoint, WeightFunc, SphereFit> mySphereFit;


namespace ttk {

  namespace pointCloud {

      int Primitive::projectionDistanceOnNormal(const float* coords,
                                               const float* normal,
                                               Scalar& dist) const {
      VectorType q;
      q << coords[0], coords[1], coords[2];

      VectorType n;
      n << normal[0], normal[1], normal[2];

      return projectionDistanceOnNormal(q, n, dist);


      }
      int Primitive::projectionDistanceOnNormal(const VectorType& q,
                                               const VectorType& n,
                                               Scalar& dist) const {
      const VectorType lq = toLocalBasis(q);

      constexpr Scalar epsilon = Eigen::NumTraits<Scalar>::dummy_precision();


      Scalar a,b,c;
      Scalar potential, norm;
      VectorType grad;

      Scalar t=0;

      int ret = 1;


      switch(m_type){
        case PrimitiveType::PLANE:
          c = m_uc + lq.dot(m_ul);
          b = m_ul.dot(n);

          if(b>epsilon){
            t = - c/b;
          }else{
            t = std::numeric_limits<Scalar>::max();
            ret = 0;
          }
          dist = std::abs(t*n.norm());
          break;

        case PrimitiveType::SPHERE:
          // grad = n + Scalar(2.f) * m_uq * lq;
          a = m_uq * n.norm();
          b = m_ul.dot(n) + Scalar(2)*m_uq*lq.dot(n);
          c = m_uc + lq.dot(m_ul) + m_uq * lq.squaredNorm();
          // norm = grad.norm();

          if(abs(a)<epsilon){
            if(abs(b)<epsilon){
              if(abs(c)<epsilon){
                t = 0;
              }else{
                t = std::numeric_limits<Primitive::Scalar>::max();
                ret = 0;
              }
            }else{
              t = -c/b;
            }
          }else{
            Scalar delta = b*b - 4*a*c;
            if(delta >= 0){
              t = (- b + sqrt(delta)) / (2 * a);
            }else{ // no solution
              t = std::numeric_limits<Primitive::Scalar>::max();
              ret = 0;
            }
          } 
          dist = std::abs(t*n.norm());
          break;

        case PrimitiveType::QUADRIC:

          a = n.dot(m_Uq * n);
          b = m_ul.dot(n) + n.dot(m_Uq*lq) + lq.dot(m_Uq*n);
          c = m_uc + lq.dot(m_ul) + lq.dot(m_Uq * lq);

          if(abs(a)<epsilon){
            if(abs(b)<epsilon){
              if(abs(c)<epsilon){
                t = 0;
              }else{
                t = std::numeric_limits<Primitive::Scalar>::max();
                ret = 0;
              }
            }else{
              t = -c/b;
            }
          }else{
            Scalar delta = b*b - 4*a*c;
            if(delta >= 0){
              t = (- b + sqrt(delta)) / (2 * a);
            }else{ // no solution
              t = std::numeric_limits<Primitive::Scalar>::max();
              ret = 0;
            }
          } 
          dist = std::abs(t*n.norm());
          break;
      }
      return ret;
    }

      int Primitive::projectionDistanceAndNormal(const VectorType& q, Scalar& dist, VectorType& normal) const {
      const VectorType lq = toLocalBasis(q);

      constexpr Scalar epsilon = Eigen::NumTraits<Scalar>::dummy_precision();


      Scalar a,b,c;
      Scalar potential, norm;
      VectorType grad;

      Scalar t=0;

      int ret = 1;


      switch(m_type){
        case PrimitiveType::PLANE:
          grad = m_ul;
          potential = m_uc + lq.dot(m_ul);
          norm = grad.norm();
          b = norm*norm;

          if(b>epsilon){
            t = - potential/norm;
          }
          dist = std::abs(t);
          normal = grad.normalized();
          break;

        case PrimitiveType::SPHERE:
          grad = m_ul + Scalar(2.f) * m_uq * lq;
          potential = m_uc + lq.dot(m_ul) + m_uq * lq.squaredNorm();
          norm = grad.norm();

          t =  (- norm + sqrt(norm*norm - Scalar(4) * m_uq * potential)) / (Scalar(2) * m_uq);
          dist = std::abs(t);
          normal = grad.normalized();
          break;

        case PrimitiveType::QUADRIC:

          grad = m_ul + Scalar(2.f) * m_Uq * lq;
          a = grad.transpose()*m_Uq * grad;
          b = grad.squaredNorm();
          c = m_uc + lq.dot(m_ul) + lq.transpose()* m_Uq * lq;
          norm = grad.norm();


          if(abs(a)<epsilon){
            if(abs(b)<epsilon){
              t = 0;
            }else{
              t = -c/b;
            }
          }else{
            Scalar delta = b*b - 4*a*c;
            if(delta >= 0){
              t = (- b + sqrt(delta)) / (2 * a);
            }else{ // no solution
              t = std::numeric_limits<Primitive::Scalar>::max();
              ret = 0;
            }
          } 
          dist = std::abs(t*norm);
          normal = grad.normalized();
          break;
      }
      return ret;
    }

    Primitive::Scalar Primitive::projectionDistance(const VectorType& q) const {
      const VectorType lq = toLocalBasis(q);

      constexpr Scalar epsilon = Eigen::NumTraits<Scalar>::dummy_precision();


      Scalar a,b,c;
      Scalar potential, norm;
      VectorType grad;

      Scalar t=0;


      switch(m_type){
        case PrimitiveType::PLANE:
          grad = m_ul;
          potential = m_uc + lq.dot(m_ul);
          norm = grad.norm();
          b = norm*norm;

          if(b>epsilon){
            t = - potential/norm;
          }
          return std::abs(t);
          break;

        case PrimitiveType::SPHERE:
          grad = m_ul + Scalar(2.f) * m_uq * lq;
          potential = m_uc + lq.dot(m_ul) + m_uq * lq.squaredNorm();
          norm = grad.norm();

          t =  (- norm + sqrt(norm*norm - Scalar(4) * m_uq * potential)) / (Scalar(2) * m_uq);
          return std::abs(t);
          break;

        case PrimitiveType::QUADRIC:

          grad = m_ul + Scalar(2.f) * m_Uq * lq;
          a = grad.transpose()*m_Uq * grad;
          b = grad.squaredNorm();
          c = m_uc + lq.dot(m_ul) + lq.transpose()* m_Uq * lq;
          norm = grad.norm();


          if(abs(a)<epsilon){
            if(abs(b)<epsilon){
              t = 0;
            }else{
              t = -c/b;
            }
          }else{
            Scalar delta = b*b - 4*a*c;
            if(delta >= 0){
              t = (- b + sqrt(delta)) / (2 * a);
              // std::cout<<"t "<<t<<std::endl;
              // std::cout<<"t*norm "<<t*norm<<std::endl;
            }else{ // no solution
              t = std::numeric_limits<Primitive::Scalar>::max();
            }
          } 
          return std::abs(t*norm);
          break;
      }
    }
    void Primitive::project(const VectorType& q, VectorType& q_proj, VectorType& n_proj) const {

      const VectorType lq = toLocalBasis(q);

      constexpr Scalar epsilon = Eigen::NumTraits<Scalar>::dummy_precision();


      Scalar a,b,c;
      Scalar potential, norm;
      VectorType grad;

      Scalar t=0;

      switch(m_type){
        case PrimitiveType::PLANE:
          grad = m_ul;
          potential = m_uc + lq.dot(m_ul);
          norm = grad.norm();
          b = norm*norm;

          if(b>epsilon){
            t = - potential/(b);
          }
          break;

        case PrimitiveType::SPHERE:
          grad = m_ul + Scalar(2.f) * m_uq * lq;
          potential = m_uc + lq.dot(m_ul) + m_uq * lq.squaredNorm();
          norm = grad.norm();

          if(norm>epsilon){
            t =  (- norm + sqrt(norm*norm - Scalar(4) * m_uq * potential)) / (Scalar(2) * m_uq * norm);
          }
          break;
        case PrimitiveType::QUADRIC:
          grad = m_ul + Scalar(2.f) * m_Uq * lq;
          c = m_uc + lq.dot(m_ul) + lq.transpose()* m_Uq * lq;
          Scalar a = grad.transpose()*m_Uq * grad;
          Scalar b = grad.squaredNorm();
          Scalar delta = b*b - 4*a*c;
          if(a>epsilon and delta >= 0)
          {
            t = (- b + sqrt(delta)) / (2 * a);
          }
          break;
      }
      q_proj = toGlobalBasis(lq + t*grad);
      n_proj = grad.normalized();
    }


    int Primitive::projectionDistanceAndNormal(const Scalar* coords, Scalar& dist, VectorType& normal) const {
      VectorType q;
      q << coords[0], coords[1], coords[2];
      return projectionDistanceAndNormal(q, dist, normal);
    }

    Primitive::Scalar Primitive::projectionDistance(const Scalar* coords) const {
      VectorType q;
      q << coords[0], coords[1], coords[2];
      return projectionDistance(q);
    }


    double computeNDDistanceSquared(const double* coordsi, const double* coordsj, int dim){
      double dist2 = 0;
      for(int k=0; k<dim; k++){
        dist2 += (coordsi[k]-coordsj[k])*(coordsi[k]-coordsj[k]);
      }
      return dist2;
    }

    double computeNDDistance(const double* coordsi, const double* coordsj, int dim){
      return std::sqrt(computeNDDistanceSquared(coordsi, coordsj, dim));
    }


    double computeNDDistanceLocalBasis(const double* fitCoordsi,
        const double* fitCoordsj,
        const float* mlsCoordsi,
        const float* mlsCoordsj,
        const float* globalBasisCenter,
        int dim){
      if(dim!=5){
        std::cout<<"/!\\ Local metrics only available in R5"<<std::endl;
      }
      return compute5DDistanceLocalBasis(fitCoordsi, fitCoordsj, mlsCoordsi, mlsCoordsj, globalBasisCenter);
    }
    double computeNDDistanceHelper(const double* fitCoordsi,
        const double* fitCoordsj,
        const float* mlsCoordsi,
        const float* mlsCoordsj,
        const float* globalBasisCenter, int dim){
      return computeNDDistance(fitCoordsi, fitCoordsj, dim);
    }



    double compute5DDistanceSquared(const double* coordsi, const double* coordsj){
      double dist2 = 0;
      for(int k=0; k<5; k++){
        dist2 += (coordsi[k]-coordsj[k])*(coordsi[k]-coordsj[k]);
      }
      return dist2;
    }

    double compute5DDistance(const double* coordsi, const double* coordsj){
      return std::sqrt(compute5DDistanceSquared(coordsi, coordsj));
    }


    double compute5DDistanceLocalBasis(const double* fitCoordsi,
        const double* fitCoordsj,
        const float* mlsCoordsi,
        const float* mlsCoordsj,
        const float* globalBasisCenter){
      float diff[3];
      double newFitCoordsi[5];
      double newFitCoordsj[5];
      for(int k=0; k<3; k++){
        diff[k] = globalBasisCenter[k] - (mlsCoordsi[k]+mlsCoordsj[k])/2;
        // std::cout<<"g "<<globalBasisCenter[k]<<" "<<mlsCoordsi[k]<<" "<<mlsCoordsj[k]<<" "<<diff[k]<<std::endl;
      }

      changeBasis(fitCoordsi, diff, newFitCoordsi);
      changeBasis(fitCoordsj, diff, newFitCoordsj);

      return compute5DDistance(newFitCoordsi, newFitCoordsj);
    }
    double compute5DDistanceHelper(const double* fitCoordsi,
        const double* fitCoordsj,
        const float* mlsCoordsi,
        const float* mlsCoordsj,
        const float* globalBasisCenter){
      return compute5DDistance(fitCoordsi, fitCoordsj);
    }

    void changeBasis(const double* inputFitCoords,
        const float* diffCenter,
        double* newFitCoords){
      // float diff[3];
      // for(int k=0; k<3; k++){
      //   diff[k] = inputBasisCenter[k]-newBasisCenter[k];
      // }

      float diffCenter2 = diffCenter[0]*diffCenter[0] 
        + diffCenter[1]*diffCenter[1] 
        + diffCenter[2]*diffCenter[2];

      float ul_dot_diffCenter = inputFitCoords[1]*diffCenter[0] 
        + inputFitCoords[2]*diffCenter[1] 
        + inputFitCoords[3]*diffCenter[2];

      // new uc
      newFitCoords[0] = inputFitCoords[0] /*uc*/
        -  ul_dot_diffCenter
        + inputFitCoords[1] * diffCenter2;

      // std::cout<<"new uc "<<newFitCoords[0]<<" "<<inputFitCoords[0]<<" "<<ul_dot_diffCenter<<" "<<inputFitCoords[1]<<" "<<diffCenter2<<std::endl;
      // new uq
      newFitCoords[1] = inputFitCoords[1];

      // new ul
      for(int k=0; k<3; k++){
        newFitCoords[1+k] = inputFitCoords[1+k] - 2*inputFitCoords[4]*diffCenter[k];
      }
      // std::cout<<" diff ";
      // for(int i=0; i<3; i++){
      //   std::cout<<" "<<diffCenter[i];
      // }
      // std::cout<<std::endl;
      // std::cout<<" fit";
      // for(int i=0; i<5; i++){
      //   std::cout<<" "<<inputFitCoords[i];
      // }
      // std::cout<<" changed to ";
      // for(int i=0; i<5; i++){
      //   std::cout<<" "<<newFitCoords[i];
      // }
      // std::cout<<std::endl;
    }
  }
}
