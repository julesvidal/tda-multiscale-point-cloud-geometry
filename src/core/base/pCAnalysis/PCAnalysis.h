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
#include "DataTypes.h"
#include <Eigen/Core>
#include <Sampling.h>
#include <PCAnalysisUtils.h>
#include <RIMLS.h>
#include <Debug.h>
#include <Triangulation.h>
#include <boost/mp11/list.hpp>
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


  class Primitive : public Debug {
    public:
      enum class PrimitiveType{PLANE=0, SPHERE=1, QUADRIC=2};
      enum {Dim=3};
      using Scalar = float;
      using VectorType = Eigen::Matrix<Scalar, Dim, 1>;
      using MatrixType = Eigen::Matrix<Scalar, Dim, Dim>;

    private:

      PrimitiveType m_type{PrimitiveType::SPHERE};
      Scalar m_uc{};
      VectorType m_ul{};
      Scalar m_uq{};
      MatrixType m_Uq{};

      std::vector<SimplexId> m_point_indices{};
      VectorType m_basisCenter{};

    public:
      Primitive();

      template <typename T>
      Primitive(T* params, PrimitiveType type, const VectorType& basisCenter): m_basisCenter{basisCenter} {

        constexpr Scalar epsilon = Eigen::NumTraits<Scalar>::dummy_precision();

        m_type = type;
        m_uc = Scalar(params[0]);
        m_ul << Scalar(params[1]), Scalar(params[2]), Scalar(params[3]);
        if(type==PrimitiveType::SPHERE){
          m_uq = Scalar(params[4]);
          if(std::abs(m_uq) < epsilon ){
            m_type = PrimitiveType::PLANE;
          }
        }else if(type==PrimitiveType::QUADRIC){
          m_Uq << params[4], params[7], params[8],
                  params[7], params[5], params[9],
                  params[8], params[9], params[6];
          if(m_Uq.isApprox(MatrixType::Zero())){
            m_type = PrimitiveType::PLANE;
          }
        }
      }

      inline Scalar f(float* point) const {
        VectorType p;
        p << point[0], point[1], point[2];
        VectorType lp = toLocalBasis(p);
        switch(m_type){
          case PrimitiveType::PLANE:
            return m_uc + m_ul.dot(lp);
            break;
          case PrimitiveType::SPHERE:
            return m_uc + m_ul.dot(lp) + m_uq*lp.squaredNorm();
            break;
          case PrimitiveType::QUADRIC:
            return m_uc + m_ul.dot(lp) + lp.dot(m_Uq*lp);
            break;
          default:
            return 0;
        }
      }

      inline Scalar kappa() const {
        switch(m_type){
          case PrimitiveType::PLANE:
            return 0;
            break;
          case PrimitiveType::SPHERE:
            return 2*std::abs(m_uq);
            break;
          case PrimitiveType::QUADRIC:
            return 2*std::abs(m_Uq.col(0)[0]+m_Uq.col(1)[1]+m_Uq.col(2)[2])/3;
            break;
          default:
            return 0;
        }
      }
      inline std::string type_str() const {
        switch(m_type){
          case PrimitiveType::PLANE:
            return "plane";
            break;
          case PrimitiveType::SPHERE:
            return "sphere";
            break;
          case PrimitiveType::QUADRIC:
            return "quadric";
            break;
          default:
            return "";
        }
      }
      inline bool isQuadric() const {
        return m_type==PrimitiveType::QUADRIC;
      }
      inline bool isSphere() const {
        return m_type==PrimitiveType::SPHERE;
      }
      inline bool isPlane() const {
        return m_type==PrimitiveType::PLANE;
      }

      inline VectorType toLocalBasis(const VectorType& p) const {
        return p-m_basisCenter;
      }
      inline VectorType toGlobalBasis(const VectorType& p) const {
        return p+m_basisCenter;
      }
      void addPoint(SimplexId id){
        m_point_indices.push_back(id);
      }
      SimplexId getNumberOfPoints() const {
        return m_point_indices.size();
      }
      void project(const VectorType& q, VectorType& q_proj, VectorType& n_proj) const ;

      Scalar projectionDistance(const VectorType& q) const ;
      Scalar projectionDistance(const Scalar* coords) const;

      int projectionDistanceAndNormal(const VectorType& q, Scalar& dist, VectorType& normal) const;
      int projectionDistanceAndNormal(const Scalar* coords, Scalar& dist, VectorType& normal) const;
      int projectionDistanceOnNormal(const VectorType& q, const VectorType& n, Scalar& dist) const;
      int projectionDistanceOnNormal(const float* coords, const float* normal, Scalar& dist) const;

      const std::vector<SimplexId>& point_indices() const {
        return m_point_indices;
      }

      void reserve(size_t n){
        m_point_indices.reserve(n);
      }

inline std::string params_str() const {
        std::stringstream msg;
        msg <<"uc "<<m_uc<<"\nul "<<m_ul.transpose();
        msg<<"\nuq "<<m_uq<<"\n";
        return msg.str();
      }
      void print(){
        std::stringstream msg;
        msg<<" Primitive: "<<type_str()<<"\n";
        msg<<" params: "<<params_str()<<"\n";
        SimplexId n = getNumberOfPoints();
        msg<<" nb points "<<n<<"\n";
        for(SimplexId i=0; i<n; i++){
          msg<<" "<<m_point_indices[i];
        }
        std::cout<<msg.str()<<std::endl;
      }

  };

  class PrimitiveSegmentation : public Debug{
    private:
      using Scalar = Primitive::Scalar;
      using VectorType = Primitive::VectorType;
      std::vector<Primitive> m_primitives{};
      std::vector<SimplexId>  m_primitive_growth{};

      std::vector<SimplexId>  m_labels{};
      std::vector<Scalar>  m_scales{};
      Scalar m_tolerance{0.01};
      


    public:
      // PrimitiveSegmentation(){};
      //

      inline void setTolerance(const Scalar& tol){
        m_tolerance=tol;
      }
      inline Scalar getTolerance() const{
        return m_tolerance;
      }
      inline Scalar getScale(const SimplexId iPrim) const{
        return m_scales[iPrim];
      }
      

      void reset(){
        m_primitives.clear();
        m_labels.clear();
      }
      void setNumberOfPoints(const SimplexId n){
        m_labels.resize(n, -1);
      }

      int addPrimitive(const Primitive& primitive, const Scalar scale){
        m_primitives.push_back(primitive);
        m_primitive_growth.push_back(0);
        m_scales.push_back(scale);
        return m_primitives.size()-1;
      }
      const Primitive& getPrimitive(int i) const {
        return m_primitives[i];
      }
      Primitive& getPrimitive(int i){
        return m_primitives[i];
      }

      SimplexId inline getNumberOfPrimitives() const{
        return m_primitives.size();
      }

      void addPointToPrimitive(SimplexId iPoint, SimplexId iPrim){
        m_labels[iPoint] = iPrim;
        getPrimitive(iPrim).addPoint(iPoint);
      }
      void addPointToPrimitive(SimplexId iPoint, SimplexId iPrim, const float* coordinates){
          Scalar maxDist = getTolerance()*getScale(iPrim);
          Scalar projDist = getPrimitive(iPrim).projectionDistance(coordinates+3*iPoint);
          if(projDist <= maxDist){
            addPointToPrimitive(iPoint, iPrim);
          }
      }
      void addPointToPrimitive(SimplexId iPoint, SimplexId iPrim, const float* coordinates, const float* normals){
        Scalar maxDist = getTolerance()*getScale(iPrim);
        Scalar projDist;
        int ret = getPrimitive(iPrim).projectionDistanceOnNormal(coordinates+3*iPoint, normals+3*iPoint, projDist);
        const float* normal = normals+3*iPoint;
        if(ret){ //projection yielded a solution
          if(projDist <= maxDist){
            addPointToPrimitive(iPoint, iPrim);
          }
        }
      }

      template <typename T>
      void growPrimitive(SimplexId iPrim, const float* coordinates, const float* normals, const T& tree){
        // std::cout<<"growing primitive "<<iPrim<<std::endl;
        Scalar maxDist = getTolerance()*getScale(iPrim);
        auto& prim = getPrimitive(iPrim);
        if(prim.kappa()>1e-2){
          return;
        }

        // prim.print_params();
        // std::cout<<"maxDist, "<<maxDist<<std::endl;
        const auto& indices = prim.point_indices();
        // for(const auto& iPoint : indices){
        //   std::cout<<" point "<<iPoint<<std::endl;
        // }

        for(SimplexId i=0; i<indices.size(); i++){
          SimplexId iPoint = indices[i];
          m_labels[iPoint] = iPrim;
        }
        std::vector<SimplexId> toAdd{};

        SimplexId last_growth = m_primitive_growth[iPrim];
        for(SimplexId i=last_growth; i<indices.size(); i++){
          // std::cout<<" point "<<iPoint<<std::endl;
          // std::cout<<" outta "<<tree.point_count()<<std::endl;
          SimplexId iPoint = indices[i];
          auto tree_query = tree.k_nearest_neighbors(iPoint, 10);
          for(auto iNeigh : tree_query){
            if(m_labels[iNeigh] != iPrim){ // neighbor not assigned
              Scalar projDist;
              const float* normal = normals+3*iNeigh;
              const float* coords =coordinates+3*iNeigh;
              int ret = prim.projectionDistanceOnNormal(coords, normal, projDist);
              if(ret){ //projection yielded a solution
               // Scalar np = projNormal.coeff(0)*normal[0] + projNormal.coeff(1)*normal[1] + projNormal.coeff(2)*normal[2];
               // np = sqrt(1-np*np);
                if(projDist <= maxDist){
                  toAdd.push_back(iNeigh);
                }
              }   
            }
          }
        }
        TTK_PSORT(this->threadNumber_, toAdd.begin(), toAdd.end());
        auto last = std::unique(toAdd.begin(), toAdd.end());
        toAdd.erase(last, toAdd.end());

        prim.reserve(toAdd.size());
        m_primitive_growth[iPrim] = prim.getNumberOfPoints();

        for(auto id : toAdd){
          addPointToPrimitive(id, iPrim);
        }

      }

      void print() const {
        std::stringstream msg;
        msg<<" Primitives Segmentation \n";
        msg<<" nb points "<<m_labels.size()<<"\n";
        SimplexId nPrim = getNumberOfPrimitives();
        msg<<" nb primitives "<<nPrim<<"\n";
        for(SimplexId i=0; i<nPrim; i++){
          const auto& prim = getPrimitive(i);
          msg<<" - "<<i<<": "<<prim.params_str()<<" ("<<prim.type_str()<<")"<<prim.getNumberOfPoints()<<" points\n";
        }
        std::cout<<msg.str()<<std::endl;

      }


  };

  class PCAnalysis : virtual public Debug {

  enum {Dim = 3};
  using Scalar = float;
  using VectorType = Eigen::Matrix<Scalar, Dim, 1>;
  using MatrixType = Eigen::Matrix<Scalar, Dim, Dim>;

  protected:
    int nbIterationsMLS_{5};
    double *geomVarOutputArray_{};
    double *k1OutputArray_{};
    double *k2OutputArray_{};
    double *curvOutputArray_{};
    float *etaOutputArray_{};
    ScaleSpaceSampling<MyPoint::Scalar, MyPoint> scaleSpaceSampler_{};
    bool aggregateFields_{false};

  public:
    PCAnalysis();


    inline void setAggregateFields(bool data){
      aggregateFields_ = data;
    }

    inline void setGeomVarOutputArray(double* data){
      geomVarOutputArray_ = data;
    }

    inline void setPrincipalCurvatureOutputArrays(double* k1, double* k2){
      k1OutputArray_ = k1;
      k2OutputArray_ = k2;
    }

    inline void setCurvOutputArray(double* data){
     curvOutputArray_ = data;
    }

    inline void setEtaOutputArray(float* data){
     etaOutputArray_ = data;
    }

    inline ScaleSpaceSampling<MyPoint::Scalar, MyPoint>* sampler(){
      return &scaleSpaceSampler_;
    }


    template<typename Fit, typename AccelerationStruct>
    int test_fit_kdtree(double tmax, Fit& _fit, Scalar* coordinates, Scalar* normals, const VectorType& _p, AccelerationStruct& _tree)
    {
      // Scalar tmax = 0.1;
      //
      
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

    template<typename Fit>
    void test_fit(Fit& _fit, std::vector<MyPoint>& _vecs, const VectorType& _p)
    {
      Scalar tmax = 0.1;

      // Set a weighting function instance
      _fit.setWeightFunc(WeightFunc(tmax));

      // Set the evaluation position
      _fit.init(_p);

      // Fit plane (method compute handles multipass fitting
      _fit.compute( _vecs.begin(), _vecs.end() );

      if(_fit.isStable())
      {
            std::cout << "Value of the scalar field at the initial point: "
                << _p.transpose()
                << " is equal to " << _fit.potential(_p)
                << std::endl;

            std::cout << "It's gradient at this place is equal to: "
                << _fit.primitiveGradient(_p).transpose()
                << std::endl;

            std::cout << "The initial point " << _p.transpose()              << std::endl
                << "Is projected at   " << _fit.project(_p).transpose() << std::endl;

            /* std::cout << "Value of the surface variation: " */
            /*     << _fit.surfaceVariation() */
            /*     << std::endl; */
      }
    }


    template<class Fit, class CoordsType>
    int executeSingleFit(double _tmax, CoordsType* _coordinates, CoordsType* _normals, int _nbPoints, VectorType& _p, SamplingGrid& _grid, double* _outputArray){

      Timer tm{};

      Fit fit;

      Ponca::KdTree<MyPoint> tree{};
      buildKDTree(_coordinates, _normals, _nbPoints, tree);

      computeMLS(_tmax, fit, _coordinates, _normals, _p, tree);

      singleFitToGridPotential(fit, _grid, _outputArray);

      printMsg("Total", 1, tm.getElapsedTime());

      return 1;
    }

    template<class Fit, class CoordsType>
    int executeFitAllPointsToGrid(double _tmax, CoordsType* _coordinates, CoordsType* _normals, int _nbPoints, SamplingGrid& _grid, double* _outputArray){

    Timer tm{};

    Ponca::KdTree<MyPoint> tree{};
    buildKDTree(_coordinates, _normals, _nbPoints, tree);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(size_t i=0; i < _grid.size(); i++){
          int I = i%_grid.dimX;
          int J = (i/_grid.dimX)%_grid.dimY;
          int K = (i/_grid.dimX)/_grid.dimY;
          double pX = _grid.originX + I*_grid.dX;
          double pY = _grid.originY + J*_grid.dY;
          double pZ = _grid.originZ + K*_grid.dZ;

          MyPoint::VectorType p;
          // p << pX, pY;//, pZ;
          p << pX, pY, pZ;

          Fit fit;
          // if(!(i%1000)){
          //   printMsg("Fitting point "+ std::to_string(i));
          // }
          computeMLS(_tmax, fit, _coordinates, _normals, p, tree);

          if(fit.isStable()){
            geomVarOutputArray_[i] = fit.geomVar();
            curvOutputArray_[i] = fit.kappa();
            _outputArray[i] = fit.potential(p);

            MyPoint::VectorType n;
            n=fit.eta_normalized();
            etaOutputArray_[3*i+0] = n.coeff(0);
            etaOutputArray_[3*i+1] = n.coeff(1);
            etaOutputArray_[3*i+2] = n.coeff(2);

          }else{
            geomVarOutputArray_[i] = std::numeric_limits<double>::quiet_NaN();
            _outputArray[i] = std::numeric_limits<double>::quiet_NaN();
            curvOutputArray_[i] = std::numeric_limits<double>::quiet_NaN();
            etaOutputArray_[3*i+0] = std::numeric_limits<double>::quiet_NaN();
            etaOutputArray_[3*i+1] = std::numeric_limits<double>::quiet_NaN();
            etaOutputArray_[3*i+2] = std::numeric_limits<double>::quiet_NaN();
          }
      }
      printMsg("Sampled scale = " + std::to_string(_tmax), 1, tm.getElapsedTime());
      return 1;
    }

    template<class Fit, class CoordsType>
    int executeFitAllPointsToGridAtScaleLevel(int _scaleLevel, CoordsType* _coordinates, CoordsType* _normals, int _nbPoints, SamplingGrid& _grid, double* _outputArray){

      std::cout<<"execyte fit all points at scale lvl"<<std::endl;
    printMsg("Samping scale " + std::to_string(_scaleLevel)+ "(t="+std::to_string(sampler()->scale(_scaleLevel))+")");
    printMsg("Sampling " +std::to_string(_grid.size()) + " points");

    Timer tm{};

    double tmax = sampler()->scale(_scaleLevel);
    Ponca::KdTree<MyPoint>& tree = sampler()->kdTree(_scaleLevel);
    unsigned char unstable = 0;


#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(size_t i=0; i < _grid.size(); i++){
          int I = i%_grid.dimX;
          int J = (i/_grid.dimX)%_grid.dimY;
          int K = (i/_grid.dimX)/_grid.dimY;
          double pX = _grid.originX + I*_grid.dX;
          double pY = _grid.originY + J*_grid.dY;
          double pZ = _grid.originZ + K*_grid.dZ;

          MyPoint::VectorType p;
          // p << pX, pY;//, pZ;
          p << pX, pY, pZ;

          Fit fit;
          // if(!(i%1000)){
          //   printMsg("Fitting point "+ std::to_string(i));
          // }
          bool isStable = computeMLS(tmax, fit, _coordinates, _normals, p, tree);

          if(isStable){
            // k1OutputArray_[i] = fit.k1();
            // k2OutputArray_[i] = fit.k2();
            geomVarOutputArray_[i] = fit.geomVar();
            curvOutputArray_[i] = fit.kappa();
            _outputArray[i] = fit.potential(p);
          }else{
            unstable = 255;
            if(!aggregateFields_){
            k1OutputArray_[i] = std::numeric_limits<double>::quiet_NaN();
            k2OutputArray_[i] = std::numeric_limits<double>::quiet_NaN();
            geomVarOutputArray_[i] = std::numeric_limits<double>::quiet_NaN();
            _outputArray[i] = std::numeric_limits<double>::quiet_NaN();
            curvOutputArray_[i] = std::numeric_limits<double>::quiet_NaN();
            }
          }
      }

    if(unstable){
        printWrn("At least one fit was not stable");
    }
      printMsg("Sampled scale = " + std::to_string(tmax), 1, tm.getElapsedTime());
      return 1;
    }

  template<typename Fit>
  int singleFitToGridPotential(Fit& _fit, SamplingGrid& _grid, double* outputArray){

    Timer tm{};

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(size_t i=0; i < _grid.size(); i++){
          int I = i%_grid.dimX;
          int J = (i/_grid.dimX)%_grid.dimY;
          int K = (i/_grid.dimX)/_grid.dimY;
          double pX = _grid.originX + I*_grid.dX;
          double pY = _grid.originY + J*_grid.dY;
          double pZ = _grid.originZ + K*_grid.dZ;

          MyPoint::VectorType p;
          // p << pX, pY; // pZ;
          p << pX, pY, pZ;
          
          // std::cout<<" point "<<pX<<" "<<pY<<" "<<pZ<<std::endl;
          // printMsg(" filling value "+ std::to_string(i));
          outputArray[i] = _fit.potential(p);
      }

      printMsg("Volume sampling", 1, tm.getElapsedTime(), threadNumber_);
      return 1;
  }




  template<typename CoordsType>
    void buildKDTree(CoordsType* _coordinates, CoordsType* _normals, int _nbPoints, Ponca::KdTree<MyPoint>& _tree){

      Timer tm{};
      struct VN{
          Scalar* coordinates;
          Scalar* normals;
          int n;

          int size() const {
            return n;
          }
      };

      _tree.build(VN({_coordinates, _normals, _nbPoints}), [&](VN buffer, Ponca::KdTree<MyPoint>::PointContainer &out){
                    out.reserve(buffer.size());
                    for(int i = 0; i<buffer.size(); i++){
                      out.push_back(MyPoint(buffer.coordinates, buffer.normals, i));
                    }});
      printMsg("kd-tree built", 1, tm.getElapsedTime());
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
      int go_on = stable && (i < 5);

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
        go_on = stable && (i<10);
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
      int go_on = stable && (i < 5);

      VectorType p = _p0;


      // mls loop
      while(go_on){

        p = _fit.project(p);

        // fit point p
        test_fit_kdtree(_tmax, _fit, _coordinates, _normals, p, _tree);
        stable = _fit.isStable();
        i++;
        go_on = stable && (i<10);
      } // end mls loop

      _p_proj=p;

      return stable;
    } 


    template<typename Fit, typename CoordsType>
    int computeMLS(double _tmax, Fit& _fit, CoordsType* _coordinates, CoordsType* _normals, const VectorType _p0, Ponca::KdTree<MyPoint>& _tree){

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
        go_on = stable && (i< nbIterationsMLS_);
      } // end mls loop

      return stable;
    } 
      

    // adapted from thibault's code in ScaleSampling.h

    void sampleScaleLog(std::vector<Scalar>& _scales, Scalar _scaleMin, Scalar _scaleMax, int _nbSamples, Scalar base){

      std::stringstream ss;
      ss << "Sampling "<<_nbSamples<<" scales btw "<<_scaleMin<<" and "<<_scaleMax<<" (base "<<base<<" )";
      printMsg(ss.str());

      _scales.clear();
      _scales.resize(_nbSamples);
      _scales[0] = _scaleMin;

      if(std::abs(base-1)<1e-3) {
        // uniform sampling
        Scalar delta0 = (_scaleMax-_scaleMin)/(_nbSamples-1);
        for(int i=1; i<_nbSamples; ++i) {
          _scales[i] = _scaleMin + delta0*i;
        }

      } else {
        // log sampling
        double lMin = log(_scaleMin)/log(10);
        double lMax = log(_scaleMax)/log(10);
        Scalar delta = (lMax-lMin) / (_nbSamples-1);
        for(int i=0; i<_nbSamples; ++i) {
          _scales[i] = std::pow(base, lMin + delta*i);
        }
      }
    }

    template<class Fit, class CoordsType>
    int executeFitPointList(double _tmax, CoordsType* _coordinates, CoordsType* _normals, int _nbPoints, CoordsType* _outputCoordinates, int _outputNbPoints, double* _outputArray){

    Timer tm{};

    Ponca::KdTree<MyPoint> tree{};
    buildKDTree(_coordinates, _normals, _nbPoints, tree);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(size_t i=0; i < _outputNbPoints; i++){

          MyPoint::VectorType p;
          // p << _outputCoordinates[3*i], _outputCoordinates[3*i+1];//, _outputCoordinates[3*i+2];
          p << _outputCoordinates[3*i], _outputCoordinates[3*i+1], _outputCoordinates[3*i+2];

          Fit fit;
          // if(!(i%1000)){
          //   printMsg("Fitting point "+ std::to_string(i));
          // }
          computeMLS(_tmax, fit, _coordinates, _normals, p, tree);

          if(fit.isStable()){
            geomVarOutputArray_[i] = fit.geomVar();
            curvOutputArray_[i] = fit.kappa_normalized();
            _outputArray[i] = fit.potential(p);
          }else{
            geomVarOutputArray_[i] = std::numeric_limits<double>::quiet_NaN();
            _outputArray[i] = std::numeric_limits<double>::quiet_NaN();
            curvOutputArray_[i] = std::numeric_limits<double>::quiet_NaN();
          }
      }
      // printMsg("Done scale = " + std::to_string(_tmax), 1, tm.getElapsedTime());
      return 1;
    }

    template<class Fit, class CoordsType>
    int projectPointListToScale(double _tmax, CoordsType* _coordinates, CoordsType* _normals, int _nbPoints, CoordsType* _outputCoordinates, std::vector<MyPoint::Scalar> _scales) {

    Timer tm{};

    Ponca::KdTree<MyPoint> tree{};
    buildKDTree(_coordinates, _normals, _nbPoints, tree);


#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(size_t i=0; i < _nbPoints; i++){

          MyPoint::VectorType p;
          // p << _outputCoordinates[3*i], _outputCoordinates[3*i+1];//, _outputCoordinates[3*i+2];
          p << _outputCoordinates[3*i], _outputCoordinates[3*i+1], _outputCoordinates[3*i+2];

          for(auto t : _scales){


              Fit fit;

              computeMLS(t, fit, _coordinates, _normals, p, tree);
              
              p = fit.project(p);
              _outputCoordinates[3*i] = p.coeff(0,0);
              _outputCoordinates[3*i+1] = p.coeff(1,0);
              _outputCoordinates[3*i+2] = p.coeff(2,0);
          }
      }
      return 1;
    }

template<class CoordsType>
  int fitNormals(double _t, CoordsType* _coordinates, CoordsType* _normals, int _nbPoints, CoordsType* _smoothed_normals){

      Timer tm{};


      Ponca::KdTree<MyPoint> tree{};
      buildKDTree(_coordinates, _normals, _nbPoints, tree);
      
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(int i=0; i<_nbPoints; i++) {
        myOrientedSphereFit fit;

        MyPoint::VectorType p;
        p << _coordinates[3*i], _coordinates[3*i+1], _coordinates[3*i+2];

        computeMLS(_t, fit, _coordinates, _normals, p, tree);
        if(!fit.isStable()){
          printWrn("fit not stable");
        }else{
          MyPoint::VectorType n = fit.primitiveGradient(p);
          double norm = n.coeff(0)*n.coeff(0)+n.coeff(1)*n.coeff(1)+n.coeff(2)*n.coeff(2);
          norm = sqrt(norm);
          n<< n.coeff(0)/norm, n.coeff(1)/norm, n.coeff(2)/norm;
          _smoothed_normals[3*i] = n.coeff(0);
          _smoothed_normals[3*i+1] = n.coeff(1);
          _smoothed_normals[3*i+2] = n.coeff(2);

          // std::cout<<"Point "<<i<<" , coords "<<p.transpose()<<", normals "<<n.transpose()<<std::endl;
        }
      }

      return 1;
  }

template<class Fit, class FitRIMLS, class SamplerType>
bool computeFitCoordinatesPlane(float* coordinates, float* normals, SimplexId
    nbPoints, const int FitDim, double* fitCoordinates, int* isStablePtr, float* mls3DCoordinates, float* mlsNormals, SamplerType&
    sampler, int iScale, std::vector<float>& resultBarycenter,
    std::vector<SimplexId>& sampling_indexes,
    float* mlsAllCoordinates,
    SimplexId BaryMaxId=-1,
    bool useNormalization=true, 
    bool useRIMLS = false){

    printMsg("compute fit coordinates: Plane");

    Timer tm{};

    using VectorType = pointCloud::MyPoint::VectorType;


    bool useKnnGraph = sampler.isUseKnnGraph();
    auto& tree = sampler.kdTree(iScale);
    auto& knnGraph = sampler.knnGraph();
    auto scale = sampler.scale(iScale);

    VectorType barycenter;
    if(BaryMaxId<0){
      BaryMaxId=nbPoints;
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
    }else{
      barycenter << coordinates[3*BaryMaxId+0], coordinates[3*BaryMaxId+1], coordinates[3*BaryMaxId+2]; 
      resultBarycenter.resize(3);
      resultBarycenter[0] = coordinates[3*BaryMaxId+0];
      resultBarycenter[1] = coordinates[3*BaryMaxId+1];
      resultBarycenter[2] = coordinates[3*BaryMaxId+2];
    }


    double tm_fitCoords=tm.getElapsedTime();
    unsigned char unstable = 0;
    // std::srand(std::time(0));
    // for each point, compute the fit in the barycenter basis
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId I = 0; I < sampling_indexes.size(); I++){
      
      SimplexId iPoint = sampling_indexes[I];

      VectorType p, p_proj;
      p << mlsAllCoordinates[3*iPoint], mlsAllCoordinates[3*iPoint+1], mlsAllCoordinates[3*iPoint+2];
      Fit pointFit;

      bool isStable{false}; 

      RIMLS<FitRIMLS,FitRIMLS> rimls{};
      if(useRIMLS){

        rimls.set_scale(scale);
        rimls.set_reweighting_step(8);
        rimls.compute(coordinates, normals, p, tree, p_proj);

        isStable = rimls.stable();


      }else{
        if(useKnnGraph){
          isStable = computeMLSKnnGraph(scale, pointFit, coordinates, normals, p, iPoint, knnGraph, p_proj);
        }else{
          isStable = computeMLS(scale, pointFit, coordinates, normals, p, tree, p_proj);
        }
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
        if(useRIMLS){
          auto rimlsFit = rimls.fit();
          mlsNormal = rimlsFit.primitiveGradient();
        }else{
          mlsNormal = pointFit.primitiveGradient();
        }

        for(int iDim=0; iDim<3; iDim++){
          mls3DCoordinates[3*I+iDim]=p_proj.coeff(iDim);
          mlsNormals[3*I+iDim]=mlsNormal.coeff(iDim);
          mlsAllCoordinates[3*iPoint+iDim] = p_proj.coeff(iDim);
        }

        
        bool normalize = useNormalization;
        MyPoint::VectorType ul;
        MyPoint::Scalar uc;

        switch(FitDim){
          default:

            if(useRIMLS){
              auto& rimlsFit = rimls.fit();
              uc=rimlsFit.potential(barycenter);
              ul=rimlsFit.primitiveGradient();
            }else{
              uc=pointFit.potential(barycenter);
              ul=pointFit.primitiveGradient();
            }
            // changeBasis(p, barycenter, uc, ul);

            ul.normalize();

            // uc
            fitCoordinates[FitDim*I+0]=uc;

            // ul
            fitCoordinates[FitDim*I+1]=ul.coeff(0);
            fitCoordinates[FitDim*I+2]=ul.coeff(1);
            fitCoordinates[FitDim*I+3]=ul.coeff(2);

            normalize=useNormalization;

            break;
        }


        // normalize;
        if(normalize){
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
        printWrn("At least one fit was not stable");
    }

    printMsg("Computed fit "+std::to_string(FitDim)+"D coordinates with scale "+std::to_string(scale), 1, tm.getElapsedTime()-tm_fitCoords);
    return true;
}

template<class Fit, class FitRIMLS, class SamplerType>
bool computeFitCoordinatesSphere(float* coordinates, float* normals, SimplexId
    nbPoints, const int FitDim, double* fitCoordinates, int* isStablePtr, float* mls3DCoordinates, float* mlsNormals, SamplerType&
    sampler, int iScale, std::vector<float>& resultBarycenter,
    std::vector<SimplexId>& sampling_indexes,
    float* mlsAllCoordinates,
    SimplexId BaryMaxId=-1,
    bool useNormalization=true, 
    bool useRIMLS = false){

    printMsg("compute fit coordinates");


    if constexpr (std::is_same_v<Fit, pointCloud::FitEllipsoid>){
      std::cout<<" ELLIPSOID FIT"<<std::endl;
    }

    Timer tm{};

    using VectorType = pointCloud::MyPoint::VectorType;


    bool useKnnGraph = sampler.isUseKnnGraph();
    auto& tree = sampler.kdTree(iScale);
    auto& knnGraph = sampler.knnGraph();
    auto scale = sampler.scale(iScale);

    VectorType barycenter;
    if(BaryMaxId<0){
      BaryMaxId=nbPoints;
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
    }else{
      barycenter << coordinates[3*BaryMaxId+0], coordinates[3*BaryMaxId+1], coordinates[3*BaryMaxId+2]; 
      resultBarycenter.resize(3);
      resultBarycenter[0] = coordinates[3*BaryMaxId+0];
      resultBarycenter[1] = coordinates[3*BaryMaxId+1];
      resultBarycenter[2] = coordinates[3*BaryMaxId+2];
    }

  // std::cout<<" fittin: "<<std::endl;
  // for(auto i : sampling_indexes){
  //   std::cout<<" "<<i;
  // }
  // std::cout<<std::endl;

    double tm_fitCoords=tm.getElapsedTime();
    unsigned char unstable = 0;
    // std::srand(std::time(0));
    // for each point, compute the fit in the barycenter basis
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId I = 0; I < sampling_indexes.size(); I++){
      
      SimplexId iPoint = sampling_indexes[I];

      VectorType p, p_proj;
      p << mlsAllCoordinates[3*iPoint], mlsAllCoordinates[3*iPoint+1], mlsAllCoordinates[3*iPoint+2];
      Fit pointFit;

      bool isStable{false}; 

      RIMLS<FitRIMLS, FitRIMLS> rimls{};
      if(useRIMLS){

        rimls.set_scale(scale);
        rimls.set_reweighting_step(8);
        rimls.compute(coordinates, normals, p, tree, p_proj);

        isStable = rimls.stable();


      }else{
        if(useKnnGraph){
          isStable = computeMLSKnnGraph(scale, pointFit, coordinates, normals, p, iPoint, knnGraph, p_proj);
        }else{
          isStable = computeMLS(scale, pointFit, coordinates, normals, p, tree, p_proj);
          // computeRIMLS(scale, pointFit, coordinates, normals, p, tree, p_proj);
        }
      }

      for(int iFitDim=0; iFitDim<FitDim; iFitDim++){
        fitCoordinates[FitDim*I+iFitDim]=0;
      }

      if(!isStable){
        unstable = 255;
        isStablePtr[I] = 0;
      }else{

        isStablePtr[I]=1;

        if(std::isnan(pointFit.kmin())){
          std::cout<<"FIT RESULT  "<<pointFit.m_uc<<" "<<pointFit.m_ul.transpose()<<" "<<pointFit.m_uq<<std::endl;
        }
        
        VectorType mlsNormal;
        if(useRIMLS){
          auto rimlsFit = rimls.fit();
          mlsNormal = rimlsFit.primitiveGradient();
        }else{
          mlsNormal = pointFit.primitiveGradient();
        }

        for(int iDim=0; iDim<3; iDim++){
          mls3DCoordinates[3*I+iDim]=p_proj.coeff(iDim);
          mlsNormals[3*I+iDim]=mlsNormal.coeff(iDim);
          mlsAllCoordinates[3*iPoint+iDim] = p_proj.coeff(iDim);
        }

        if(FitDim!=3){

          if(useRIMLS){
            rimls.fit().changeBasis(barycenter);
          } else{
          pointFit.changeBasis(barycenter);
          }
        }


        
        bool normalize = useNormalization;
        MyPoint::VectorType ul;
        MyPoint::Scalar uc;
        MyPoint::Scalar uq;
        MyPoint::Scalar kmin, kmax;

        switch(FitDim){
          case 1:
            uq = pointFit.m_uq;
            fitCoordinates[FitDim*I+0]=uq/2;
            normalize=false;
            break;
          case 2:

            if(useRIMLS){

              kmin = rimls.fit().kmin();
              kmax = rimls.fit().kmax();
            }else{

              kmin = pointFit.kmin();
              kmax = pointFit.kmax();
            }
            if(std::abs(kmin)>std::abs(kmax)){
              std::swap(kmin,kmax);
            }

            fitCoordinates[FitDim*I+0]=kmin;
            fitCoordinates[FitDim*I+1]=kmax;

            normalize=false;
            break;

          case 3:

            if(useRIMLS){

              kmin = rimls.fit().kmin();
              kmax = rimls.fit().kmax();
              uc = rimls.fit().m_uc;
            }else{

              kmin = pointFit.kmin();
              kmax = pointFit.kmax();
              uc = pointFit.m_uc;
            }
            if(std::abs(kmin)>std::abs(kmax)){
              std::swap(kmin,kmax);
            }

            fitCoordinates[FitDim*I+0]=uc;
            fitCoordinates[FitDim*I+1]=kmin;
            fitCoordinates[FitDim*I+2]=kmax;

            normalize=false;
            break;

          default:

            if(useRIMLS){
              auto& rimlsFit = rimls.fit();
              uc=rimlsFit.m_uc;
              ul=rimlsFit.m_ul;
              uq = rimlsFit.m_uq;
            }else{
              uc=pointFit.m_uc;
              ul=pointFit.m_ul;
              uq = pointFit.m_uq;
            }

            // uc
            fitCoordinates[FitDim*I+0]=uc;

            // ul
            fitCoordinates[FitDim*I+1]=ul.coeff(0);
            fitCoordinates[FitDim*I+2]=ul.coeff(1);
            fitCoordinates[FitDim*I+3]=ul.coeff(2);

            // uq
            fitCoordinates[FitDim*I+4]=uq;

            normalize=useNormalization;

            break;
        }

          
  
        // normalize;
        if(normalize){
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
        printWrn("At least one fit was not stable");
    }

    printMsg("Computed fit "+std::to_string(FitDim)+"D coordinates with scale "+std::to_string(scale), 1, tm.getElapsedTime()-tm_fitCoords);
    return true;
}

template<class Fit, class FitRIMLS, class SamplerType>
bool computeFitCoordinatesEllipsoid(float* coordinates, float* normals, SimplexId
    nbPoints, const int FitDim, double* fitCoordinates, int* isStablePtr, float* mls3DCoordinates, float* mlsNormals, SamplerType&
    sampler, int iScale, std::vector<float>& resultBarycenter,
    std::vector<SimplexId>& sampling_indexes,
    float* mlsAllCoordinates,
    SimplexId BaryMaxId=-1,
    bool useNormalization=true, 
    bool useRIMLS = false){

    printMsg("compute fit coordinates");


    if constexpr (std::is_same_v<Fit, pointCloud::FitEllipsoid>){
      std::cout<<" ELLIPSOID FIT"<<std::endl;
    }

    Timer tm{};

    using VectorType = pointCloud::MyPoint::VectorType;


    bool useKnnGraph = sampler.isUseKnnGraph();
    auto& tree = sampler.kdTree(iScale);
    auto& knnGraph = sampler.knnGraph();
    auto scale = sampler.scale(iScale);

    VectorType barycenter;
    if(BaryMaxId<0){
      BaryMaxId=nbPoints;
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
    }else{
      barycenter << coordinates[3*BaryMaxId+0], coordinates[3*BaryMaxId+1], coordinates[3*BaryMaxId+2]; 
      resultBarycenter.resize(3);
      resultBarycenter[0] = coordinates[3*BaryMaxId+0];
      resultBarycenter[1] = coordinates[3*BaryMaxId+1];
      resultBarycenter[2] = coordinates[3*BaryMaxId+2];
    }

  // std::cout<<" fittin: "<<std::endl;
  // for(auto i : sampling_indexes){
  //   std::cout<<" "<<i;
  // }
  // std::cout<<std::endl;

    double tm_fitCoords=tm.getElapsedTime();
    unsigned char unstable = 0;
    // std::srand(std::time(0));
    // for each point, compute the fit in the barycenter basis
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId I = 0; I < sampling_indexes.size(); I++){
      
      SimplexId iPoint = sampling_indexes[I];

      VectorType p, p_proj;
      p << mlsAllCoordinates[3*iPoint], mlsAllCoordinates[3*iPoint+1], mlsAllCoordinates[3*iPoint+2];
      Fit pointFit;

      bool isStable{false}; 

      RIMLS<FitRIMLS, FitRIMLS> rimls{};
      if(useRIMLS){

        rimls.set_scale(scale);
        rimls.set_reweighting_step(8);
        rimls.compute(coordinates, normals, p, tree, p_proj);

        isStable = rimls.stable();


      }else{
        if(useKnnGraph){
          isStable = computeMLSKnnGraph(scale, pointFit, coordinates, normals, p, iPoint, knnGraph, p_proj);
        }else{
          isStable = computeMLS(scale, pointFit, coordinates, normals, p, tree, p_proj);
          // computeRIMLS(scale, pointFit, coordinates, normals, p, tree, p_proj);
        }
      }

      for(int iFitDim=0; iFitDim<FitDim; iFitDim++){
        fitCoordinates[FitDim*I+iFitDim]=0;
      }

      if(!isStable){
        unstable = 255;
        isStablePtr[I] = 0;
      }else{

        isStablePtr[I]=1;

        if(std::isnan(pointFit.kmin())){
          std::cout<<"FIT RESULT  "<<pointFit.m_uc<<" "<<pointFit.m_ul.transpose()<<" "<<pointFit.m_uq<<std::endl;
        }
        
        VectorType mlsNormal;
        if(useRIMLS){
          auto rimlsFit = rimls.fit();
          mlsNormal = rimlsFit.primitiveGradient();
        }else{
          mlsNormal = pointFit.primitiveGradient();
        }

        for(int iDim=0; iDim<3; iDim++){
          mls3DCoordinates[3*I+iDim]=p_proj.coeff(iDim);
          mlsNormals[3*I+iDim]=mlsNormal.coeff(iDim);
          mlsAllCoordinates[3*iPoint+iDim] = p_proj.coeff(iDim);
        }

        // if(FitDim!=10){
        //   if(useRIMLS){
        //     rimls.fit().changeBasis(barycenter);
        //   } else{
        //   pointFit.changeBasis(barycenter);
        //   }
        // }


        
        bool normalize = useNormalization;
        MyPoint::VectorType ul;
        MyPoint::Scalar uc;
        MatrixType Uq;

        switch(FitDim){

          case 2:
            {

            if(useRIMLS){
              Uq = rimls.fit().m_uq;
            }else{
              Uq = pointFit.m_uq;
            }

            Scalar kmin = -2*Uq.col(0)[0];
            Scalar kmax = -2*Uq.col(1)[1];
            if(std::abs(kmin)>std::abs(kmax)){
              std::swap(kmin,kmax);
            }
            fitCoordinates[FitDim*I+0]=kmin;
            fitCoordinates[FitDim*I+1]=kmax;
            }
            break;

          case 10:


            if(useRIMLS){
              auto& rimlsFit = rimls.fit();
              uc=rimlsFit.m_uc;
              ul=rimlsFit.m_ul;
              Uq = rimlsFit.m_uq;
            }else{

              uc=pointFit.m_uc;
              ul=pointFit.m_ul;
              Uq = pointFit.m_uq;
            }

            changeBasis(p, barycenter, uc, ul, Uq);

            // uc
            fitCoordinates[FitDim*I+0]=uc;

            // ul
            fitCoordinates[FitDim*I+1]=ul.coeff(0);
            fitCoordinates[FitDim*I+2]=ul.coeff(1);
            fitCoordinates[FitDim*I+3]=ul.coeff(2);

            // Uq
            fitCoordinates[FitDim*I+4]=Uq.col(0)[0];
            fitCoordinates[FitDim*I+5]=Uq.col(1)[1];
            fitCoordinates[FitDim*I+6]=Uq.col(2)[2];
            fitCoordinates[FitDim*I+7]=Uq.col(1)[0];
            fitCoordinates[FitDim*I+8]=Uq.col(2)[0];
            fitCoordinates[FitDim*I+9]=Uq.col(2)[1];

            normalize = useNormalization;
            break;

          default:
            printErr("dim1=10");
            break;
        }

          
  
        // normalize;
        if(normalize){
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
        printWrn("At least one fit was not stable");
    }

    printMsg("Computed fit "+std::to_string(FitDim)+"D coordinates with scale "+std::to_string(scale), 1, tm.getElapsedTime()-tm_fitCoords);
    return true;
}



template<class Fit, class SamplerType>
bool computeDistanceMatrix(float* coordinates, float* normals, SimplexId
    nbPoints, float* fitCoordinates, int* isStablePtr, float* mls3DCoordinates, SamplerType&
    sampler, int iScale, std::vector<float>& resultBarycenter, SimplexId BaryMaxId, std::vector<std::vector<double>>& distMat, std::vector<SimplexId>& sampling_indexes, bool UsePlane=false, bool useNormalization=true, bool useProjectiveNorm=false){

    printMsg("compute fit coords and distance matrix");
    if(UsePlane){
    std::cout<<"using planes"<<std::endl;
    }else{
    std::cout<<"using spheres"<<std::endl;
    }

    Timer tm{};

    std::cout<<"here 0"<<std::endl;
    std::cout<<" size "<<sampling_indexes.size()<<std::endl;
    distMat.resize(sampling_indexes.size(), std::vector<double>(sampling_indexes.size(), 0));
    // for(SimplexId I = 0; I < sampling_indexes.size(); I++){
    //   distMat[I].resize()
    // }

    std::cout<<"here 2"<<std::endl;

    using VectorType = pointCloud::MyPoint::VectorType;


    bool useKnnGraph = sampler.isUseKnnGraph();
    auto& tree = sampler.kdTree(iScale);
    auto& knnGraph = sampler.knnGraph();
    auto scale = sampler.scale(iScale);

    VectorType barycenter;
    if(BaryMaxId<0){
      BaryMaxId=nbPoints;
      std::cout<<"Basis set to Full barycenter"<<std::endl;

      // find barycenter
      double tm_bary=tm.getElapsedTime();
      float bx=0, by=0, bz=0;
      for(SimplexId i=0; i<nbPoints; i++){
        bx+=coordinates[3*i];
        by+=coordinates[3*i+1];
        bz+=coordinates[3*i+2];
      }
      barycenter << bx/nbPoints, by/nbPoints, bz/nbPoints;
      printMsg("computed barycenter", 1, tm.getElapsedTime()-tm_bary);
      std::cout<<" bary "<<barycenter.transpose()<<std::endl;
      resultBarycenter.resize(3);
      resultBarycenter[0]=bx/nbPoints;
      resultBarycenter[1]=by/nbPoints;
      resultBarycenter[2]=bz/nbPoints;
    }else{
      barycenter << coordinates[3*BaryMaxId+0], coordinates[3*BaryMaxId+1], coordinates[3*BaryMaxId+2]; 
      resultBarycenter.resize(3);
      resultBarycenter[0] = coordinates[3*BaryMaxId+0];
      resultBarycenter[1] = coordinates[3*BaryMaxId+1];
      resultBarycenter[2] = coordinates[3*BaryMaxId+2];
    }


    double tm_fitCoords=tm.getElapsedTime();
    unsigned char unstable = 0;
    std::vector<Fit> fits(sampling_indexes.size());
    std::cout<<"here 3"<<std::endl;

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
        isStable = computeMLSKnnGraph(scale, fits[I], coordinates, normals, p, iPoint, knnGraph, p_proj);
      }else{
        isStable = computeMLS(scale, fits[I], coordinates, normals, p, tree, p_proj);
      }
      if(!isStable){
        unstable = 255;
        isStablePtr[iPoint] = 0;
        for(int iDim=0; iDim<5; iDim++){
          fitCoordinates[5*I+iDim]=0;
        }
      }else{

        MyPoint::VectorType ul;
        MyPoint::Scalar uc, uq;


        fits[I].changeBasis(barycenter);

        ul=fits[I].m_ul;
        uc=fits[I].m_uc;
        uq=fits[I].m_uq;//*std::abs(fits[I].m_uq);
        double norm = std::sqrt(ul.coeff(0)*ul.coeff(0)+ul.coeff(1)*ul.coeff(1)+ul.coeff(2)*ul.coeff(2) + uc*uc + uq*uq);

        for(int iDim=0; iDim<3; iDim++){
          mls3DCoordinates[3*iPoint+iDim]=p_proj.coeff(iDim);
        }



        fitCoordinates[5*I+2]=ul.coeff(0)/norm;
        fitCoordinates[5*I+3]=ul.coeff(1)/norm;
        fitCoordinates[5*I+4]=ul.coeff(2)/norm;

        // TAU (distance from barycenter to sphere)
        fitCoordinates[5*I]=uc/norm;

        // KAPPA (mean curvature)
        if(UsePlane){
          fitCoordinates[5*I+1]=0;
        }else{
          fitCoordinates[5*I+1]=uq/norm;
        }

        isStablePtr[iPoint]=1;

      }
    }

        // fill distMat
      for(SimplexId I = 0; I < sampling_indexes.size(); I++){

        SimplexId iPoint = sampling_indexes[I];

        VectorType p_i;
        p_i << coordinates[3*iPoint], coordinates[3*iPoint+1], coordinates[3*iPoint+2];


        MyPoint::VectorType ul;
        MyPoint::Scalar uc, uq;
        double norm;

        for(SimplexId J = I+1; J < sampling_indexes.size(); J++){
          
          SimplexId jPoint = sampling_indexes[J];

          VectorType p_j;
          p_j << coordinates[3*jPoint], coordinates[3*jPoint+1], coordinates[3*jPoint+2];

          // get center between p_j and p
          VectorType center;
          center << (coordinates[3*iPoint]+coordinates[3*jPoint])/2, (coordinates[3*iPoint+1]+coordinates[3*jPoint+1])/2, (coordinates[3*iPoint+2]+coordinates[3*jPoint+2])/2;

          fits[I].changeBasis(center);
          fits[J].changeBasis(center);

          ul=fits[I].m_ul;
          uc=fits[I].m_uc;
          uq=fits[I].m_uq;//*std::abs(fits[I].m_uq);
          norm = ul.coeff(0)*ul.coeff(0)+ul.coeff(1)*ul.coeff(1)+ul.coeff(2)*ul.coeff(2) + uc*uc + uq*uq;

          MyPoint::VectorType ul_j;
          MyPoint::Scalar uc_j, uq_j;
          ul_j=fits[J].m_ul;
          uc_j=fits[J].m_uc;
          uq_j=fits[J].m_uq;//*std::abs(fits[J].m_uq);

          auto norm_j = ul_j.coeff(0)*ul_j.coeff(0)+ul_j.coeff(1)*ul_j.coeff(1)+ul_j.coeff(2)*ul_j.coeff(2) + uc_j*uc_j + uq_j*uq_j;

          double dist=0;


          if(useProjectiveNorm){
            dist = ul.coeff(0)*ul_j.coeff(0)
                +  ul.coeff(1)*ul_j.coeff(1)
                +  ul.coeff(2)*ul_j.coeff(2)
                + uc*uc_j
                + uq*uq_j;
            if(std::isnan(dist)){
              std::cout<<"first nan"<<std::endl;
            }

            dist = dist*dist/(norm*norm_j);
            if(std::isnan(dist)){
              std::cout<<"second nan"<<std::endl;
            }
            if(dist>1){
              dist=1;
            }
            double t=dist;
            dist = std::acos(std::sqrt(dist));
            if(std::isnan(dist)){
              std::cout<<"third nan "<<t<<std::endl;
            }

          }else{
          norm = std::sqrt(norm);
          norm_j = std::sqrt(norm_j);
            if(!useNormalization){
              norm  = 1;
              norm_j =1;
            }

            dist = (ul.coeff(0)/norm-ul_j.coeff(0)/norm_j)*(ul.coeff(0)/norm-ul_j.coeff(0)/norm_j)
                      + (ul.coeff(1)/norm-ul_j.coeff(1)/norm_j)*(ul.coeff(1)/norm-ul_j.coeff(1)/norm_j)
                      + (ul.coeff(2)/norm-ul_j.coeff(2)/norm_j)*(ul.coeff(2)/norm-ul_j.coeff(2)/norm_j)
                      + (uc/norm-uc_j/norm_j)*(uc/norm-uc_j/norm_j)
                      + (uq/norm-uq_j/norm_j)*(uq/norm-uq_j/norm_j);
            dist = std::sqrt(dist);
          }
          distMat[I][J] = dist;
          distMat[J][I] = dist;
        }
      }


    if(unstable){
        printWrn("At least one fit was not stable, expect segfault");
    }
    // if(!one_stable){
    //     printWrn("Not fit was stable :()");

    // }

    printMsg("Computed fit 5D dist mat with scale "+std::to_string(scale), 1, tm.getElapsedTime()-tm_fitCoords);


    return true;
}

void changeBasis(const MyPoint::VectorType& oldBasisCenter,
                 MyPoint::VectorType& newBasisCenter,
                 MyPoint::Scalar& uc,
                 MyPoint::VectorType& ul,
                 MyPoint::MatrixType& Uq){
    auto diff = oldBasisCenter - newBasisCenter;
    uc = uc - ul.dot(diff) + diff.dot(Uq*diff);
    ul = ul - Uq*diff - Uq.transpose()*diff;
}

void changeBasis(const MyPoint::VectorType& oldBasisCenter,
                 const MyPoint::VectorType& newBasisCenter,
                 MyPoint::Scalar& uc,
                 MyPoint::VectorType& ul){
    auto diff = oldBasisCenter - newBasisCenter;
    uc = uc - ul.dot(diff);
}

void changeBasis(const MyPoint::VectorType& oldBasisCenter,
                 const MyPoint::VectorType& newBasisCenter,
                 MyPoint::Scalar& uc,
                 MyPoint::VectorType& ul,
                 MyPoint::Scalar& uq){
    auto diff = oldBasisCenter - newBasisCenter;
    uc = uc - ul.dot(diff) + diff.dot(uq*diff);
    ul = ul - Scalar(2.)*uq*diff;
}

  }; // PCAnalysis class
  
template<class DataPoint>
void getNeighborhood(SimplexId queryPointId, double scale, const float* coordinates,
    const Ponca::KdTree<DataPoint>& tree, std::vector<SimplexId>& neighborhood){

    typename DataPoint::VectorType p;
    p << coordinates[3*queryPointId], coordinates[3*queryPointId+1], coordinates[3*queryPointId+2];
    auto tree_query = tree.range_neighbors(p,scale);
    neighborhood.clear();
    for(auto id : tree_query){
      neighborhood.push_back(id);
    }
}
template<class DataPoint>
void getNeighborhoodKnnGraph(SimplexId queryPointId, double scale,
    const Ponca::KnnGraph<DataPoint>& tree, std::vector<SimplexId>& neighborhood){

    auto tree_query = tree.range_neighbors(queryPointId,scale);
    neighborhood.clear();
    for(auto id : tree_query){
      neighborhood.push_back(id);
    }
}

} // namespace pointCloud
} // namespace ttk
