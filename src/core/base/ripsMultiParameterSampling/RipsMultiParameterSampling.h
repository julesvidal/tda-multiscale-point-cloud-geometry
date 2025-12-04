/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::RipsMultiParameterSampling
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %RipsMultiParameterSampling class that computes for each vertex of a
/// triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'RipsMultiParameterSampling'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

#include <algorithm>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>
#include <numeric>

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <PCAnalysisUtils.h>
#include <Sampling.h>
#include "../../base/ripsPersistence/nanoflann.hpp"
#include "DataTypes.h"
#include "PersistenceDiagramUtils.h"
#include <RipsPersistence.h>
#include <RipsTree.h>
#include <PCAnalysis.h>

#include <Eigen/Eigen>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/Map.h>

using namespace ttk;
using PrimitiveType = pointCloud::Primitive::PrimitiveType;

namespace ttk {


  class Timings: public Debug{
    public:
      enum class STEP{INIT=0,
                      SUBSAMPLING,
                      POISSON,
                      FITTING,
                      EDGES,
                      MST,
                      EPSILON_SAMPLING,
                      DENSITY,
                      PERSISTENCE,
                      RIPS_TREE,
                      RIPS_ARCS, 
                      PROCESS_ARCS,
                      GROW_COVERAGE,
                      TOTAL,
                      count};

    private:
      std::vector<double> m_timings{};
      Timer m_tm{};

     const int m_nbSteps = static_cast<int>(STEP::count);

    public:

      void init(){
        m_timings.resize(m_nbSteps);
        std::fill(m_timings.begin(), m_timings.end(), 0.0);
        m_tm.reStart();
      }

      void clock(const STEP step){
        double t0 = m_timings[static_cast<int>(step)];
        double t1 = m_tm.getElapsedTime();
        m_timings[static_cast<int>(step)] = t1 - t0;
      }

      void save(const std::string& filename){
          std::ofstream myfile;
          myfile.open(filename, std::ofstream::out | std::ofstream::app);

          std::stringstream line;
          for(int i = 0; i<m_nbSteps; i++){
            line <<m_timings[i];
            if(i<m_nbSteps-1){
              line<<", ";
            }else{
              line<<"\n";
            }
          }
          myfile << line.str();
      }
  };


  struct ComponentGridTuple {
      SimplexId ccId;
      SimplexId iRatio;
      SimplexId iThreshold;

      ComponentGridTuple(SimplexId c, SimplexId iratio, SimplexId ithreshold) : ccId(c), iRatio(iratio), iThreshold(ithreshold){}

      bool operator==(const ComponentGridTuple &rhs) const {
        return this->iRatio == rhs.iRatio 
              && this->iThreshold == rhs.iThreshold 
              && this->ccId == rhs.ccId;
      }

      bool operator<(const ComponentGridTuple &rhs) const {
        return this->iRatio < rhs.iRatio 
              or (this->iRatio == rhs.iRatio  and  
                  (this->iThreshold > rhs.iThreshold or(this->iThreshold == rhs.iThreshold && this->ccId < rhs.ccId)));
      }
    };

  /**
   * The RipsMultiParameterSampling class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class RipsMultiParameterSampling : virtual public Debug {

    private:

    protected:
      std::vector<RipsTree> trees_{};
      RipsTree treeMS_{};
      DiagramType diagramMS_{};
      DiagramType ripsDiagramMS_{};
      // std::vector<double> ratios_{};
      // std::vector<double> thresholds_{};
      std::vector<std::vector<double>> persistenceGrid_;
      // std::map<SimplexId,std::pair<SimplexId, double>> parents2D_{};
      std::vector<std::vector<SimplexId>> children2D_{};
      std::vector<SimplexId> parents2D_{};
      std::vector<double> persistence2d_{};
      std::vector<SimplexId>parentList_{};
      std::vector<SimplexId>parentOrder_{};
      std::vector<std::vector<SimplexId>> fathers_{};
      std::vector<std::vector<SimplexId>> sons_{};
      std::vector<SimplexId> globalToLocal_{};
      std::vector<SimplexId> globalToFullCompact_{};
      std::vector<SimplexId> fullCompactToGlobal_{};
      SimplexId nbComponents_{};
      double JaccardConstraint{0.5};
      int NbCCMax{1};
      int SizeMin{1};
      double DiagMax{-1};
      double DiagMaxDisplay{-1};

      // ScaleSpaceSampling<pointCloud::MyPoint::Scalar> sampler_{};

      bool UseNF{false};
      bool HierarchicalSampling{true};
      double ScaleSamplingMin{0.01};
      double ScaleSamplingMax{1};
      double ScaleSamplingBase{10};
      int ScaleSamplingCount{5};
      double SamplingFactor{0.1};

    public:
      RipsMultiParameterSampling();
      double computeJaccardIndex(const std::vector<SimplexId>& cc1, const std::vector<SimplexId>& cc2) const;
      void findMatchingComponent(SimplexId ccId,
          SimplexId inputTreeId, double inputThreshold,
          SimplexId targetTreeId, SimplexId targetThresholdId, double targetThreshold,
          SimplexId nCol,
          SimplexId& matchingCC,
          SimplexId& localCCId,
          double& jaccard);

      SimplexId getFullGridId(SimplexId ccId, SimplexId iRow, SimplexId iCol, SimplexId nCol, SimplexId n)const;
      void getGridTuple(SimplexId fullGridId, SimplexId nCol, SimplexId n, SimplexId& ccId, SimplexId& iRatio, SimplexId& iThreshold)const;
      SimplexId localToGlobal(SimplexId ccId, SimplexId iRow, SimplexId iCol, SimplexId nCol);
      void compute2DParents(std::vector<std::vector<double>> persistenceGrid);
      void compute2DParentsOneToOne(std::vector<std::vector<double>> persistenceGrid, const std::vector<double>& ratioVec);
      void buildGlobalToLocal(std::vector<std::vector<double>> persistenceGrid, const std::vector<double>& ratioVec);
      void findPersistentComponents(const std::vector<std::vector<double>>& persistenceGrid, const std::vector<double>& ratioVec, std::vector<double>& persistence2d);
      void propagate2dParent(SimplexId parentId, SimplexId ccId, SimplexId originRow, SimplexId originCol, SimplexId targetRow, SimplexId targetCol, const std::vector<std::vector<double>>& persistenceGrid, const std::vector<double>& ratioVec);
      void outputGraph(std::string outputFileName, const std::vector<std::vector<double>>& persistenceGrid);
void writeCoords(std::string outputFileName,
                                    const SimplexId nb,
                                    std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                    std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                    std::vector<std::vector<SimplexId>>& sampling_indexes,
                                    std::vector<SimplexId>& pointsAtScaleStartIdx,
                                    std::vector<std::vector<double>>& fitCoordinates,
                                    int dim,
                                    const float* coordinates) const;
      void writeMatrixAndCoords(std::string outputFileName,
          const SimplexId nbPointsMS,
                                    std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                    std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                    std::vector<std::vector<SimplexId>>& sampling_indexes,
                                    std::vector<SimplexId>& pointsAtScaleStartIdx,
                                    std::vector<std::vector<double>>& fitCoordinates,
                                    std::vector<std::vector<float>>& mls3DCoordinates,
                                    std::vector<float>& pointCloudBarycenter,
                                    std::function<double(const double*, const double*, const float*, const float*, const float*, int)> distanceFunction,
                                    int dim,
                                    const float* coordinates) const;
      // void writeMatrixAndCoords(std::string outputFileName,
      //     const std::vector<std::vector<double>>& mat,
      //     const std::vector<SimplexId>& localToGlobalIdx,
      //     const std::vector<SimplexId>& filteredToLocalIdx,
      //     const float* coordinates) const ;
      bool isIncluded(SimplexId ccId0, SimplexId iRatio0, SimplexId iThreshold0, SimplexId ccId1, SimplexId iRatio1, SimplexId iThreshold1) const;
      bool isPointSetIncluded(const std::vector<SimplexId>& cc1, const std::vector<SimplexId>& cc2) const;
      void coordsToDistanceMatrix(float* coords,
          int dim,
          std::vector<SimplexId> sampling,
          std::vector<std::vector<double>>& distMat);

      template <typename FieldType>
        void computeSingleComponentMean(FieldType* field, std::vector<FieldType>& mean, int dim, SimplexId n)const ;

      template <typename FieldType>
        void computeSinglePrimitiveType(FieldType* normalsPtr, std::vector<FieldType>& meanNormalsVec, SimplexId nbPoints, SimplexId& primitiveType, std::vector<FieldType>& covVec)const;

      template <typename T>
        double pointDistance(T* p0, T* p1, int dim);

      template <typename T>
        double pointSquaredDistance(T* p0, T* p1, int dim);


      template <typename EpsSamplerType>
        int retrieveSizeMSVector(EpsSamplerType& epsilonSamplers,
            std::vector<SimplexId>& sizes_ms,
            std::vector<SimplexId>& pointsAtScaleStartIdx,
            std::vector<std::vector<SimplexId>>& filteredIdxToEpsSampling,
            std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId );

      template <typename EpsSamplerType>
        int performEpsilonSampling(EpsSamplerType& epsilonSamplers,
            double epsilon,
            const SimplexId epsilonNbMaxPoint,
            SimplexId& nbPointsMS,
            std::vector<std::vector<double>>& fitCoordsVec,
            std::vector<SimplexId>& pointsAtScaleStartIdx,
            std::vector<std::vector<SimplexId>>& filteredToLocalIdVec,
            std::vector<std::vector<SimplexId>>& filteredIdxToEpsSampling,
            std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId );


template <typename EpsSampler>
int getPrimitives(pointCloud::PrimitiveSegmentation& seg,
                                              const RipsTree& tree,
                                              SimplexId ccId,
                                              const PrimitiveType primitiveType,
                                              double threshold,
                                              float maxDist,
                                              const float* coordinates,
                                              const float* normals,
                                              std::vector<std::vector<double>>& fitCoordinates,
                                              const std::vector<float>& scales,
                                              const std::vector<float>& pcBarycenter,
                                              std::vector<std::vector<SimplexId>>& sampling_indexes,
                                              std::vector<std::vector<SimplexId>>& filteredToLocalIdVec,
                                              std::vector<SimplexId>& pointsAtScaleStartIdx,
                                              std::vector<EpsSampler> epsilonSampler,
                                              std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId){

  int dim{};
  switch(primitiveType){
    case PrimitiveType::PLANE: //plane
      dim = 4;
      break;
    case PrimitiveType::SPHERE: // sphere
      dim=5;
      break;
    case PrimitiveType::QUADRIC: //quadric
      dim=10;
      break;
    default:
      printErr("Wrong primitive type");
      return -1;
  }

  SimplexId nbScales = sampling_indexes.size();
  if(epsilonIdToFilteredId.size()!=nbScales or filteredToLocalIdVec.size()!=nbScales or epsilonSampler.size()!=nbScales or pointsAtScaleStartIdx.size()!=nbScales){
    std::cout<<"SIZE DOES NOT MATCH SCALE NUMBER"<<std::endl;
  }

  SimplexId centroid = tree.getComponentSafe(ccId).birthId();
  int centroid_scale = findScaleId(centroid, pointsAtScaleStartIdx);
  SimplexId centroid_epsId = centroid-pointsAtScaleStartIdx[centroid_scale]; // epsilon id
  SimplexId centroid_filteredId = epsilonIdToFilteredId[centroid_scale][centroid_epsId];
  SimplexId centroid_localId = filteredToLocalIdVec[centroid_scale][centroid_filteredId];

  float scale = scales[centroid_scale];

  double* fitParams = &fitCoordinates[centroid_scale][dim*centroid_localId];
  Primitive::VectorType basisCenter;
  basisCenter << Primitive::Scalar(pcBarycenter[0]),Primitive::Scalar(pcBarycenter[1]),Primitive::Scalar(pcBarycenter[2]);
  Primitive prim(fitParams,
                 primitiveType,
                basisCenter);
  int primId = seg.addPrimitive(prim, scale);


  std::vector<SimplexId> pointIds;
  tree.getComponentPointIdsAtThreshold(ccId, threshold, pointIds, 0);

  std::vector<SimplexId> nodeIds;
  tree.getComponentNodeIdsAtThreshold(ccId, threshold, nodeIds, 0);


  SimplexId nbPoints = pointIds.size();
  if(nodeIds.size()!=nbPoints){
    std::cout<<"MISMATCH SIZE "<<nodeIds.size()<<" vs "<<nbPoints<<std::endl;
  }

  std::vector<SimplexId> toAdd{};
  // std::vector<std::vector<SimplexId>> localIds(nbScales);

  for(int i=0; i<nbPoints; i++){
    int iScale = findScaleId(pointIds[i], pointsAtScaleStartIdx);
    SimplexId epsId = pointIds[i]-pointsAtScaleStartIdx[iScale]; // epsilon id
    SimplexId filteredId = epsilonIdToFilteredId[iScale][epsId];

    SimplexId localId = filteredToLocalIdVec[iScale][filteredId];
    SimplexId globalId = sampling_indexes[iScale][localId];

    toAdd.push_back(globalId);


    // find children in epsilon sampling
    auto& children = epsilonSampler[iScale].getChildren(epsId);

    for(SimplexId iChild=1; iChild<children.size(); iChild++){
      auto child_filteredId = children[iChild];
      SimplexId child_localId = filteredToLocalIdVec[iScale][child_filteredId]; // localId
      SimplexId child_globalId = sampling_indexes[iScale][child_localId];
      toAdd.push_back(child_globalId);
    }
  }

  TTK_PSORT(this->threadNumber_, toAdd.begin(), toAdd.end());
  auto last = std::unique(toAdd.begin(), toAdd.end());
  toAdd.erase(last, toAdd.end());

  seg.getPrimitive(primId).reserve(toAdd.size());
  for(const SimplexId& iPoint : toAdd){
    seg.addPointToPrimitive(iPoint, primId, coordinates, normals);
  }
  return 0;
}

  };// RipsMultiParameterSampling class
                                                                                                               //

template <typename FieldType>
void RipsMultiParameterSampling::computeSingleComponentMean(FieldType* field, std::vector<FieldType>& mean, int dim, SimplexId n)const {

    mean.resize(dim);
    for(SimplexId iPoint=0; iPoint<n; iPoint++){
      for(int i=0; i<dim; i++){
        mean[i] += field[dim*iPoint+i]/n;
      }
    }
}


template <typename FieldType>
void RipsMultiParameterSampling::computeSinglePrimitiveType(FieldType* normalsPtr, std::vector<FieldType>& meanNormalsVec, SimplexId nbPoints, SimplexId& primitiveType, std::vector<FieldType>& covVec) const{

    int dim=3;
    std::vector<FieldType> covVector(dim*dim,0);

    for(SimplexId iPoint=0; iPoint<nbPoints; iPoint++){
      for(int i=0; i<dim; i++){
        for(int j=i; j<dim; j++){
          covVector[i*dim+j] += (normalsPtr[dim*iPoint+i]-meanNormalsVec[i])*(normalsPtr[dim*iPoint+j]-meanNormalsVec[j])/nbPoints;
        }
      }
    }

    for(int i=0; i<dim; i++){
      for(int j=0; j<i; j++){
        covVector[i*dim+j] = covVector[j*dim+i];
    }}

    Eigen::Matrix<FieldType, -1, -1, 0, 6, 6> covMatrix = Eigen::Map<Eigen::Matrix<FieldType, -1, -1>, Eigen::RowMajor>(covVector.data(), dim, dim);
    auto ev = covMatrix.eigenvalues();

    // find max of eigenvalues
    FieldType maxEv = std::numeric_limits<FieldType>::min();

    for(int i=0; i<dim; i++){
      if(ev[i].real()>maxEv){
        maxEv = ev[i].real();
      }
    }


    primitiveType = 0;
    covVec.resize(dim);

    for(int i=0; i<dim; i++){
      covVec[i] = ev[i].real();
      if(ev[i].real()<0.1){
        primitiveType+=1;
      }
    }
}

template <typename T>
double ttk::RipsMultiParameterSampling::pointSquaredDistance(T* p0, T* p1, int dim){
  double d=0;
  for(int i =0; i<dim; i++){
    // std::cout<<p0[i]<<" "<<p1[i]<<std::endl;
    d+=(p1[i]-p0[i])*(p1[i]-p0[i]);
  }
  return d;
}

template <typename T>
double ttk::RipsMultiParameterSampling::pointDistance(T* p0, T* p1, int dim){
  return sqrt(pointSquaredDistance(p0, p1, dim));
}
template <typename EpsSamplerType>
int ttk::RipsMultiParameterSampling::performEpsilonSampling(EpsSamplerType& epsilonSamplers,
                                                              double epsilon,
                                                              const SimplexId epsilonNbMaxPoint,
                                                               SimplexId& nbPointsMS,
                                                               std::vector<std::vector<double>>& fitCoordsVec,
                                                               std::vector<SimplexId>& pointsAtScaleStartIdx,
                                                               std::vector<std::vector<SimplexId>>& filteredToLocalIdVec,
                                                               std::vector<std::vector<SimplexId>>& filteredIdxToEpsSampling,
                                                               std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId ){
    nbPointsMS = 0;
    int nbScales = epsilonSamplers.size();

    double current_epsilon=epsilon;
    double max_epsilon = 0.0;
    SimplexId min_nb_reps = std::numeric_limits<SimplexId>::max();

    Timer tm_epsilon_sampling{};

    for(int iScale=0; iScale<nbScales; iScale++){
      pointsAtScaleStartIdx[iScale] = nbPointsMS;

      SimplexId nbFilteredPointsAtScale = filteredToLocalIdVec[iScale].size();

      printMsg("Epsilon Sampling for scale "+std::to_string(iScale));
      auto& epsSampler = epsilonSamplers[iScale];
      epsSampler.setEpsilon(current_epsilon);
      epsSampler.setNbRepMax(epsilonNbMaxPoint);
      int rng_seed=iScale;
      epsSampler.buildSampling(nbFilteredPointsAtScale, fitCoordsVec[iScale].data(), filteredToLocalIdVec[iScale], epsilonIdToFilteredId[iScale], rng_seed);
      
      std::cout<<"epsilon "<<epsSampler.epsilon()<<std::endl;
      current_epsilon=std::max(epsSampler.epsilon(), epsilon);

      // build filtered to epsilon
      filteredIdxToEpsSampling[iScale].resize(filteredToLocalIdVec[iScale].size());
      for(SimplexId id_eps=0; id_eps <epsilonIdToFilteredId[iScale].size(); id_eps++){
        filteredIdxToEpsSampling[iScale][epsilonIdToFilteredId[iScale][id_eps]] = id_eps;
      }

      SimplexId nbEpsilonPointsAtScale = epsilonIdToFilteredId[iScale].size();
      nbPointsMS += nbEpsilonPointsAtScale;

      if(epsSampler.epsilon()>max_epsilon){
        max_epsilon = epsSampler.epsilon();
      }
      if(nbEpsilonPointsAtScale<min_nb_reps){
        min_nb_reps = nbEpsilonPointsAtScale;
      }
    }

    std::cout<<" Epsilon Sampling done in "<<tm_epsilon_sampling.getElapsedTime()<<"s."<<std::endl;
    std::cout<<" max_epsilon "<<max_epsilon<<",   min_nb_reps "<<min_nb_reps<<std::endl;
    return 0;
}

template <typename EpsSamplerType>
int ttk::RipsMultiParameterSampling::retrieveSizeMSVector(EpsSamplerType& epsilonSamplers,
                                                               std::vector<SimplexId>& sizes_ms,
                                                               std::vector<SimplexId>& pointsAtScaleStartIdx,
                                                               std::vector<std::vector<SimplexId>>& filteredIdxToEpsSampling,
                                                               std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId ){

      int nbScales = epsilonSamplers.size();

      // this should be able to be done in parallel
      for(int iScale=0; iScale<nbScales; iScale++){
        auto& epsReps = epsilonSamplers[iScale].representant();
        for(SimplexId i=0; i<epsReps.size(); i++){
          SimplexId rep_filtered = epsReps[i]; // monoscale full filtered index
          SimplexId rep_ms = filteredIdxToEpsSampling[iScale][rep_filtered] + pointsAtScaleStartIdx[iScale];
          sizes_ms[rep_ms] += 1;
        }
      }
      return 0;
}

}// namespace ttk
