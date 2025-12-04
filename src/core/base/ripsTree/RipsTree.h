/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::RipsTree
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %RipsTree class that computes for each vertex of a
/// triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'RipsTree'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include "DataTypes.h"
#include <Debug.h>

#include <PersistenceDiagramUtils.h>
#include <Triangulation.h>
#include <cstddef>
#include <sstream>

#include <Eigen/Eigen>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/Map.h>

using namespace ttk;

namespace ttk {

  struct RipsArc {
    SimplexId ccId{-1};
    double birth{-1};
    double death{-1};
    SimplexId birthId{-1};
    SimplexId deathId{-1};
    int type{-1};
    int treeId{-1};
    double diag{-1};
    double autosimilarity{-1};
  };

  struct RipsAugmentedArc {
    SimplexId ccId{-1};
    double birth{-1};
    double death{-1};
    SimplexId branchId{-1};
  };


  class RipsCC {
    private:
      SimplexId size_{1};
      SimplexId birthId_{};
      SimplexId deathId_{};
      SimplexId pairIdentifier_{};
      double death_{};
      double birth_{};
      SimplexId parent_{};
      std::vector<SimplexId> children_{};
      std::vector<double> diag_aabb_5D_{};
      SimplexId startNodeId_{};
      SimplexId endNodeId_{};
      SimplexId parent2D_{-1};
      SimplexId filteredId_{-1};
      // ComponentGridTuple parent2D_{-1,-1,-1};

    public:
      RipsCC(){}
      void copy(PersistencePair& p, SimplexId id = -1){
        birthId_=p.birth.id;
        deathId_=p.death.id;
        size_=p.cc_size;
        pairIdentifier_=id;
        death_=p.death.sfValue;
        birth_=p.birth.sfValue;
      }

      SimplexId getParent2D()const{
       return parent2D_;
      }

      void setParent2D(SimplexId data){
        parent2D_=data;
      }


      // void setParent2D(SimplexId c, SimplexId iratio, SimplexId ithreshold){
      //   parent2D_.ccId=c;
      //   parent2D_.iRatio=iratio;
      //   parent2D_.iThreshold=ithreshold;
      // }
      // ComponentGridTuple& getParent2D(){
      //   return parent2D_;
      // }

      // void setParent2D(const ComponentGridTuple& p){
      //   parent2D_.ccId=p.ccId;
      //   parent2D_.iRatio=p.iRatio;
      //   parent2D_.iThreshold=p.iThreshold;
      // }

      void setStartNodeId(SimplexId start){
        startNodeId_=start;
      }
      void setEndNodeId(SimplexId end){
        endNodeId_=end;
      }
      inline SimplexId size() const {
        return size_;
      }
      inline SimplexId deathId() const {
        return deathId_;
      }
      inline SimplexId birthId() const {
        return birthId_;
      }
      inline double birth() const {
        return birth_;
      }
      inline double death() const {
        return death_;
      }
      inline SimplexId parent() const{
        return parent_;
      }
      inline SimplexId startNodeId() const {
        return startNodeId_;
      }
      inline SimplexId endNodeId() const{
        return endNodeId_;
      }
      // parent helpers
      void setParent(SimplexId parent){
        parent_=parent;
      }

      // children helpers
      void clearChildren(){
        children_.clear();
      }

      void set_diag_aabb_5D(std::vector<double>& diag){
        diag_aabb_5D_ = std::move(diag);
      }
      const std::vector<double>& diag_aabb_5D() const{
        return diag_aabb_5D_;
      }

      const double getDiag() const{
        if(diag_aabb_5D_.size()>0){
          return diag_aabb_5D_[0];
        }else{
          return 0.0;
        }
      }

      void addChild(SimplexId child){
        children_.push_back(child);
      }

      SimplexId getChild(size_t i) const{
        if(i>=children_.size()){
          return -1;
        }
        return children_[i];
      }

      const std::vector<SimplexId>& getChildren() const{
        return children_;
      }

      SimplexId numberOfChildren() const{
        return children_.size();
      }

      std::string string() const{
        std::stringstream msg;
        msg<<" ids: "<<startNodeId_<<"-"<<endNodeId_<<", death: "<<death_<<", final size: "<<size_<<", parent: "<<parent_;
        return msg.str();
      }

      SimplexId filteredId(){
        return filteredId_;
      }
      void setFilteredId(SimplexId data){
        filteredId_ = data;
      }


      ~RipsCC(){
        children_.clear();
        diag_aabb_5D_.clear();
      }
  };

  class RipsTree : virtual public Debug {

  private:
    std::vector<SimplexId> tree_roots_{};
    std::vector<RipsCC> components_{};
    std::vector<SimplexId> nodeIds_{};
    std::vector<SimplexId> pointIdToNodeId_{};
    std::vector<RipsArc> ripsArcs_{};
    std::vector<RipsArc> ripsArcs2_{};
    DiagramType arcDiagram_{};
    double maxDeath_{};
    SimplexId filteredNumberOfComponents_{-1};
    std::vector<SimplexId> filteredIdToNodeId_{};
    SimplexId sizeMin_{-1};

  public:
    RipsTree();


    void clear(){
      components_.clear();
      nodeIds_.clear();
      tree_roots_.clear();
      pointIdToNodeId_.clear();

    }

    const std::vector<SimplexId>& getMergeOrder() const {
      return nodeIds_;
    }

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    }

    RipsCC& getComponent(SimplexId ccId) {
      return components_[ccId];
    }
    const RipsCC& getComponentSafe(SimplexId ccId) const {
      return components_[ccId];
    }

  void buildArcMT(std::vector<SimplexId>& nodes, std::vector<SimplexId>& arcs);
  void buildRipsTree(DiagramType& diagram);
  void buildRipsArcs(SimplexId sizeMin);
void buildRipsArcsWithBirth( const std::vector<double>& births);
  void buildRipsArcsWithBirth(SimplexId sizeMin,
                      double diagMax,
                      const std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                      const std::vector<std::vector<SimplexId>>& filteredToLocalId,
                      const std::vector<SimplexId>& pointsAtScaleStartIdx,
                      const std::vector<std::vector<double>>& fitCoordsVec,
                      const std::vector<double>& births,
                      const int dim);
  void buildRipsArcs(SimplexId sizeMin,
                      double diagMax,
                      const std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                      const std::vector<std::vector<SimplexId>>& filteredToLocalId,
                      const std::vector<SimplexId>& pointsAtScaleStartIdx,
                      const std::vector<std::vector<double>>& fitCoordsVec,
                      const int dim);
  void buildArcDiagram();
  SimplexId fillIndices(SimplexId ccId, SimplexId currentId);
  SimplexId sizeAtThreshold(SimplexId ccId, double threshold) const;
  SimplexId augmentedSizeAtThreshold(SimplexId ccId, double threshold, const SimplexId sizeMin) const;
  void buildPointIdToNodeId();
  SimplexId pointIdToNodeId(SimplexId pointId, double threshold) const;

  inline SimplexId fullSize(SimplexId ccId) const{
    return components_[ccId].endNodeId()-components_[ccId].startNodeId()+1;
  }


  void findContainingComponents(std::vector<SimplexId>& pointsIds, double threshold, std::vector<SimplexId>& ccIds) const;

  int getComponentPointIdsAtThreshold(SimplexId ccId, double threshold, std::vector<SimplexId>& points, const SimplexId sizeMin = 1) const;
  int getComponentNodeIdsAtThreshold(SimplexId ccId, double threshold, std::vector<SimplexId>& points, const SimplexId sizeMin=1) const;
  void getFullComponentNodeIds(SimplexId ccId, std::vector<SimplexId>& points) const;
  void getClusteringAtThreshold(double threshold, SimplexId sizeMin, std::vector<SimplexId>& sampling_indexes);
  void computeArcAutoSimilarity(const std::vector<std::vector<double>>& distMat, const std::vector<SimplexId>& sampling);
  // double computeJaccardIndex(const std::vector<SimplexId>& cc1, const std::vector<SimplexId>& cc2) const;
  // double computeJaccardIndex(SimplexId cc1, SimplexId cc2) const;

  template <typename FieldType>
  void computeMeanFieldAtThreshold(SimplexId ccId, double threshold, FieldType* field, int dim) const;
  template <typename FieldType>
  void computeCovarianceMatrixAtThreshold(SimplexId ccId, double threshold, FieldType* field, FieldType* meanField, int dim, SimplexId& primitiveResult) const;

  template <typename FieldType>
  void computeFullMeanField(SimplexId ccId, FieldType* field, int dim) const;

  template <typename FieldType>
  void dispatchScalarFullComponent(SimplexId ccId, FieldType* fieldArray, FieldType value) const;
  template <typename FieldType>
  void dispatchScalarComponentAtThreshold(SimplexId ccId, double threshold, FieldType* fieldArray, FieldType value, const SimplexId sizeMin = 1) const;

  SimplexId getNumberOfFilteredComponents(SimplexId sizeMin);

  // printing methods
  void printNodeIds() const;
  void printSizes(double threshold) const;

  inline SimplexId getPointId(SimplexId ccId) const{
    return components_[ccId].birthId();
  }
  inline SimplexId getNumberOfPoints()const{
    return components_.size();
  }
  inline SimplexId getFullSize(SimplexId ccId)const{
    return components_[ccId].size();
  }
  inline double getPersistence(SimplexId ccId)const{
    return components_[ccId].death();
  }
  inline SimplexId getSize(SimplexId ccId)const{
    return components_[ccId].size();
  }

  inline double maxDeath() const{
    return maxDeath_;
  }

  const std::vector<RipsArc>& getRipsArcs2() const{
    return ripsArcs2_;
  }
  const std::vector<RipsArc>& getRipsArcs() const{
    return ripsArcs_;
  }

  const DiagramType& getArcDiagram() const{
    return arcDiagram_;
  }

  SimplexId getParent2D(SimplexId ccId)const{
    return components_[ccId].getParent2D();
  }
  void setParent2D(SimplexId ccId, SimplexId p){
    components_[ccId].setParent2D(p);
  }

  SimplexId filteredIdToNodeId(SimplexId filteredId)const{
    std::cout<<"size of filtered index "<<filteredId<<"/"<<filteredIdToNodeId_.size()<<std::endl;
    return filteredIdToNodeId_[filteredId];
  }

void buildAABBDiagonals(const std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                  const std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                  const std::vector<SimplexId>& pointsAtScaleStartIdx,
                                  const std::vector<std::vector<double>>& fitCoordsVec,
                                  const int dim);

void compute5DAABBDiagonals(SimplexId iCC,
                             const std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                             const std::vector<std::vector<SimplexId>>& filteredToLocalId,
                             const std::vector<SimplexId>& pointsAtScaleStartIdx,
                             const std::vector<std::vector<double>>& fitCoordsVec,
                             const int dim,
                             std::vector<double>& output_aabb);
}; //RipsTree class

template <typename FieldType>
void RipsTree::computeCovarianceMatrixAtThreshold(SimplexId ccId, double threshold, FieldType* field, FieldType* meanField, int dim, SimplexId& primitiveResult) const{

    std::vector<SimplexId> componentPointIds{};
    getComponentPointIdsAtThreshold(ccId, threshold, componentPointIds);

    SimplexId repId = components_[ccId].birthId();
    SimplexId nbPoints = componentPointIds.size();

    std::vector<FieldType> covVector(dim*dim,0);


    std::vector<FieldType> mean(dim, 0);
    for(int i=0; i<dim; i++){
      mean[i] = meanField[dim*repId+i];
    }


    for(SimplexId& id : componentPointIds){
      for(int i=0; i<dim; i++){
        for(int j=i; j<dim; j++){
          covVector[i*dim+j] += (field[dim*id+i]-mean[i])*(field[dim*id+j]-mean[j])/nbPoints;
        }
      }
    }
    for(int i=0; i<dim; i++){
      for(int j=0; j<i; j++){
        covVector[i*dim+j] = covVector[j*dim+i];
    }}

    // std::cout<<"Cov vector matrix: "<<std::endl;
    // for(int i=0; i<dim; i++){
    //   for(int j=0; j<dim; j++){
    //     auto c = covVector[i*dim+j];
    //     std::cout<<" "<<c;
    //   }
    //   std::cout<<std::endl;
    // }

    Eigen::Matrix<FieldType, -1, -1, 0, 6, 6> covMatrix = Eigen::Map<Eigen::Matrix<FieldType, -1, -1>, Eigen::RowMajor>(covVector.data(), dim, dim);
    auto ev = covMatrix.eigenvalues();

    primitiveResult = 0;
    for(int i=0; i<dim; i++){
      if(ev[i].real()<0.1){
        primitiveResult+=1;
      }
    }
}

template <typename FieldType>
void RipsTree::computeMeanFieldAtThreshold(SimplexId ccId, double threshold, FieldType* field, int dim) const{

    std::vector<SimplexId> componentPointIds{};
    getComponentPointIdsAtThreshold(ccId, threshold, componentPointIds);

    std::vector<FieldType>  mean(dim, 0);

    for(SimplexId& id : componentPointIds){
      for(int i=0; i<dim; i++){
        mean[i] += field[dim*id+i]/componentPointIds.size();
      }
    }

    for(SimplexId& id : componentPointIds){
      for(int i=0; i<dim; i++){
        field[dim*id+i]=mean[i];
      }
    }
}

// template <typename FieldType>
// void RipsTree::computeFullMeanField(SimplexId ccId, FieldType* field, int dim) const{

//     std::vector<SimplexId> componentPoints{};
//     getFullComponentPoints(ccId, componentPoints);

//     std::vector<FieldType>  mean(dim, 0);

//     for(SimplexId& cc : componentPoints){
//       SimplexId id = getPointId(cc);
//       for(int i=0; i<dim; i++){
//         mean[i] += field[dim*id+i]/componentPoints.size();
//       }
//     }
//     for(SimplexId& cc : componentPoints){
//       SimplexId id = getPointId(cc);
//       for(int i=0; i<dim; i++){
//         field[dim*id+i]=mean[i];
//       }
//     }
// }

template <typename FieldType>
void RipsTree::dispatchScalarComponentAtThreshold(SimplexId ccId, double threshold, FieldType* fieldArray, FieldType value, const SimplexId  sizeMin) const{

    std::vector<SimplexId> componentPointIds{};
    getComponentPointIdsAtThreshold(ccId, threshold, componentPointIds, sizeMin);

    for(SimplexId& id : componentPointIds){
      fieldArray[id]=value;
    }
}

// template <typename FieldType>
// void RipsTree::dispatchScalarFullComponent(SimplexId ccId, FieldType* fieldArray, FieldType value) const{

//     std::vector<SimplexId> componentPoints{};
//     getFullComponentPoints(ccId, componentPoints);

//     for(SimplexId& cc : componentPoints){
//       SimplexId id = getPointId(cc);
//       // std::cout<<" cc "<<cc<<" is dispatched value "<<value<<" from "<<ccId<<std::endl;
//       fieldArray[id]=value;
//     }
// }

int findScaleId(SimplexId ripsBirthId, const std::vector<SimplexId>& pointsAtScaleStartIdx);

double AABBToDiag(std::vector<double>& aabb, int dim);
void compute5DAABB(const std::vector<SimplexId>& local_ids, const std::vector<double>& fitCoords, const int dim, std::vector<double>& outputAABB);
double compute5DAABBDiagonal(const std::vector<SimplexId>& local_ids, const std::vector<double>& fitCoords, const int dim);
void compute5DAABBMultiscale(const std::vector<std::vector<SimplexId>>& local_ids, const std::vector<std::vector<double>>& fitCoordsVec, const int dim, std::vector<double>& outputAABB);
double compute5DAABBDiagonalMultiscale(const std::vector<std::vector<SimplexId>>& local_ids, const std::vector<std::vector<double>>& fitCoordsVec, const int dim);
double compute5DAABBDiagonalMultiscale(const std::vector<SimplexId>& multiscale_ids,
                                            const std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                            const std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                            const std::vector<SimplexId>& pointsAtScaleStartIdx,
                                            const std::vector<std::vector<double>>& fitCoordsVec,
                                            const int dim);

}// namespace ttk
