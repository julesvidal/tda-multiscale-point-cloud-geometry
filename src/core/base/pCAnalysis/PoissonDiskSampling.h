#pragma once

//ttk includes
#include <Debug.h>
#include <Shuffle.h>

// Ponca includes
#include "DataTypes.h"
#include <Ponca/SpatialPartitioning>

//std includes
#include <cmath>
#include <PCAnalysisUtils.h>
#include <numeric>
#include <random>

using namespace ttk;
using namespace pointCloud;


template <class Scalar, typename Point>
class ScaleSpaceSampling : public ttk::Debug{


  public:
  ScaleSpaceSampling(){
    this->setDebugMsgPrefix("ScaleSpaceSampling");
  };


  // access
  public:
  void setSamplingFactor(double factor){
    m_factor = factor;
  }
  double samplingFactor(){
    return m_factor;
  }
  inline Ponca::KdTree<Point>& kdTree(int lvl){
    return m_trees[lvl >=0 ? lvl : 0];
  }
  inline Ponca::KdTree<Point>& kdTree(){
    return kdTree(m_curLevel);
  }
  inline Ponca::KnnGraph<Point>& knnGraph(){
    return *m_KnnGraph;
  }
  inline Scalar scale(){
    return scale(m_curLevel);
  }
  inline Scalar scale(int lvl){
    return m_scales[lvl];
  }
  inline const std::vector<Scalar>& scales(){
    return m_scales;
  }
  inline int level(int i){
    return m_levels[i];
  }
  inline int scale_count() const{
    return m_scales.size();
  }
  inline bool is_sub_sampled(int iPoint, int level){
    return (level <= m_levels[iPoint]);
  };

  inline bool isUseKnnGraph() const{
    return useKnnGraph_;
  }
  inline void setUseKnnGraph(bool data){
    useKnnGraph_=data;
  }


  void print(){
    this->printMsg("ScaleSpaceSampling:");
    std::cout<<" Scales : ";
    for(auto s : m_scales){
      std::cout<<" "<<s;
    }
    std::cout<<"\n Poisson Disk levels:\n";
    for(int i=0; i< m_scales.size(); i++){
      std::cout<<" level "<<i<<": ";
      for(auto l : m_levels){
        if(l>=i) std::cout<<" "<<l;
        else std::cout<<" _";
      }
      std::cout<<std::endl;
    }
      std::cout<<std::endl;
  }

  public:

  public:
  void sampleScaleLog(Scalar _scaleMin, Scalar _scaleMax, int _nbSamples, double base = 10);

  void poissonDiskSampling(int _nbPoints, Scalar* _coordinates, Scalar* _normals);

  void buildKDTree(int level, Scalar* _coordinates, Scalar* _normals, int _nbPoints);
  void buildKnnGraph();

  void filterUnstablePoints(int level, Scalar* _coordinates, int _nbPoints, std::vector<SimplexId>& densitySampling, int* isStable);



  private:

  int m_curLevel{0};
  std::vector<int> m_levels{};
  std::vector<Scalar> m_scales{};
  double m_factor{0.1};
  std::vector<Ponca::KdTree<Point>> m_trees{};

  bool useKnnGraph_{false};
  Ponca::KnnGraph<Point>* m_KnnGraph{};
};

template <class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::sampleScaleLog(Scalar _scaleMin, Scalar _scaleMax, int _nbSamples, double base){

      std::stringstream ss;
      ss << "Sampling "<<_nbSamples<<" scales btw "<<_scaleMin<<" and "<<_scaleMax<<" (base "<<base<<" )";
      printMsg(ss.str());

      m_scales.clear();
      m_scales.resize(_nbSamples);
      m_scales[0] = _scaleMin;
      if(_nbSamples==1){
        return;
      }

      if(std::abs(base)<1e-3) {
        // uniform sampling
        Scalar delta0 = (_scaleMax-_scaleMin)/(_nbSamples-1);
        for(int i=1; i<_nbSamples; ++i) {
          m_scales[i] = _scaleMin + delta0*i;
        }

      } else {
        // log sampling
        double lMin = log(_scaleMin)/log(base);
        double lMax = log(_scaleMax)/log(base);
        Scalar delta = (lMax-lMin) / (_nbSamples-1);
        for(int i=0; i<_nbSamples; ++i) {
          m_scales[i] = std::pow(base, lMin + delta*i);
        }
      }
}


template<class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::poissonDiskSampling(int _nbPoints, Scalar* _coordinates, Scalar* _normals){
  
  printMsg("Start Poisson Disk Sampling Ponca");

  Timer tm{};

  if(m_scales.empty()){
    printErr("No scales detected");
    return;
  }
  m_trees.resize(m_scales.size());

  // if(!m_tree){
  //   printErr("KD-Tree not built");
  //   return;
  // }

  m_levels.resize(_nbPoints, -1);
  std::vector<bool> free(_nbPoints, true);

  // shuffle indices for random sampling
  std::vector<int> shuffledIndices(_nbPoints);
  std::iota(shuffledIndices.begin(), shuffledIndices.end(), 0);

  // shuffle them using the seed 0
  std::mt19937 random_engine{};
  random_engine.seed(0);
  // use the Fisher-Yates algorithm instead of std::shuffle, whose
  // results are platform-dependent
  ttk::shuffle(shuffledIndices, random_engine);

  for(int iScale = 0; iScale<m_scales.size(); iScale++){
    Scalar disk_radius = m_factor * m_scales[iScale];

    // build kd tree of previous level
    buildKDTree(iScale-1, _coordinates, _normals, _nbPoints);


    // random sampling on available points
    for(int iPoint=0; iPoint<_nbPoints; iPoint++){
      int iShuffledPoint = shuffledIndices[iPoint];
      m_levels[iShuffledPoint] = iScale; // sample this point
      auto tree_query = kdTree(iScale-1).range_neighbors(iShuffledPoint, disk_radius);
      for(int iNeigh : tree_query){
        free[iNeigh] = false; // mark point as unavailable for further sampling
      }
    }
  }

  }

  // build tree for last level
  buildKDTree(m_scales.size()-1, _coordinates, _normals, _nbPoints);


  printMsg("Poisson Disk Sampling done", 1, tm.getElapsedTime(), threadNumber_);
}


template <class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::buildKDTree(int level, Scalar* _coordinates, Scalar* _normals, int _nbPoints){


  Timer tm{};
  struct VN{
      Scalar* coordinates;
      Scalar* normals;
      int n;

      int size() const {
        return n;
      }
  };

  std::vector<int> sampling{};

  if(level==-1){
    sampling.resize(_nbPoints);
    std::iota(sampling.begin(), sampling.end(),0);
  }else{
    sampling.reserve(_nbPoints);
    for(int i = 0; i<_nbPoints; i++){
      if(is_sub_sampled(i,level)){
        sampling.push_back(i);
      }
    }
  }

  typedef typename Ponca::KdTree<Point>::PointContainer PtContainer;
  typedef typename Ponca::KdTree<Point>::IndexContainer IdxContainer;

  kdTree(level).buildWithSampling(VN({_coordinates, _normals, _nbPoints}), IdxContainer(sampling),
      [=](VN buffer, PtContainer &out){
                out.reserve(buffer.size());
                for(int i = 0; i<buffer.size(); i++){
                out.push_back(Point(buffer.coordinates, buffer.normals, i));
                }});
  printMsg("kd-tree built at level "+std::to_string(level) + ", point count "+std::to_string(kdTree(level).point_count()) + ", index count "+std::to_string(kdTree(level).sample_count()) , 1, tm.getElapsedTime());
}



template <class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::buildKnnGraph(){
  m_KnnGraph = new Ponca::KnnGraph<Point>(kdTree(-1),3);
}


template <class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::filterUnstablePoints(int level, Scalar* _coordinates, int _nbPoints, std::vector<SimplexId>& outputSubSampling, int* isStable){

  const std::vector<SimplexId>& currentSampling = kdTree(level).samples();
  std::vector<int> sampling{};
  SimplexId nbPointSampling = currentSampling.size();
  outputSubSampling.reserve(nbPointSampling);


  std::cout<<" Filtering unstable points"<<std::endl;
  sampling.reserve(nbPointSampling);
  for(SimplexId I = 0; I<nbPointSampling; I++){
    SimplexId i = currentSampling[I];
    if(is_sub_sampled(i,level) and isStable[I]){
      outputSubSampling.push_back(I);
    }
  }
}
