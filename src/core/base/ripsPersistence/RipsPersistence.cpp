#include "DataTypes.h"
#include "Debug.h"
#include "PCAnalysis.h"
#include "PCAnalysisUtils.h"
#include "Sampling.h"
#include <RipsPersistence.h>
#include <algorithm>
#include <boost/smart_ptr/detail/local_counted_base.hpp>
#include <cmath>
#include <functional>
#include <limits>
#include <sstream>
#include <algorithm>
#include <string>
#include <numeric>
#include <gmpxx.h>


using namespace ttk;

RipsPersistence::RipsPersistence() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("RipsPersistence");
  this->setDebugLevel(3);
}

// void RipsPersistence::computeEdgeVectorAtScaleDistanceMatrix(std::vector<std::vector<double>>& distMat, std::vector<SimplexId>& sampling_indexes, std::vector<SimplexId>& densitySampling, float* density){

//   Timer tm{};
//   printMsg("Computing exact list of edges from DISTANCE MATRIX, w/ density subsampling");

//   SimplexId nbPoints = densitySampling.size();
//   printMsg(std::to_string(nbPoints)+" points here");
//   printMsg("size of sampling vec : "+std::to_string(sampling_indexes.size()));

//   LongSimplexId size = (LongSimplexId)nbPoints*((LongSimplexId)nbPoints-1)/2;
//   edges_.clear();
//   edges_.resize(size);


// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(1)
// #endif // TTK_ENABLE_OPENMP
//   for(LongSimplexId I=0; I<(LongSimplexId)nbPoints-1; I++){
      
//       LongSimplexId i = densitySampling[I];


//     for(LongSimplexId J=I+1; J<(LongSimplexId)nbPoints; J++){
      
//       LongSimplexId j = densitySampling[J];


//       LongSimplexId id = J-1 + I*(2*nbPoints-I-3)/2;

//       // float minD = std::min(density[sampling_indexes[i]],density[sampling_indexes[j]]);
//       // float deltaD = std::abs(density[sampling_indexes[i]]-density[sampling_indexes[j]]);
//       // float maxD = std::max(density[sampling_indexes[i]],density[sampling_indexes[j]]);

//       // if(minD>1e12){
//       //   minD = 1;
//       //   maxD = 1;
//       //   deltaD=1;
//       // }
//       double dist = distMat[i][j];
//       if(std::isnan(dist)){
//         std::cout<<"dist is nan"<<std::endl;
//       }
//       // if(deltaD>10){
//       //   dist*=deltaD;
//       // }
//       // std::cout<<" dist is  "<<dist<<std::endl;
//       edges_[id] = std::make_tuple(I, J, dist);
//     }
//   }

//   std::cout<<" Base Method: "<<edges_.size()<<" edges in "<<tm.getElapsedTime()<<" s."<<std::endl;

//   const auto cmp = [=](const edgeTuple &e0, const edgeTuple &e1){
//     return compareEdgeTuple(e0, e1);
//   };

//   TTK_PSORT(this->threadNumber_, edges_.begin(), edges_.end(), cmp);
//   printMsg("Done, found ");
// }

// void RipsPersistence::computeEdgeVectorND(SimplexId nbPoints, float* coordinates, int dim){

//   Timer tm{};

//   LongSimplexId size = (LongSimplexId)nbPoints*((LongSimplexId)nbPoints-1)/2;
//   edges_.clear();
//   edges_.resize(size);

//   printMsg("Computing exact list of edges in "+std::to_string(dim)+"D");

// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(1)
// #endif // TTK_ENABLE_OPENMP
//   for(LongSimplexId i=0; i<(LongSimplexId)nbPoints-1; i++){

//       double ni = 0;
//       for(int k=0; k<dim; k++){
//         ni += coordinates[dim*i+k]*coordinates[dim*i+k];
//       }


//     for(LongSimplexId j=i+1; j<(LongSimplexId)nbPoints; j++){


//       // index in the upper triangular distance matrix
//       LongSimplexId id = j-1 + i*(2*nbPoints-i-3)/2;

//       double distCoords2 = 0; // ((xj-xi)*(xj-xi) + (yj-yi)*(yj-yi) + (zj-zi)*(zj-zi));
//       for(int k=0; k<dim; k++){
//         distCoords2 += (coordinates[dim*i+k]-coordinates[dim*j+k])*(coordinates[dim*i+k]-coordinates[dim*j+k]);
//       }
//       double dist = std::sqrt(distCoords2);

//       double nj=0;
//       for(int k=0; k<dim; k++){
//         nj += coordinates[dim*j+k]*coordinates[dim*j+k];
//       }

//       double dist2=0;
//       for(int k=0; k<dim; k++){
//         dist2 += coordinates[dim*i+k]*coordinates[dim*j+k];
//       }
//       double dist3=dist2*dist2/(ni*nj);
//       double dist4 = dist2;
//       // std::sqrt(dist3);
//       if(dist4>1){
//         dist4=1;
//       }else if(dist4<-1){
//         dist4=-1;
//       }
//       double dist5 = std::acos(dist4);

//       // std::cout<<"distances: proj "<<dist5<<" , euclidean "<<dist<<" "<<ni<<" "<<nj<<std::endl;
//       if(useProjectiveNorm_){
//         dist = dist5;
//       }
//       // if(std::isnan(dist5)){
//       //   std::cout<<"BIG NAN dist "<<dist<<" "<<ni<<" "<<nj<<std::endl;
//       //   std::cout<<"        dists "<<dist<<" "<<dist2<<" "<<dist3<<" "<<dist4<<" "<<dist5<<std::endl;
//       // }

//       if((ni<=1e-6 or nj<=1e-6) and useProjectiveNorm_){
//       }else{
//         edges_[id] = std::make_tuple(i, j, dist);
//       }
//     }
//   }

//   std::cout<<" Base Method: "<<edges_.size()<<" edges in "<<tm.getElapsedTime()<<" s."<<std::endl;

//   const auto cmp = [=](const edgeTuple &e0, const edgeTuple &e1){
//     return compareEdgeTuple(e0, e1);
//   };

//   TTK_PSORT(this->threadNumber_, edges_.begin(), edges_.end(), cmp);
//   printMsg("Done, found ");
// }
// void RipsPersistence::computeEdgeVector(SimplexId nbPoints, float* coordinates, float* normals, double ratio){

//   Timer tm{};

//   LongSimplexId size = (LongSimplexId)nbPoints*((LongSimplexId)nbPoints-1)/2;
//   edges_.clear();
//   edges_.resize(size);

//   printMsg("Computing list of edges, with ratio "+std::to_string(ratio));

// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif // TTK_ENABLE_OPENMP
//   for(LongSimplexId i=0; i<(LongSimplexId)nbPoints-1; i++){
//     double xi, yi, zi, nxi, nyi, nzi;
//     xi = coordinates[3*i];
//     yi = coordinates[3*i+1];
//     zi = coordinates[3*i+2];
//     nxi = normals[3*i];
//     nyi = normals[3*i+1];
//     nzi = normals[3*i+2];


//     for(LongSimplexId j=i+1; j<(LongSimplexId)nbPoints; j++){

//       double xj, yj, zj, nxj, nyj, nzj;
//       xj = coordinates[3*j];
//       yj = coordinates[3*j+1];
//       zj = coordinates[3*j+2];
//       nxj = normals[3*j];
//       nyj = normals[3*j+1];
//       nzj = normals[3*j+2];

//       // index in the upper triangular distance matrix
//       LongSimplexId id = j-1 + i*(2*nbPoints-i-3)/2;

//       double distCoords2 = ((xj-xi)*(xj-xi) + (yj-yi)*(yj-yi) + (zj-zi)*(zj-zi));
//       double distNorm2 = ((nxj-nxi)*(nxj-nxi) + (nyj-nyi)*(nyj-nyi) + (nzj-nzi)*(nzj-nzi));
//       double dist = std::sqrt((distCoords2+ratio*ratio*distNorm2)/((1+ratio)*(1+ratio)));

//       // if(1){
//       //   printMsg("Distance exact");
//       // std::cout<<" distance "<<distCoords2+ratio*distNorm2<<" "<<distCoords2<<" "<<distNorm2<<" "<<ratio*distNorm2<<std::endl;
//       // std::cout<<"   coords "<<i<<" "<<xi<<" "<<yi<<" "<<zi<<" "<<nxi<<" "<<nyi<<" "<<nzi<<std::endl;
//       // std::cout<<"   coords "<<j<<" "<<xj<<" "<<yj<<" "<<zj<<" "<<nxj<<" "<<nyj<<" "<<nzj<<std::endl;
//       // }

//       edges_[id] = std::make_tuple(i, j, dist);
//     }
//   }

//   std::cout<<" Base Method: "<<edges_.size()<<" edges in "<<tm.getElapsedTime()<<" s."<<std::endl;

//   const auto cmp = [=](const edgeTuple &e0, const edgeTuple &e1){
//     return compareEdgeTuple(e0, e1);
//   };

//   TTK_PSORT(this->threadNumber_, edges_.begin(), edges_.end(), cmp);
//   printMsg("Done, found ");
// }

    void RipsPersistence::sortPersistenceDiagram(SimplexId nbPoints){

      auto cmp = [=](const PersistencePair &a, const PersistencePair &b) {
        double pa = a.death.sfValue - a.birth.sfValue;
        double pb = b.death.sfValue - b.birth.sfValue;

        if(pa!=pb){
          return pa>pb;
        }

        if(a.birth.sfValue!=a.birth.sfValue){
          return a.birth.sfValue < b.birth.sfValue;
        }

        if(a.death.id != b.death.id){
          return !compareEdgeTuple(edges_[a.death.id], edges_[b.death.id]);
        }

        return a.birth.id < b.birth.id;

        // if(pa > pb){
        //   return true;
        // }
        // if(pa < pb){
        //   return true;
        // }
        // if(a.death.sfValue > b.death.sfValue){
        //   return true;
        // }
        // if(a.death.sfValue < b.death.sfValue){
        //   return false;
        // }
        // if(a.cc_size>b.cc_size){
        //   return true;
        // }
        // if(a.cc_size<b.cc_size){
        //   return false;
        // }
        // // if(!compareEdgeTuple(edges_[a.death.id], edges_[b.death.id]) and !compareEdgeTuple(edges_[b.death.id], edges_[a.death.id])){
        //   // std::cout<<"BIG BIG ISSUE"<<std::endl;
        //   // auto& e0 = edges_[a.death.id];
        //   // auto& e1 = edges_[b.death.id];
        //   // std::cout<<" e0 "<<std::get<0>(e0)<<" "<<std::get<1>(e0)<<" "<<std::get<2>(e0)<<std::endl;
        //   // std::cout<<" e1 "<<std::get<0>(e1)<<" "<<std::get<1>(e1)<<" "<<std::get<2>(e1)<<std::endl;
        //   // if(std::isnan(std::get<2>(e1))){
        //   //   std::cout<<"detected nan"<<std::endl;
        //   // }
        // // }
        // if(a.death.id == b.death.id){
        //   return a.birth.id<b.birth.id;
        // }
        // return !compareEdgeTuple(edges_[a.death.id], edges_[b.death.id]);
      };


        std::cout<<"test0, diagram size "<<diagram_.size()<<std::endl;

    std::sort(diagram_.begin(), diagram_.end(), cmp);

    // diagram_[0].death.sfValue = 1.1*diagram_[1].death.sfValue;

        std::cout<<"test1"<<std::endl;

  SimplexId size = nbPoints>0 ? nbPoints : diagram_.size();
  std::vector<SimplexId> goodIndex(size);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < diagram_.size(); ++i) {
    if(diagram_[i].birth.id>=goodIndex.size()){
      std::cout<<"error size"<<std::endl;
    }
    goodIndex[diagram_[i].birth.id] = i;
  }

        std::cout<<"test2"<<std::endl;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < diagram_.size(); ++i) {
    // auto& cid = diagram_[i].cc_parent;
    diagram_[i].cc_parent = goodIndex[diagram_[i].cc_parent];
  }
        // std::cout<<"test3"<<std::endl;
  // edges_.clear();
  // std::vector<edgeTuple>(edges_).swap(edges_);
  std::cout<<"sorted persistence diagram"<<std::endl;
    }

void RipsPersistence::setInitialSizesAndRepresentants(std::vector<SimplexId>& sizes,
                                                      std::vector<SimplexId>& reps){
  if(sizes.size()!=reps.size()){
    printErr("Sizes Mismatch in sizes and reps init");
    return;
  }
  sizes_ = std::move(sizes);
  reps_ = std::move(reps);
}
void RipsPersistence::setInitialSizesCopy(const std::vector<SimplexId>& sizes){
  sizes_ = sizes;
}
void RipsPersistence::setInitialSizes(std::vector<SimplexId>& sizes){
  sizes_ = std::move(sizes);
}
void RipsPersistence::setEdges(std::vector<edgeTuple>& edges){
  edges_ = std::move(edges);
}

    void RipsPersistence::computePersistenceDiagram(double maxEdgeLength,
                                                    SimplexId nbPoints){
			std::vector<edgeTuple> dummyVec{};
			computePersistenceDiagram(maxEdgeLength, nbPoints, dummyVec);
			dummyVec.clear();
		}

    void RipsPersistence::computePersistenceDiagramWithBirth(double maxEdgeLength,
                                                    SimplexId nbPoints,
                                                    std::vector<edgeTuple>& mergingEdges,
                                                    const std::vector<double>& point_births){

      std::cout<<"compute persistence diagram with births"<<std::endl;
      printMsg("Computing H0 persistence, "+std::to_string(nbPoints)+" points");
      std::cout<<"number of edges: "<<edges_.size()<<std::endl;
      std::cout<<"max Edge Length "<<maxEdgeLength<<std::endl;
      diagram_.clear();

      if(sizes_.size()!=nbPoints){
        printMsg("Initializing sizes");
        sizes_.resize(nbPoints,1);
        std::fill(sizes_.begin(), sizes_.end(), 1);
      }else{
        printMsg("Using precomputed inital sizes");

      }
      if(reps_.size()!=nbPoints){
        printMsg("Initializing representants");
        reps_.resize(nbPoints);
        std::iota(reps_.begin(), reps_.end(), 0);
      }else{
        printMsg("Using precomputed inital representants");
      }

      Timer tm{};

      // get representative of current extremum
      const auto getRep = [&](SimplexId v) -> SimplexId {
        auto r = reps_[v];
        while(r != v) {
          v = r;
          r = reps_[v];
        }
        return r;
      };

      // merge connected components of two extrema
      const auto merge = [&](SimplexId cc0, SimplexId cc1) -> SimplexId {
       auto r0 = getRep(cc0);
       auto r1 = getRep(cc1);

       if(r0!=r1){
        if(point_births[r0]<point_births[r1]){
          std::swap(r0, r1);
        }
        // r0 dies at the benefit of r1
        reps_[r0]=r1;
        sizes_[r1]+=sizes_[r0];

        // return the birth id
        return r0;
       }
       return -1;
      };
      const auto merge_smallest = [&](SimplexId cc0, SimplexId cc1) -> SimplexId {
       auto r0 = getRep(cc0);
       auto r1 = getRep(cc1);

       if(r0!=r1){
        if(point_births[r0]>point_births[r1]){
          std::swap(r0, r1);
        }
        // r0 dies at the benefit of r1
        reps_[r0]=r1;
        sizes_[r1]+=sizes_[r0];

        // return the birth id
        return r0;
       }
       return -1;
      };

      const auto addPair = [&](SimplexId birthId, SimplexId edgeId, double radius, bool isFinite){
        // if(point_births[birthId]>radius){
          // std::cout<<"THIS SHOULD NOT HAPPEN"<<std::endl;
          // std::cout<<point_births[birthId]<<" "<<radius<<std::endl;
          // std::cout<<birthId<<" "<<edgeId<<std::endl;
          // std::cout<<" "<<std::get<0>(edges_[edgeId])<<" "<<std::get<1>(edges_[edgeId])<<" "<<std::get<2>(edges_[edgeId])<<std::endl;
        // }
          diagram_.emplace_back(PersistencePair{
              CriticalVertex{birthId, CriticalType::Local_minimum, point_births[birthId],{}},
              CriticalVertex{edgeId, CriticalType::Saddle1, radius,{}},
              0, isFinite, sizes_[birthId], reps_[birthId],{}});
      };

      mergingEdges.clear();

      // sequential loop 
      for(SimplexId i=0; i<edges_.size(); i++){
        auto e = edges_[i];
        double edgeLength = std::get<2>(e);
        // if(edgeLength>maxEdgeLength and maxEdgeLength>0){
        //   break;
        // }
        SimplexId c0 = std::get<0>(e);
        SimplexId c1 = std::get<1>(e);

        SimplexId dyingCC = merge(c0,c1);
        if(dyingCC != -1){
         addPair(dyingCC,i, edgeLength, true);
         mergingEdges.push_back(e);
        }
      }


      std::vector<SimplexId> toMerge{};
      for(SimplexId i=0; i<nbPoints; i++){
        if(reps_[i]==i){
          toMerge.push_back(i);
        }
      }

      std::cout<<"FOUND "<<toMerge.size()<<" connected components"<<std::endl;

      for(auto& cc : toMerge){
        SimplexId deathId{};
        double death{};
        if(edges_.empty()){
          printWrn("No edges");
          deathId = -cc;
          death = 42;
        }else{
          deathId =  edges_.size()-1;
          death = std::get<2>(edges_.back())+1e-5;
        }
        addPair(cc, deathId, death, false);
        // std::cout<<" add infinite pair cc "<<cc<<" size "<<sizes_[cc]<<std::endl;
      }

      printMsg("Found "+std::to_string(diagram_.size())+" persistence pairs",1,tm.getElapsedTime(),threadNumber_);
    }
    void RipsPersistence::computePersistenceDiagram(double maxEdgeLength,
                                                    SimplexId nbPoints,
                                                    std::vector<edgeTuple>& mergingEdges){

      std::cout<<"compute persistence diagram"<<std::endl;
      printMsg("Computing H0 persistence, "+std::to_string(nbPoints)+" points");
      std::cout<<"number of edges: "<<edges_.size()<<std::endl;
      std::cout<<"max Edge Length "<<maxEdgeLength<<std::endl;
      diagram_.clear();

      if(sizes_.size()!=nbPoints){
        printMsg("Initializing sizes");
        sizes_.resize(nbPoints,1);
        std::fill(sizes_.begin(), sizes_.end(), 1);
      }else{
        printMsg("Using precomputed inital sizes");

      }
      if(reps_.size()!=nbPoints){
        printMsg("Initializing representants");
        reps_.resize(nbPoints);
        std::iota(reps_.begin(), reps_.end(), 0);
      }else{
        printMsg("Using precomputed inital representants");
      }

      Timer tm{};

      // get representative of current extremum
      const auto getRep = [&](SimplexId v) -> SimplexId {
        auto r = reps_[v];
        while(r != v) {
          v = r;
          r = reps_[v];
        }
        return r;
      };

      // merge connected components of two extrema
      const auto merge = [&](SimplexId cc0, SimplexId cc1) -> SimplexId {
       auto r0 = getRep(cc0);
       auto r1 = getRep(cc1);

       if(r0!=r1){
        if(sizes_[r0]>sizes_[r1]){
          std::swap(r0, r1);
        }
        // r0 dies at the benefit of r1
        reps_[r0]=r1;
        sizes_[r1]+=sizes_[r0];

        // return the birth id
        return r0;
       }
       return -1;
      };
      const auto merge_smallest = [&](SimplexId cc0, SimplexId cc1) -> SimplexId {
       auto r0 = getRep(cc0);
       auto r1 = getRep(cc1);

       if(r0!=r1){
        if(sizes_[r0]<sizes_[r1]){
          std::swap(r0, r1);
        }
        // r0 dies at the benefit of r1
        reps_[r0]=r1;
        sizes_[r1]+=sizes_[r0];

        // return the birth id
        return r0;
       }
       return -1;
      };

      const auto addPair = [&](SimplexId birthId, SimplexId edgeId, double radius, bool isFinite){
        // std::cout<<"add pair with radius "<<radius<<std::endl;
        diagram_.emplace_back(PersistencePair{
            CriticalVertex{birthId, CriticalType::Local_minimum, 0,{}},
            CriticalVertex{edgeId, CriticalType::Local_maximum, radius,{}},
            0, isFinite, sizes_[birthId], reps_[birthId],{}});
      };

      mergingEdges.clear();

      // sequential loop 
      for(SimplexId i=0; i<edges_.size(); i++){
        auto e = edges_[i];
        double edgeLength = std::get<2>(e);
        if(edgeLength>maxEdgeLength and maxEdgeLength>0){
          break;
        }
        SimplexId c0 = std::get<0>(e);
        SimplexId c1 = std::get<1>(e);

        SimplexId dyingCC = merge(c0,c1);
        if(dyingCC != -1){
         addPair(dyingCC,i, edgeLength/2, true);
         mergingEdges.push_back(e);
        }
      }

      // merge non-connected components
      // SimplexId biggest = getRep(0);
      // SimplexId smallest = getRep(0);
      std::vector<SimplexId> toMerge{};
      for(SimplexId i=0; i<nbPoints; i++){
        if(reps_[i]==i){
          // std::cout<<"FOUND non connected component "<<i<<" of size "<<sizes_[i]<<std::endl;
          // if(sizes_[i]>sizes_[biggest]){
          //   biggest = getRep(i);
          // }
          // if(sizes_[i]<sizes_[smallest]){
          //   smallest = getRep(i);
          // }
          toMerge.push_back(i);
          // std::cout<<"     merge with "<<getRep(0)<<std::endl;
          // nbMerges++;
          // addPair(getRep(i), edges_.size()-nbMerges, std::get<2>(edges_[edges_.size()-nbMerges]),false);
          // SimplexId dyingCC = merge(i,getRep(0));
          // std::cout<<"     dying is "<<dyingCC<<std::endl;
          // if(dyingCC!=-1){
          //   addPair(dyingCC, edges_.size()-nbMerges, std::get<2>(edges_[edges_.size()-nbMerges]),false);
          // }else{
          // }
        }
      }

      int nbMerges=1;
      std::cout<<"FOUND "<<toMerge.size()<<" non connected components"<<std::endl;
      // for(SimplexId cc : toMerge){
      //   // if(cc == smallest){
      //   //   continue;
      //   // }
      //   // if(1){
      //   //   // merge to biggest
      //   //   std::cout<<"MERGING non connected component "<<cc<<" of size "<<sizes_[cc]<<std::endl;
      //   //   std::cout<<"     merge with "<<biggest<<" of size "<<sizes_[biggest]<<std::endl;
      //   //   SimplexId dyingCC = merge(cc,biggest);
      //   //   // std::cout<<"     dying is "<<dyingCC<<std::endl;
      //   //   if(dyingCC!=-1){
      //   //     nbMerges++;
      //   //     double pers = sizes_[dyingCC]<=1 ? 0 : std::get<2>(edges_[edges_.size()-nbMerges]);
      //   //     addPair(dyingCC, edges_.size()-nbMerges, pers,true);

      //   //   }
      //   // }else{
      //   //   // merge to smallest 
      //   //   std::cout<<"MERGING non connected component "<<cc<<" of size "<<sizes_[cc]<<std::endl;
      //   //   std::cout<<"     merge with "<<smallest<<" of size "<<sizes_[smallest]<<std::endl;
      //   //   SimplexId dyingCC = merge(cc,smallest);
      //   //   // std::cout<<"     dying is "<<dyingCC<<std::endl;
      //   //   if(dyingCC!=-1){
      //   //     nbMerges++;
      //   //     double pers = sizes_[dyingCC]<=1 ? 0 : std::get<2>(edges_[edges_.size()-nbMerges]);
      //   //     addPair(dyingCC, edges_.size()-nbMerges, pers, false);
      //   //   }
      //   // }
      // }

      // merge smallest and biggest
      // merge(biggest,smallest);
      // nbMerges++;
      // addPair(smallest, edges_.size()-nbMerges, std::get<2>(edges_[edges_.size()-nbMerges]), false);

      // Add last global pair
      // addPair(getRep(biggest),edges_.size()-1,std::get<2>(edges_.back())/2*10,false);
      //
      for(auto& cc : toMerge){
        SimplexId deathId{};
        double death{};
        if(edges_.empty()){
          printWrn("No edges");
          deathId = -cc;
          death = 42;
        }else{
          deathId =  edges_.size()-1;
          death = std::get<2>(edges_.back())+1e-5;
        }
        addPair(cc, deathId, death, false);
        // std::cout<<" add infinite pair cc "<<cc<<" size "<<sizes_[cc]<<std::endl;
      }

      printMsg("Found "+std::to_string(diagram_.size())+" persistence pairs",1,tm.getElapsedTime(),threadNumber_);
    }

// void RipsPersistence::computeEdgeVectorWithSamplingNF(ScaleSpaceSampling<pointCloud::MyPoint6D::Scalar, pointCloud::MyPoint6D>& sampler,
//                                                     SimplexId nbPoints, float* coordinates, float* normals, double ratio){
//   printMsg("Compute Rips edges with 6D sampling - NanoFlann");
//   Timer tm{};

//   sampler.printTreesInfosNF();
//   // auto pt = sampler.getPTest();
//   edges_.clear();
//   int nScale = sampler.scale_count();
//   for(int iScale=nScale-1; 0<iScale; iScale--){

//     LongSimplexId nbPointsAtScale = sampler.sampling(iScale).size();
//     double distMax = iScale == nScale-1 ? std::numeric_limits<double>::max() : sampler.scale(iScale+1)*sampler.samplingFactor();
//     auto ptree = sampler.kdTreeNF(iScale);
//     // sampler.printTreesInfosNF();
//     Timer tm_scale{};

//   for(LongSimplexId id=0; id<sampler.sampling(iScale).size(); id++){


//     LongSimplexId i = sampler.sampling(iScale)[id];
//     double xi, yi, zi, nxi, nyi, nzi;
//     xi = coordinates[3*i];
//     yi = coordinates[3*i+1];
//     zi = coordinates[3*i+2];
//     nxi = normals[3*i];
//     nyi = normals[3*i+1];
//     nzi = normals[3*i+2];


//     std::vector<nanoflann::ResultItem<unsigned int, float>> ret_matches;
//     const float queryPoint[6] = {(float)xi,(float)yi,(float)zi,(float)nxi,(float)nyi,(float)nzi};
//     SimplexId nMatches = ptree->radiusSearch(queryPoint, distMax*distMax, ret_matches);

//     for(auto r : ret_matches){
//       LongSimplexId j = sampler.sampling(iScale)[r.first];

//       if(i<j){
//         edges_.emplace_back(i, j, r.second);
//       }
//     }

//     int d_progress = (int)nbPointsAtScale*0.03+1;
//     if(id%d_progress==0){ // write progress
//       std::stringstream msg;
//       msg<<"Scale "<<iScale<<" ("<<nbPointsAtScale<<", "<<distMax<<"), size "<<ptree->size_;
//       printMsg(msg.str(), (double)id/nbPointsAtScale, tm_scale.getElapsedTime(), threadNumber_, debug::LineMode::REPLACE);

//     }
//   }
//     std::stringstream msg;
//     msg<<"Scale "<<iScale<<" ("<<nbPointsAtScale<<", "<<distMax<<"), size "<<ptree->size_;
//     printMsg(msg.str(), 1, tm_scale.getElapsedTime(), threadNumber_);
//   }

//   if(1){
//     LongSimplexId nbPointsAtScale = nbPoints;
//     double distMax = nScale<1 ? std::numeric_limits<double>::max() : sampler.scale(1)*sampler.samplingFactor();
//     auto ptree = sampler.kdTreeNF(0);

//     Timer tm_scale;
//     // sampler.printTreesInfosNF();
//     // std::cout<<"Last scale "<<0<<" ("<<nbPointsAtScale<<"points, dist "<<distMax<<"), size "<<ptree->size_<<" address "<<ptree<<std::endl;

//   for(LongSimplexId i=0; i<nbPointsAtScale; i++){

//     double xi, yi, zi, nxi, nyi, nzi;
//     xi = coordinates[3*i];
//     yi = coordinates[3*i+1];
//     zi = coordinates[3*i+2];
//     nxi = normals[3*i];
//     nyi = normals[3*i+1];
//     nzi = normals[3*i+2];


//     std::vector<nanoflann::ResultItem<unsigned int, float>> ret_matches;
//     const float queryPoint[6] = {(float)xi,(float)yi,(float)zi,(float)nxi,(float)nyi,(float)nzi};
//     SimplexId nMatches = ptree->radiusSearch(queryPoint, distMax*distMax, ret_matches);

//     bool isFull=false;
//     for(auto r : ret_matches){
//       LongSimplexId j = r.first;

//       if(i<j){
//         edges_.emplace_back(i, j, r.second);
//         isFull = edges_.size()==edges_.max_size();
//         if(isFull) break;
//       }
//     }
//     if(isFull){
//       printWrn("Max size reached, stopping adding edges");
//       break;
//     }
//     int d_progress = (int)nbPointsAtScale*0.03+1;
//     if(i%d_progress==0){ // write progress
//       std::stringstream msg;
//       msg<<"Scale "<<0<<" ("<<nbPointsAtScale<<", "<<distMax<<"), size "<<ptree->size_;
//       printMsg(msg.str(), (double)i/nbPointsAtScale, tm_scale.getElapsedTime(), threadNumber_, debug::LineMode::REPLACE);

//     }
//   }

//     std::stringstream msg;
//     msg<<"Scale "<<0<<" ("<<nbPointsAtScale<<", "<<distMax<<"), size "<<ptree->size_;
//     printMsg(msg.str(), 1, tm_scale.getElapsedTime(), threadNumber_);
//   }

//   std::stringstream msg;
//   msg<<"Hierarchical Method NF : "<<edges_.size()<<" edges";
//   printMsg(msg.str(), 1, tm.getElapsedTime(), threadNumber_);

//   const auto cmp = [=](const edgeTuple &e0, const edgeTuple &e1){
//     return compareEdgeTuple(e0, e1);
//   };
//   //   return (std::get<2>(e0) < std::get<2>(e1)) || ( (std::get<2>(e0) == std::get<2>(e1)) &&
//   //                           ( (std::get<0>(e0) < std::get<0>(e1)) || ((std::get<0>(e0) == std::get<0>(e1))&&(std::get<1>(e0) < std::get<1>(e1))) )
//   //                           );
//   // };

//   TTK_PSORT(this->threadNumber_, edges_.begin(), edges_.end(), cmp);
//   printMsg("Done, found ");

// }

// void RipsPersistence::computeEdgeVectorWithSampling(ScaleSpaceSampling<pointCloud::MyPoint6D::Scalar, pointCloud::MyPoint6D>& sampler,
//                                                     SimplexId nbPoints, float* coordinates, float* normals, double ratio){

//   printMsg("Compute Rips edges with 6D sampling");
//   Timer tm{};

//   // sampler.printTreesInfosNF();
//   // auto pt = sampler.getPTest();
//   edges_.clear();
//   int nScale = sampler.scale_count();
//   for(int iScale=nScale-1; 0<iScale; iScale--){

//     double distMax = iScale == nScale-1 ? std::numeric_limits<double>::infinity() : sampler.scale(iScale+1)*sampler.samplingFactor();
//     auto& tree = sampler.kdTree(iScale);

//     Timer tm_scale{};

//     LongSimplexId nbPointsAtScale = tree.samples().size();

//     // std::cout<<" points in here at scale "<<iScale<<": "<<tree.sample_count()<<" "<<tree.samples().size()<<std::endl;
//   // for(LongSimplexId id=0; id<nbPointsAtScale; id++){
//   //   LongSimplexId i = tree.samples()[id];
//   //   std::cout<<" "<<i;
//   // }
//   // std::cout<<std::endl;
//     // if(iScale==nScale-1){
//     //   for(LongSimplexId id=0; id<nbPointsAtScale; id++){
//     //     LongSimplexId i = tree.samples()[id];
//     //     double xi, yi, zi, nxi, nyi, nzi;
//     //     xi = coordinates[3*i];
//     //     yi = coordinates[3*i+1];
//     //     zi = coordinates[3*i+2];
//     //     nxi = normals[3*i];
//     //     nyi = normals[3*i+1];
//     //     nzi = normals[3*i+2];

//     //     for(LongSimplexId jd=0; jd<nbPointsAtScale; jd++){
//     //       LongSimplexId j = tree.samples()[jd];
//     //       double xj, yj, zj, nxj, nyj, nzj;
//     //       xj = coordinates[3*j];
//     //       yj = coordinates[3*j+1];
//     //       zj = coordinates[3*j+2];
//     //       nxj = normals[3*j];
//     //       nyj = normals[3*j+1];
//     //       nzj = normals[3*j+2];


//     //       double distCoords2 = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi) + (zj-zi)*(zj-zi);
//     //       double distNorm2 = (nxj-nxi)*(nxj-nxi) + (nyj-nyi)*(nyj-nyi) + (nzj-nzi)*(nzj-nzi);
//     //       double dist = std::sqrt((distCoords2 + distNorm2*ratio*ratio)/((1+ratio)*(1+ratio)));

//     //       if(i<j){
//     //         edges_.emplace_back(i, j, dist);
//     //         // std::cout<<"ADDED "<<std::endl;
//     //         //   std::cout<<" distance "<<distCoords2+ratio*distNorm2<<" "<<distCoords2<<" "<<distNorm2<<" "<<ratio*distNorm2<<std::endl;
//     //         //   std::cout<<"   coords "<<i<<" "<<xi<<" "<<yi<<" "<<zi<<" "<<nxi<<" "<<nyi<<" "<<nzi<<std::endl;
//     //         //   std::cout<<"   coords "<<j<<" "<<xj<<" "<<yj<<" "<<zj<<" "<<nxj<<" "<<nyj<<" "<<nzj<<std::endl;
//     //       }
//     //     }
//     //   }
//     //   continue;
//     // }

//   for(LongSimplexId id=0; id<nbPointsAtScale; id++){

//     LongSimplexId i = tree.samples()[id];

//     double xi, yi, zi, nxi, nyi, nzi;
//     xi = coordinates[3*i];
//     yi = coordinates[3*i+1];
//     zi = coordinates[3*i+2];
//     nxi = normals[3*i];
//     nyi = normals[3*i+1];
//     nzi = normals[3*i+2];



//     auto tree_query =tree.range_neighbors(i,distMax);

//     // int nMatches=0;
//     // for(LongSimplexId j : tree_query){nMatches++;}
//     // std::cout<<"query from "<<i<<" : "<<nMatches<<" matches"<<std::endl;

//     for(LongSimplexId j : tree_query){
//       // LongSimplexId j = tree.samples()[jIndex];
//       // std::cout<<"   - "<<j<<std::endl;
//       double xj, yj, zj, nxj, nyj, nzj;
//       xj = coordinates[3*j];
//       yj = coordinates[3*j+1];
//       zj = coordinates[3*j+2];
//       nxj = normals[3*j];
//       nyj = normals[3*j+1];
//       nzj = normals[3*j+2];


//       double distCoords2 = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi) + (zj-zi)*(zj-zi);
//       double distNorm2 = (nxj-nxi)*(nxj-nxi) + (nyj-nyi)*(nyj-nyi) + (nzj-nzi)*(nzj-nzi);
//       double dist = std::sqrt((distCoords2 + distNorm2*ratio*ratio)/((1+ratio)*(1+ratio)));

//       if(i<j){
//         edges_.emplace_back(i, j, dist);
//       }
//     }

//     int d_progress = (int)nbPointsAtScale*0.03+1;
//     if(id%d_progress==0){ // write progress
//       std::stringstream msg;
//       msg<<"Scale "<<iScale<<" ("<<nbPointsAtScale<<", "<<distMax<<"), size "<<tree.sample_count();
//       printMsg(msg.str(), (double)id/nbPointsAtScale, tm_scale.getElapsedTime(), threadNumber_, debug::LineMode::REPLACE);

//     }
//   }
//     std::stringstream msg;
//     msg<<"Scale "<<iScale<<" ("<<nbPointsAtScale<<", "<<distMax<<"), size "<<tree.sample_count();
//     printMsg(msg.str(), 1, tm_scale.getElapsedTime(), threadNumber_);
//   }

//   if(1){

//     double distMax =  nScale<1 ? std::numeric_limits<double>::infinity() : sampler.scale(1)*sampler.samplingFactor()*1.01;
//     auto& tree = sampler.kdTree(0);

//     Timer tm_scale{};

//     LongSimplexId nbPointsAtScale = tree.samples().size();

//   for(LongSimplexId id=0; id<nbPointsAtScale; id++){

//     LongSimplexId i = tree.samples()[id];

//     double xi, yi, zi, nxi, nyi, nzi;
//     xi = coordinates[3*i];
//     yi = coordinates[3*i+1];
//     zi = coordinates[3*i+2];
//     nxi = normals[3*i];
//     nyi = normals[3*i+1];
//     nzi = normals[3*i+2];


//     auto tree_query =tree.range_neighbors(i,distMax);

//       // int nMatches=0;
//     // for(LongSimplexId j : tree_query){nMatches++;}

//     for(LongSimplexId j : tree_query){
//       // LongSimplexId j = tree.samples()[jIndex];
//       double xj, yj, zj, nxj, nyj, nzj;
//       xj = coordinates[3*j];
//       yj = coordinates[3*j+1];
//       zj = coordinates[3*j+2];
//       nxj = normals[3*j];
//       nyj = normals[3*j+1];
//       nzj = normals[3*j+2];


//       double distCoords2 = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi) + (zj-zi)*(zj-zi);
//       double distNorm2 = (nxj-nxi)*(nxj-nxi) + (nyj-nyi)*(nyj-nyi) + (nzj-nzi)*(nzj-nzi);
//       double dist = std::sqrt((distCoords2 + distNorm2*ratio*ratio)/((1+ratio)*(1+ratio)));

//       if(i<j){
//         edges_.emplace_back(i, j, dist);
//       // if(1){
//       //   printWrn("Distance is weirdly big, might be a bug");
//       // std::cout<<" distance "<<distCoords2+ratio*distNorm2<<" "<<distCoords2<<" "<<distNorm2<<" "<<ratio*distNorm2<<std::endl;
//       // std::cout<<"   coords "<<i<<" "<<xi<<" "<<yi<<" "<<zi<<" "<<nxi<<" "<<nyi<<" "<<nzi<<std::endl;
//       // std::cout<<"   coords "<<j<<" "<<xj<<" "<<yj<<" "<<zj<<" "<<nxj<<" "<<nyj<<" "<<nzj<<std::endl;
//       // }
//       }
//     }

//     int d_progress = (int)nbPointsAtScale*0.03+1;
//     if(id%d_progress==0){ // write progress
//       std::stringstream msg;
//       msg<<"Scale "<<0<<" ("<<nbPointsAtScale<<", "<<distMax<<"), size "<<tree.sample_count();
//       printMsg(msg.str(), (double)id/nbPointsAtScale, tm_scale.getElapsedTime(), threadNumber_, debug::LineMode::REPLACE);

//     }
//   }
//       std::stringstream msg;
//       msg<<"Scale "<<0<<" ("<<nbPointsAtScale<<", "<<distMax<<"), size "<<tree.sample_count();
//       printMsg(msg.str(), 1, tm_scale.getElapsedTime(), threadNumber_);
//   }

//   std::stringstream msg;
//   msg<<"Hierarchical Method : "<<edges_.size()<<" edges";
//   printMsg(msg.str(), 1, tm.getElapsedTime(), threadNumber_);

//   const auto cmp = [=](const edgeTuple &e0, const edgeTuple &e1){
//     return compareEdgeTuple(e0, e1);
//   };

//   TTK_PSORT(this->threadNumber_, edges_.begin(), edges_.end(), cmp);
// }


bool RipsPersistence::compareEdgeTuple(const edgeTuple& e0, const edgeTuple& e1){
    return (std::get<2>(e0) < std::get<2>(e1)) || ( (std::get<2>(e0) == std::get<2>(e1)) &&
                            ( (std::get<0>(e0) < std::get<0>(e1)) || ((std::get<0>(e0) == std::get<0>(e1))&&(std::get<1>(e0) < std::get<1>(e1))) )
                            );
}

int RipsPersistence::findScaleId(SimplexId ripsBirthId, const std::vector<SimplexId>& pointsAtScaleStartIdx) const {
  for(int iScale = pointsAtScaleStartIdx.size()-1; iScale>=0; iScale--){
    if(ripsBirthId>=pointsAtScaleStartIdx[iScale]){
      return iScale;
    }
  }
  std::cout<<"ERROR Found no scale id from pointId"<<std::endl;
  return -1;
}



void RipsPersistence::compute3DEmbedding(SimplexId _nbPoints,
                  int _dim,
                  const std::vector<edgeTuple>& _edges,
                  std::vector<std::vector<double>>& _outputCoordinates) const {

  printMsg("Computing 3d Embedding",0,0, threadNumber_);


    size_t nb = _nbPoints;
    size_t nbnb = nb*nb;
    std::vector<double> distanceMatrix(nbnb,0);


#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(auto& e : _edges){
      // std::cout<<"edge"<<std::endl;
      size_t i = std::get<0>(e);
      size_t j = std::get<1>(e);
      double dist = std::get<2>(e);
      // std::cout<<"dist "<<dist<<std::endl;
      // std::cout<<"max size "<<distanceMatrix.max_size()<<std::endl;
      // std::cout<<"filling "<<i+_nbPoints*j<<" and "<<j+_nbPoints*i<<" outta "<<distanceMatrix.size()<<std::endl;
      // std::cout<<"i "<<i<<" j "<<j<<std::endl;
      // std::cout<<"nb "<<_nbPoints<<std::endl;
      distanceMatrix[(size_t)(i + nb*j)] = dist;
      distanceMatrix[(size_t)(j + nb*i)] = dist;
    }

    DimensionReduction dr{};
    auto method = DimensionReduction::METHOD::MDS;
    dr.setIsInputDistanceMatrix(true);
    dr.setInputMethod(method);
    dr.setDebugLevel(debugLevel_);
    if(method == DimensionReduction::METHOD::MDS){
      dr.setThreadNumber(1);
    }
    dr.setInputModulePath("/home/julius/ttk/core/base/dimensionReduction");
    dr.setInputNumberOfComponents(_dim);
    dr.execute(_outputCoordinates, distanceMatrix, _nbPoints, _nbPoints);


}

void RipsPersistence::replaceBirthIds(std::vector<SimplexId>& sampling_indexes){
  for(auto& pp : diagram_){
    pp.birth.id =  sampling_indexes[pp.birth.id];
  }
}

void RipsPersistence::sortEdgesAscending(){
  const auto cmp = [=](const edgeTuple &e0, const edgeTuple &e1){
    return compareEdgeTuple(e0, e1);
  };

  TTK_PSORT(this->threadNumber_, edges_.begin(), edges_.end(), cmp);
}

// void RipsPersistence::filterDistanceMatrixByDensity(std::vector<std::vector<double>>& distMat, float min_density, float* density, std::vector<SimplexId>& sampling_indexes){
  
//   SimplexId nbPoints = sampling_indexes.size();
//   for(int i=0; i<nbPoints; i++){
//     if(density[sampling_indexes[i]]<min_density){
//       for(int j=0; j<nbPoints; j++){
//         distMat[i][j] = (double)std::numeric_limits<double>::max();
//       }
//     }
//   }
// }

// void RipsPersistence::computeEdgeVectorMultiScaleND(const std::vector<std::vector<SimplexId>>& sampling,
//                                                     const std::vector<std::vector<double>>& fitCoordinates,
//                                                     const std::vector<std::vector<float>>& mls3DCoordinates,
//                                                     const std::vector<float>& pointCloudBarycenter,
//                                                     std::function<double(const double*, const double*, const float*, const float*, const float*)> distanceFunction,
//                                                     int dim){

//   Timer tm{};
//   printMsg("Computing exact list of edges in "+std::to_string(dim)+"D, WITH SUBSAMPLING, all scales mixed");

//   SimplexId nbScales = sampling.size();

//   edges_.clear();

//   // LongSimplexId nbEdges = nbEdgeMonoScale * nbScales + nbEdgeBtwScales*(nbScales-1);
//   LongSimplexId nbEdges = 0;
//   std::vector<SimplexId> edgesAtScaleStartIdx(nbScales,0);
//   std::vector<SimplexId> edgesBtwScalesStartIdx(nbScales,0);
//   std::vector<SimplexId> pointsAtScaleStartIdx(nbScales,0);

//   for (int iScale=0; iScale<nbScales; iScale++){
//     SimplexId nbPointsAtScale = sampling[iScale].size();
//     LongSimplexId nbEdgeAtScale = (LongSimplexId)nbPointsAtScale*((LongSimplexId)nbPointsAtScale-1)/2;
//     LongSimplexId nbEdgeWithNextScale = iScale<nbScales-1 ? (LongSimplexId)nbPointsAtScale*((LongSimplexId)sampling[iScale+1].size()) : 0;
//     std::cout<<" scale "<<iScale<<std::endl;
//     std::cout<<"   - "<<nbPointsAtScale<<" points"<<std::endl;
//     std::cout<<"   - "<<nbEdgeAtScale<<" edges intern"<<std::endl;
//     std::cout<<"   - "<<nbEdgeWithNextScale<<" edges w/ next scale"<<std::endl;

//     edgesBtwScalesStartIdx[iScale] = edgesAtScaleStartIdx[iScale] + nbEdgeAtScale;

//     if(iScale<nbScales-1){
//       edgesAtScaleStartIdx[iScale+1] = edgesBtwScalesStartIdx[iScale] + nbEdgeWithNextScale;
//       pointsAtScaleStartIdx[iScale+1] = pointsAtScaleStartIdx[iScale] + nbPointsAtScale;
//     }

//     nbEdges += nbEdgeAtScale + nbEdgeWithNextScale;
//   }

//   edges_.resize(nbEdges);


//     LongSimplexId numberOfEdges=0;
//   for (int iScale=0; iScale<nbScales; iScale++){

//     std::cout<<" Scale "<<iScale<< " : "<<std::endl;

//     SimplexId nbPointsAtScale = sampling[iScale].size();
//     LongSimplexId nbEdgeAtScale = (LongSimplexId)nbPointsAtScale*((LongSimplexId)nbPointsAtScale-1)/2;
//     LongSimplexId nbEdgeWithNextScale = iScale<nbScales-1 ? 0 : (LongSimplexId)nbPointsAtScale*((LongSimplexId)sampling[iScale+1].size());
//     // SimplexId nbPoints = sampling[0].size();
//     // LongSimplexId nbEdgeMonoScale = (LongSimplexId)nbPoints*((LongSimplexId)nbPoints-1)/2;
//     // LongSimplexId nbEdgeBtwScales = (LongSimplexId)nbPoints*((LongSimplexId)nbPoints);
//     // LongSimplexId nbEdges = nbEdgeMonoScale * nbScales + nbEdgeBtwScales*(nbScales-1);


//     // if(nbPoints!=sampling[iScale].size()){
//     //   printWrn("BIG ISSUE: different number of points in each scale");
//     // }

//     const std::vector<SimplexId>& sampling_indexes = sampling[iScale];
//     const double* coordinates = fitCoordinates[iScale].data();
//     const float* mlsCoordinates = mls3DCoordinates[iScale].data();
//     LongSimplexId numberEdgeAtScale=0;
//     LongSimplexId numberEdgeWithNextScale=0;

// // #ifdef TTK_ENABLE_OPENMP
// // #pragma omp parallel for num_threads(1)
// // #endif // TTK_ENABLE_OPENMP
// //
//     for(LongSimplexId I=0; I<(LongSimplexId)nbPointsAtScale-1; I++){

//       LongSimplexId i = sampling_indexes[I];

//       const double* coordsi = &coordinates[dim*i];

//       for(LongSimplexId J=I+1; J<(LongSimplexId)nbPointsAtScale; J++){

//         LongSimplexId j = sampling_indexes[J];

//         const double* coordsj = &coordinates[dim*j];

//         // index in the upper triangular distance matrix
//         LongSimplexId id = J-1 + I*(2*nbPointsAtScale-I-3)/2;

//         double dist = distanceFunction(coordsi, coordsj, &mlsCoordinates[3*i], &mlsCoordinates[3*j], pointCloudBarycenter.data());

//         edges_[edgesAtScaleStartIdx[iScale]+id] = std::make_tuple(I+pointsAtScaleStartIdx[iScale], J+pointsAtScaleStartIdx[iScale], dist);
//         numberEdgeAtScale++;
//         numberOfEdges++;
//       }
//     }

//     if(iScale<nbScales-1){ 

//       const std::vector<SimplexId>& sampling_indexes_next_scale = sampling[iScale+1];
//       const double* coordinates_next_scale = fitCoordinates[iScale+1].data();


//     SimplexId nbPointNextScale = sampling[iScale+1].size();

//       // GO FOR DISTANCE BETWEEN CONSECUTIVE SCALES
// // #ifdef TTK_ENABLE_OPENMP
// // #pragma omp parallel for num_threads(1)
// // #endif // TTK_ENABLE_OPENMP
//     for(LongSimplexId I=0; I<(LongSimplexId)nbPointsAtScale; I++){

//       LongSimplexId i = sampling_indexes[I];
      
//       const double* coordsi = &coordinates[dim*i];

//       for(LongSimplexId J=0; J<(LongSimplexId)nbPointNextScale; J++){

//         LongSimplexId j = sampling_indexes_next_scale[J];

//         const double* coordsj = &coordinates_next_scale[dim*j];

//         // index in the full matrix
//         LongSimplexId id = J + I*nbPointNextScale;

//         double dist = distanceFunction(coordsi, coordsj, &mlsCoordinates[3*i], &mlsCoordinates[3*j], pointCloudBarycenter.data());

//         edges_[edgesBtwScalesStartIdx[iScale] + id] = std::make_tuple(I+pointsAtScaleStartIdx[iScale], J+pointsAtScaleStartIdx[iScale+1], dist);
//         numberOfEdges++;
//         numberEdgeWithNextScale++;
//       }
//     }

//     }
//     std::cout<<"   - "<<numberEdgeAtScale<<" edges intern found"<<std::endl;
//     std::cout<<"   - "<<numberEdgeWithNextScale<<" edges w/ nest one found"<<std::endl;
//   }

//   std::cout<<" Base Method: "<<edges_.size()<<" edges in "<<tm.getElapsedTime()<<" s."<<std::endl;
//   std::cout<<" compared to : "<<numberOfEdges<<" effective edges"<<std::endl;

//   const auto cmp = [=](const edgeTuple &e0, const edgeTuple &e1){
//     return compareEdgeTuple(e0, e1);
//   };

//   TTK_PSORT(this->threadNumber_, edges_.begin(), edges_.end(), cmp);
//   printMsg("Done, found ");
// }

// void RipsPersistence::computeEdgeVectorAtScaleND(float* coordinates, int dim, std::vector<SimplexId>& sampling_indexes, std::vector<SimplexId>& densitySampling){

//   Timer tm{};
//   printMsg("Computing exact list of edges in "+std::to_string(dim)+"D, WITH SUBSAMPLING");

//   // auto& tree = sampler.kdTree(scale_index);
//   // auto& sampling_indexes = tree.samples();
//   // SimplexId nbPoints = tree.sample_count();
//   SimplexId nbPoints = densitySampling.size();
//   // std::cout<<"scale index "<<scale_index<<" vs sampler scale count "<<sampler.scale_count()<<std::endl;
//   printMsg(std::to_string(nbPoints)+" points here");
//   printMsg("size of sampling vec : "+std::to_string(sampling_indexes.size()));
//   // if(densitySampling.size()!=nbPoints){
//   //   printErr("size != nbPoints /!\\");
//   // }

//   LongSimplexId size = (LongSimplexId)nbPoints*((LongSimplexId)nbPoints-1)/2;
//   edges_.clear();
//   edges_.resize(size);


// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(1)
// #endif // TTK_ENABLE_OPENMP
//   for(LongSimplexId I=0; I<(LongSimplexId)nbPoints-1; I++){
      
//       LongSimplexId i = densitySampling[I];

//       double ni = 0;
//       for(int k=0; k<dim; k++){
//         ni += coordinates[dim*i+k]*coordinates[dim*i+k];
//       }


//     for(LongSimplexId J=I+1; J<(LongSimplexId)nbPoints; J++){
      
//       LongSimplexId j = densitySampling[J];


//       // index in the upper triangular distance matrix
//       LongSimplexId id = J-1 + I*(2*nbPoints-I-3)/2;

//       double distCoords2 = 0; // ((xj-xi)*(xj-xi) + (yj-yi)*(yj-yi) + (zj-zi)*(zj-zi));
//       for(int k=0; k<dim; k++){
//         distCoords2 += (coordinates[dim*i+k]-coordinates[dim*j+k])*(coordinates[dim*i+k]-coordinates[dim*j+k]);
//       }
//       double dist = std::sqrt(distCoords2);

//       double nj=0;
//       for(int k=0; k<dim; k++){
//         nj += coordinates[dim*j+k]*coordinates[dim*j+k];
//       }

//       double dist2=0;
//       for(int k=0; k<dim; k++){
//         dist2 += coordinates[dim*i+k]*coordinates[dim*j+k];
//       }
//       double dist3=dist2*dist2/(ni*nj);
//       double dist4 = dist2;
//       // std::sqrt(dist3);
//       if(dist4>1){
//         dist4=1;
//       }else if(dist4<-1){
//         dist4=-1;
//       }
//       double dist5 = std::acos(dist4);

//       // std::cout<<"distances: proj "<<dist5<<" , euclidean "<<dist<<" "<<ni<<" "<<nj<<std::endl;
//       if(useProjectiveNorm_){
//         dist = dist5;
//       }
//       // if(std::isnan(dist5)){
//       //   std::cout<<"BIG NAN dist "<<dist<<" "<<ni<<" "<<nj<<std::endl;
//       //   std::cout<<"        dists "<<dist<<" "<<dist2<<" "<<dist3<<" "<<dist4<<" "<<dist5<<std::endl;
//       // }

//       if((ni<=1e-6 or nj<=1e-6) and useProjectiveNorm_){
//       }else{
//         edges_[id] = std::make_tuple(I, J, dist);
//       }
//     }
//   }

//   std::cout<<" Base Method: "<<edges_.size()<<" edges in "<<tm.getElapsedTime()<<" s."<<std::endl;

//   const auto cmp = [=](const edgeTuple &e0, const edgeTuple &e1){
//     return compareEdgeTuple(e0, e1);
//   };

//   TTK_PSORT(this->threadNumber_, edges_.begin(), edges_.end(), cmp);
//   printMsg("Done, found ");
// }


// void RipsPersistence::computeAcceleratedEdgeVectorMultiScaleND(const std::vector<std::vector<SimplexId>>& sampling,
//                                                     const std::vector<std::vector<double>>& fitCoordinates,
//                                                     const std::vector<std::vector<float>>& mls3DCoordinates,
//                                                     const std::vector<float>& pointCloudBarycenter,
//                                                     std::function<double(const double*, const double*, const float*, const float*, const float*)> distanceFunction,
//                                                     int dim,
//                                                     double maxEdgeLength){

//   Timer tm{};
//   printMsg("Computing exact list of edges in "+std::to_string(dim)+"D, WITH SUBSAMPLING, all scales mixed, below "+std::to_string(maxEdgeLength));

//   SimplexId nbScales = sampling.size();

//   edges_.clear();
//   edges_.reserve(10000000);

//   if(maxEdgeLength<0){
//     maxEdgeLength = std::numeric_limits<double>::max();
//   }

//   // LongSimplexId nbEdges = nbEdgeMonoScale * nbScales + nbEdgeBtwScales*(nbScales-1);
//   LongSimplexId nbEdges = 0;
//   std::vector<SimplexId> edgesAtScaleStartIdx(nbScales,0);
//   std::vector<SimplexId> edgesBtwScalesStartIdx(nbScales,0);
//   std::vector<SimplexId> pointsAtScaleStartIdx(nbScales,0);

//   for (int iScale=0; iScale<nbScales; iScale++){
//     SimplexId nbPointsAtScale = sampling[iScale].size();
//     LongSimplexId nbEdgeAtScale = (LongSimplexId)nbPointsAtScale*((LongSimplexId)nbPointsAtScale-1)/2;
//     LongSimplexId nbEdgeWithNextScale = iScale<nbScales-1 ? (LongSimplexId)nbPointsAtScale*((LongSimplexId)sampling[iScale+1].size()) : 0;
//     std::cout<<" scale "<<iScale<<std::endl;
//     std::cout<<"   - "<<nbPointsAtScale<<" points"<<std::endl;
//     std::cout<<"   - "<<nbEdgeAtScale<<" edges intern"<<std::endl;
//     std::cout<<"   - "<<nbEdgeWithNextScale<<" edges w/ next scale"<<std::endl;

//     edgesBtwScalesStartIdx[iScale] = edgesAtScaleStartIdx[iScale] + nbEdgeAtScale;

//     if(iScale<nbScales-1){
//       edgesAtScaleStartIdx[iScale+1] = edgesBtwScalesStartIdx[iScale] + nbEdgeWithNextScale;
//       pointsAtScaleStartIdx[iScale+1] = pointsAtScaleStartIdx[iScale] + nbPointsAtScale;
//     }

//     nbEdges += nbEdgeAtScale + nbEdgeWithNextScale;
//   }

//   // edges_.resize(nbEdges);


//     LongSimplexId numberOfEdges=0;
//   for (int iScale=0; iScale<nbScales; iScale++){

//     std::cout<<" Scale "<<iScale<< " : "<<std::endl;

//     SimplexId nbPointsAtScale = sampling[iScale].size();
//     LongSimplexId nbEdgeAtScale = (LongSimplexId)nbPointsAtScale*((LongSimplexId)nbPointsAtScale-1)/2;
//     LongSimplexId nbEdgeWithNextScale = iScale<nbScales-1 ? 0 : (LongSimplexId)nbPointsAtScale*((LongSimplexId)sampling[iScale+1].size());
//     // SimplexId nbPoints = sampling[0].size();
//     // LongSimplexId nbEdgeMonoScale = (LongSimplexId)nbPoints*((LongSimplexId)nbPoints-1)/2;
//     // LongSimplexId nbEdgeBtwScales = (LongSimplexId)nbPoints*((LongSimplexId)nbPoints);
//     // LongSimplexId nbEdges = nbEdgeMonoScale * nbScales + nbEdgeBtwScales*(nbScales-1);


//     // if(nbPoints!=sampling[iScale].size()){
//     //   printWrn("BIG ISSUE: different number of points in each scale");
//     // }

//     const std::vector<SimplexId>& sampling_indexes = sampling[iScale];
//     const double* coordinates = fitCoordinates[iScale].data();
//     const float* mlsCoordinates = mls3DCoordinates[iScale].data();
//     LongSimplexId numberEdgeAtScale=0;
//     LongSimplexId numberEdgeWithNextScale=0;

// // #ifdef TTK_ENABLE_OPENMP
// // #pragma omp parallel for num_threads(1)
// // #endif // TTK_ENABLE_OPENMP
// //
//     for(LongSimplexId I=0; I<(LongSimplexId)nbPointsAtScale-1; I++){

//       LongSimplexId i = sampling_indexes[I];

//       const double* coordsi = &coordinates[dim*i];

//       for(LongSimplexId J=I+1; J<(LongSimplexId)nbPointsAtScale; J++){

//         LongSimplexId j = sampling_indexes[J];

//         const double* coordsj = &coordinates[dim*j];

//         // index in the upper triangular distance matrix
//         LongSimplexId id = J-1 + I*(2*nbPointsAtScale-I-3)/2;

//         double dist = distanceFunction(coordsi, coordsj, &mlsCoordinates[3*i], &mlsCoordinates[3*j], pointCloudBarycenter.data());

//         if(dist <= maxEdgeLength){
//           edges_.emplace_back(I+pointsAtScaleStartIdx[iScale], J+pointsAtScaleStartIdx[iScale], dist);
//           numberEdgeAtScale++;
//           numberOfEdges++;
//         }
//       }
//     }

//     if(iScale<nbScales-1){ 

//       const std::vector<SimplexId>& sampling_indexes_next_scale = sampling[iScale+1];
//       const double* coordinates_next_scale = fitCoordinates[iScale+1].data();


//     SimplexId nbPointNextScale = sampling[iScale+1].size();

//       // GO FOR DISTANCE BETWEEN CONSECUTIVE SCALES
// // #ifdef TTK_ENABLE_OPENMP
// // #pragma omp parallel for num_threads(1)
// // #endif // TTK_ENABLE_OPENMP
//     for(LongSimplexId I=0; I<(LongSimplexId)nbPointsAtScale; I++){

//       LongSimplexId i = sampling_indexes[I];
      
//       const double* coordsi = &coordinates[dim*i];

//       for(LongSimplexId J=0; J<(LongSimplexId)nbPointNextScale; J++){

//         LongSimplexId j = sampling_indexes_next_scale[J];

//         const double* coordsj = &coordinates_next_scale[dim*j];

//         // index in the full matrix
//         LongSimplexId id = J + I*nbPointNextScale;

//         double dist = distanceFunction(coordsi, coordsj, &mlsCoordinates[3*i], &mlsCoordinates[3*j], pointCloudBarycenter.data());

//         if(dist <= maxEdgeLength){
//           edges_.emplace_back(I+pointsAtScaleStartIdx[iScale], J+pointsAtScaleStartIdx[iScale+1], dist);
//           numberOfEdges++;
//           numberEdgeWithNextScale++;
//         }
//       }
//     }

//     }
//     std::cout<<"   - "<<numberEdgeAtScale<<" edges intern found"<<std::endl;
//     std::cout<<"   - "<<numberEdgeWithNextScale<<" edges w/ nest one found"<<std::endl;
//   }

//   std::cout<<" Base Method: "<<nbEdges<<" edges in "<<tm.getElapsedTime()<<" s."<<std::endl;
//   std::cout<<" compared to : "<<numberOfEdges<<" effective edges"<<std::endl;

//   const auto cmp = [=](const edgeTuple &e0, const edgeTuple &e1){
//     return compareEdgeTuple(e0, e1);
//   };

//   TTK_PSORT(this->threadNumber_, edges_.begin(), edges_.end(), cmp);
//   printMsg("Done, found ");
// }


// template<class PointND>
// void RipsPersistence::computeApproximatedEdgeVectorMultiScaleND(std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
//                                                     std::vector<std::vector<SimplexId>>& filteredToLocalId,
//                                                     std::vector<std::vector<double>>& fitCoordinates,
//                                                     std::vector<std::vector<float>>& mls3DCoordinates,
//                                                     std::vector<float>& pointCloudBarycenter,
//                                                     std::function<double(const double*, const double*, const float*, const float*, const float*)> distanceFunction,
//                                                     int dim,
//                                                     double maxEdgeLength){

//   Timer tm{};
//   printMsg("Computing exact list of edges in "+std::to_string(dim)+"D, WITH SUBSAMPLING, all scales mixed, below "+std::to_string(maxEdgeLength));

//   SimplexId nbScales = epsilonIdToFilteredId.size();

//   edges_.clear();
//   edges_.reserve(10000000);

//   if(maxEdgeLength<0){
//     maxEdgeLength = std::numeric_limits<double>::max();
//   }

//   long long int numberOfEdges = 0;


//   std::vector<SimplexId> pointsAtScaleStartIdx(nbScales,0);
//   std::vector<ScaleSpaceSampling<PointND::Scalar, PointND>> sampler5D(nbScales);
//   double diag_aabb_5D = 2;


//   SimplexId totalNbPoints=0;
//   for (int iScale=0; iScale<nbScales; iScale++){
//     totalNbPoints+=epsilonIdToFilteredId[iScale].size();
//   }
//   std::cout<<" Total nb of points : "<<totalNbPoints<<std::endl;
// // #ifdef TTK_ENABLE_OPENMP
// // #pragma omp parallel for num_threads(threadNumber_)
// // #endif // TTK_ENABLE_OPENMP
//   for (int iScale=0; iScale<nbScales-1; iScale++){
//     auto& sampler = sampler5D[iScale];
//     sampler.setDebugLevel(debugLevel_);
//     sampler.setSamplingFactor(0.1);
//     sampler.sampleScaleLog(1e-5, 2, 1);
//     sampler.poissonDiskSamplingBucketIndexing(fitCoordinates, epsilonIdToFilteredId, filteredToLocalId,  iScale, iScale+1);

//     SimplexId nbPointsAtScale = epsilonIdToFilteredId[iScale].size();
//     // std::cout<<" scale "<<iScale<<" : "<<nbPointsAtScale<<" points"<<std::endl;
//     pointsAtScaleStartIdx[iScale+1] = pointsAtScaleStartIdx[iScale] + nbPointsAtScale;

//   }


//   // { // last level
//   //   auto& sampler = sampler5D[nbScales-1];
//   //   sampler.setDebugLevel(debugLevel_);
//   //   sampler.setSamplingFactor(0.1);
//   //   sampler.sampleScaleLog(0.001*diag_aabb_5D, 1*diag_aabb_5D, 5);
//   //   sampler.poissonDiskSamplingBucketIndexing(fitCoordinates, sampling, nbScales-1);
//   // }

//   for (int iScale=0; iScale<nbScales-1; iScale++){
//     auto& sampler = sampler5D[iScale];

//     int nbLevels = sampler.scale_count();
//     bool last_scale = (iScale+1 == nbScales-1);


//     std::cout<<"\nScale "<<iScale<<std::endl;
//     for(int level=nbLevels-1; level>=nbLevels-1; level--){
//       // double distMax = level == nbLevels-1 ? std::numeric_limits<double>::infinity() : sampler.scale(level)*sampler.samplingFactor();
//       double distMax = std::numeric_limits<double>::max();
//       distMax = std::min(distMax, (double)maxEdgeLength);
//       long long unsigned int addedEdges=0;
//       double max_added_dist=0;
//       std::cout<<" \nLevel "<<level<<", distmax = "<<distMax<<std::endl;
//       auto& tree = sampler.kdTree(level);
//       LongSimplexId nbPointsAtLevel = tree.samples().size();

//       std::cout<<" tree samples "<<nbPointsAtLevel<< " (potentially "<<nbPointsAtLevel*(nbPointsAtLevel-1)<<" edges )"<<std::endl;

//       // std::cout<<" tree points: ";
//       // for(SimplexId isamp=0; isamp<tree.sample_count(); isamp++){
//       //   std::cout<<" "<<tree.samples()[isamp]<<" ("<<tree.pointDataFromSample(isamp).pos().transpose()<<")";
//       // }
//       // std::cout<<std::endl;

//       for(LongSimplexId id=0; id<nbPointsAtLevel; id++){
//         // std::cout<<" id "<<id<<std::endl;
//         if(id>=tree.samples().size()){
//           std::cout<<" error 1"<<std::endl;
//         }
//         SimplexId index_biscale_i = tree.samples()[id];
//         SimplexId index_monoscale_i, I; 
//         int idScale=-1;
//         const double* coordinates_i{}; 
//         const float* mlsCoordinates_i;
//         if(index_biscale_i < epsilonIdToFilteredId[iScale].size()){ // point in iScale
//           index_monoscale_i = index_biscale_i;
//           idScale=iScale;
//           // std::cout<<"there 0 id: "<<id<<" "<<" biscale : "<<index_biscale_i<<" monoscale :"<<index_monoscale_i<<std::endl;
//           // std::cout<<"sampling size "<<sampling[iScale].size()<<std::endl;
//         } else{ // point in iScale+1
//           index_monoscale_i = index_biscale_i - epsilonIdToFilteredId[iScale].size(); // verifier
//           idScale = iScale+1;
//         }
//         if(index_monoscale_i>=epsilonIdToFilteredId[idScale].size()){
//           std::cout<<" error 2"<<std::endl;
//         }
//         if(epsilonIdToFilteredId[idScale][index_monoscale_i]>=filteredToLocalId[idScale].size()){
//           std::cout<<" error 3"<<std::endl;
//         }
//         I = filteredToLocalId[idScale][epsilonIdToFilteredId[idScale][index_monoscale_i]]; //local id
                                                                                           
//         if(idScale>=fitCoordinates.size()){
//           std::cout<<" error 4"<<std::endl;
//         }
//         coordinates_i = fitCoordinates[idScale].data();

//         if(idScale>=mls3DCoordinates.size()){
//           std::cout<<" error 5"<<std::endl;
//         }
//         mlsCoordinates_i = mls3DCoordinates[idScale].data();

//         if(dim*I+4>=fitCoordinates[idScale].size()){
//           std::cout<<" error 6"<<std::endl;
//         }
//         const double* coordsi = &coordinates_i[dim*I];
//         if(idScale>=pointsAtScaleStartIdx.size()){
//           std::cout<<" error 7"<<std::endl;
//         }
//         const SimplexId index_multiscale_i = index_monoscale_i+pointsAtScaleStartIdx[idScale];
//         if(index_biscale_i>=tree.point_count()){
//           std::cout<<" error 8"<<std::endl;
//         }
 
//         auto tree_query =tree.range_neighbors(index_biscale_i,distMax);

//         SimplexId nMatches=0;
//         for(LongSimplexId jd : tree_query){
//           nMatches++;

//           LongSimplexId index_biscale_j = jd;
//           SimplexId index_monoscale_j, J; 
//           const double* coordinates_j;
//           const float* mlsCoordinates_j;
//           int jdScale=-1;

//           if(index_biscale_j < epsilonIdToFilteredId[iScale].size()){ // point in iScale
//             index_monoscale_j = index_biscale_j;
//             jdScale=iScale;
//           } else{ // point in iScale+1
//             index_monoscale_j = index_biscale_j - epsilonIdToFilteredId[iScale].size(); // verifier
//             jdScale=iScale+1;
//           }

//         if(index_monoscale_j>=epsilonIdToFilteredId[jdScale].size()){
//           std::cout<<" error 12"<<std::endl;
//         }
//         if(epsilonIdToFilteredId[jdScale][index_monoscale_j]>=filteredToLocalId[jdScale].size()){
//           std::cout<<" error 13"<<std::endl;
//         }

//           J = filteredToLocalId[jdScale][epsilonIdToFilteredId[jdScale][index_monoscale_j]];

//         if(jdScale>=fitCoordinates.size()){
//           std::cout<<" error 14"<<std::endl;
//         }
//         if(jdScale>=mls3DCoordinates.size()){
//           std::cout<<" error 15"<<std::endl;
//         }

//           coordinates_j = fitCoordinates[jdScale].data();
//           mlsCoordinates_j = mls3DCoordinates[jdScale].data();
//           // std::cout<<" jd here 1 "<<std::endl;
//         if(dim*J+4>=fitCoordinates[jdScale].size()){
//           std::cout<<" error 16"<<std::endl;
//         }
//           const double* coordsj = &coordinates_j[dim*J];
//         if(jdScale>=pointsAtScaleStartIdx.size()){
//           std::cout<<" error 17"<<std::endl;
//         }
//           const SimplexId index_multiscale_j = index_monoscale_j+pointsAtScaleStartIdx[jdScale];

//           bool both_second_scale = (idScale == iScale+1 and jdScale == iScale+1);

//           if(index_multiscale_i < index_multiscale_j and (!both_second_scale or last_scale)){

//         if(3*I+2>=mls3DCoordinates[idScale].size()){
//           std::cout<<" error 20"<<std::endl;
//         }
//         if(3*J+2>=mls3DCoordinates[jdScale].size()){
//           std::cout<<" error 21"<<std::endl;
//         }
//             double dist = distanceFunction(coordsi, coordsj, &mlsCoordinates_i[3*I], &mlsCoordinates_j[3*J], pointCloudBarycenter.data());

//             if(dist <= maxEdgeLength ){

//               edges_.emplace_back(index_multiscale_i, index_multiscale_j, dist);
//               // std::cout<<"                    added edge "<<index_monoscale_i+pointsAtScaleStartIdx[idScale]<<" - "<<index_monoscale_j+pointsAtScaleStartIdx[jdScale]<<" :"<<dist<<std::endl;
//               numberOfEdges++;
//               addedEdges++;
//               // if(addedEdges%10000000 == 0){
//               //   std::cout<<" already "<<addedEdges<<" edges"<<std::endl;
//               // }
//               // if(jdScale==idScale and dist<1e-4){
//               //   std::cout<<" THIS SHOULD NOT HAPPEN: d="<<dist<<" btw "<<I<<"("<<idScale<<") and "<<J<<"("<<jdScale<<")"<<std::endl;
//               // }
//               if(dist>max_added_dist){
//                 max_added_dist=dist;
//               }
//             }
//           }

//         }
//         // std::cout<<" end tree query"<<std::endl;
//         // std::cout<<" found "<<nMatches<<" matches w/ dist "<<distMax<<std::endl;
//         // std::cout<<" done jd "<<std::endl;
//       }
//       std::cout<<" added "<<addedEdges<<" edges with max dist "<<max_added_dist<<std::endl;
//     }

//     // now add edges of the lowest level
//     // if(nbLevels>1){
//     //   int level = 0;
//     for(int level=nbLevels-2; level>=0; level--){
//       double distMax = level < nbLevels-1 ? sampler.scale(level)*sampler.samplingFactor()*1.1 : std::numeric_limits<double>::max();
//       SimplexId addedEdges=0;
//       double max_added_dist=0;
//       double min_added_dist=std::numeric_limits<double>::infinity();
//       std::cout<<" Level "<<level<<", distmax = "<<distMax<<std::endl;
//       auto& tree1 = sampler.kdTree(level+1);
//       auto& tree0 = sampler.kdTree(level);
//       LongSimplexId nbPointsAtLevel1 = tree1.samples().size();
//       LongSimplexId nbPointsAtLevel0 = tree0.samples().size();
//       // std::map<SimplexId,int> isMatched{};
//       // std::cout<<" points: "<<nbPointsAtLevel1<<" and "<<nbPointsAtLevel0<<std::endl;
//       // std::cout<<" tree1: ";
//       // for(LongSimplexId id=0; id<nbPointsAtLevel1; id++){
//       //   SimplexId index_biscale_i = tree1.samples()[id];
//       //   std::cout<<" "<<index_biscale_i;
//       // }
//       // std::cout<<std::endl;
//       // std::cout<<" tree0: ";
//       // for(LongSimplexId id=0; id<nbPointsAtLevel0; id++){
//       //   SimplexId index_biscale_i = tree0.samples()[id];
//       //   std::cout<<" "<<index_biscale_i;
//       //   isMatched[index_biscale_i]=0;
//       // }
//       // std::cout<<std::endl;


//       for(LongSimplexId id=0; id<nbPointsAtLevel1; id++){
//         // std::cout<<" id "<<id<<std::endl;
//         // std::cout<<" tree samples "<<tree.samples().size()<<std::endl;
//         SimplexId index_biscale_i = tree1.samples()[id];
//         SimplexId index_monoscale_i, I; 

//         // isMatched[index_biscale_i]+=1;

//         int idScale=-1;
//         const double* coordinates_i{}; 
//         const float* mlsCoordinates_i;
//         if(index_biscale_i < epsilonIdToFilteredId[iScale].size()){ // point in iScale
//           index_monoscale_i = index_biscale_i;
//           idScale=iScale;
//           // std::cout<<"there 0 id: "<<id<<" "<<" biscale : "<<index_biscale_i<<" monoscale :"<<index_monoscale_i<<std::endl;
//           // std::cout<<"sampling size "<<sampling[iScale].size()<<std::endl;
//         } else{ // point in iScale+1
//           index_monoscale_i = index_biscale_i - epsilonIdToFilteredId[iScale].size(); // verifier
//           idScale = iScale+1;
//         }
//         I = filteredToLocalId[idScale][epsilonIdToFilteredId[idScale][index_monoscale_i]];


//         // if(index_monoscale_i+pointsAtScaleStartIdx[idScale]==176765){
//         //   std::cout<<" found 176765 in tree 1 with indices "<<index_biscale_i<<" "<<index_monoscale_i<<std::endl;
//         // }

//         // if(index_biscale_i==0 or index_biscale_i==11091){
//         //   std::cout<<"  looking at point "<<index_biscale_i<<std::endl;
//         // }
//         coordinates_i = fitCoordinates[idScale].data();
//         mlsCoordinates_i = mls3DCoordinates[idScale].data();
//         const double* coordsi = &coordinates_i[dim*I];


//         // query in the lowest level

//         // std::cout<<" query for index "<<index_biscale_i<<" : "<<std::endl;
//         auto tree_query =tree0.range_neighbors(index_biscale_i,distMax);
//         SimplexId nMatches=0;
//         for(LongSimplexId jd : tree_query){
//           nMatches++;
//           // std::cout<<"jd "<<jd<<std::endl;
//           // std::cout<<" tree samples "<<tree.samples().size()<<std::endl;
//           LongSimplexId index_biscale_j = jd;

//           // isMatched[index_biscale_j]+=1;

//           // std::cout<<" "<<index_biscale_j;
//           SimplexId index_monoscale_j, J; 
//           const double* coordinates_j;
//           const float* mlsCoordinates_j;
//           int jdScale=-1;

//           // std::cout<<" jd here 0 "<<std::endl;
//           if(index_biscale_j < epsilonIdToFilteredId[iScale].size()){ // point in iScale
//           // std::cout<<" there 0 0 "<<std::endl;
//             index_monoscale_j = index_biscale_j;
//             jdScale=iScale;
//             // std::cout<<" jd: "<<jd<<" "<<" biscale : "<<index_biscale_j<<" monoscale :"<<index_monoscale_j<<std::endl;
//             // std::cout<<"sampling size "<<sampling[iScale].size()<<std::endl;
//           // std::cout<<" here 0 0 "<<std::endl;
//           } else{ // point in iScale+1
//           // std::cout<<" there 0 1 "<<std::endl;
//             index_monoscale_j = index_biscale_j - epsilonIdToFilteredId[iScale].size(); // verifier
//             jdScale=iScale+1;
//             // std::cout<<" jd: "<<jd<<" "<<" biscale : "<<index_biscale_j<<" monoscale :"<<index_monoscale_j<<std::endl;
//             // std::cout<<"sampling size "<<sampling[iScale+1].size()<<std::endl;
//           // std::cout<<" here 0 1 "<<std::endl;
//           }
//           J = filteredToLocalId[idScale][epsilonIdToFilteredId[jdScale][index_monoscale_j]];
//             // if(index_biscale_j==11091 or index_biscale_j==0 or index_biscale_i==0 or index_biscale_i==11091){
//             // std::cout<<"                   original point "<<index_biscale_i<<std::endl;
//             // std::cout<<"                   neighbor point "<<index_biscale_j<<std::endl;
//             // }
//             // if(J+pointsAtScaleStartIdx[jdScale]==11091){
//           // std::cout<<"                   neighbor point "<<J+pointsAtScaleStartIdx[jdScale]<<std::endl;
//             // }
//             // if(J+pointsAtScaleStartIdx[jdScale]==0){
//           // std::cout<<"                   neighbor point "<<J+pointsAtScaleStartIdx[jdScale]<<std::endl;
//             // }
//         // if(I+pointsAtScaleStartIdx[idScale]==11091){
//           // std::cout<<"                   neighbor point "<<J+pointsAtScaleStartIdx[jdScale]<<std::endl;
//         // }
//         // if(I+pointsAtScaleStartIdx[idScale]==0){
//           // std::cout<<"                   neighbor point "<<J+pointsAtScaleStartIdx[jdScale]<<std::endl;
//         // }
//         // if(index_monoscale_j+pointsAtScaleStartIdx[jdScale]==176765){
//         //   std::cout<<" found 176765 in tree 0 with indices "<<index_biscale_j<<" "<<index_monoscale_j<<std::endl;
//         // }
//           coordinates_j = fitCoordinates[jdScale].data();
//           mlsCoordinates_j = mls3DCoordinates[jdScale].data();
//           // std::cout<<" jd here 1 "<<std::endl;
//           const double* coordsj = &coordinates_j[dim*J];


//           // std::cout<<" jd here 2"<<std::endl;
//           double dist = distanceFunction(coordsi, coordsj, &mlsCoordinates_i[3*I], &mlsCoordinates_j[3*J], pointCloudBarycenter.data());

//           if(dist <= maxEdgeLength){
//             // if(index_biscale_j==11091 or index_biscale_j==0 or index_biscale_i==0 or index_biscale_i==11091){
//             //   std::cout<<"                      edge "<<index_monoscale_i+pointsAtScaleStartIdx[idScale]<<" - "<<index_monoscale_j+pointsAtScaleStartIdx[jdScale]<<" :"<<dist<<std::endl;
//             // }
//             // if(I+pointsAtScaleStartIdx[idScale] ==0 or J+pointsAtScaleStartIdx[jdScale]==0){
//             // std::cout<<"                      edge "<<I+pointsAtScaleStartIdx[idScale]<<" - "<<J+pointsAtScaleStartIdx[jdScale]<<" :"<<dist<<std::endl;
//             // }
//             // if(I+pointsAtScaleStartIdx[idScale] ==11091 or J+pointsAtScaleStartIdx[jdScale]==11091){
//             // std::cout<<"                      edge "<<I+pointsAtScaleStartIdx[idScale]<<" - "<<J+pointsAtScaleStartIdx[jdScale]<<" :"<<dist<<std::endl;
//             // }
//             edges_.emplace_back(index_monoscale_i+pointsAtScaleStartIdx[idScale], index_monoscale_j+pointsAtScaleStartIdx[jdScale], dist);
//             numberOfEdges++;
//             addedEdges++;
//             if(dist>max_added_dist){
//               max_added_dist=dist;
//             }
//             if(dist<min_added_dist){
//               min_added_dist=dist;
//             }
//           }else{
//             // std::cout<<"edge "<<dist<<" discarded"<<std::endl;
//           }

//           // std::cout<<" done jd "<<std::endl;
//         }
//         // std::cout<<std::endl;
//         // std::cout<<" found "<<nMatches<<" matches w/ dist "<<distMax<<std::endl;
//       }
//       // for(LongSimplexId id=0; id<nbPointsAtLevel0; id++){
//       //   SimplexId index_biscale_i = tree0.samples()[id];
//       //   if(isMatched[index_biscale_i]==0){
//       //     std::cout<<"  ERROR: POINT "<<index_biscale_i<<" was not matched to any point in tree1"<<std::endl;
//       //   }
//       // }
//       std::cout<<" added "<<addedEdges<<" edges with max dist "<<max_added_dist<<" and min dist "<<min_added_dist<<std::endl;
//     }
//   }



//   std::cout<<" Found: "<<numberOfEdges<<" edges in "<<tm.getElapsedTime()<<" s."<<std::endl;
//   // std::cout<<" compared to : "<<numberOfEdges<<" effective edges"<<std::endl;

//   const auto cmp = [=](const edgeTuple &e0, const edgeTuple &e1){
//     return compareEdgeTuple(e0, e1);
//   };

//   TTK_PSORT(this->threadNumber_, edges_.begin(), edges_.end(), cmp);
//   printMsg("Done, found ");
// }

// void RipsPersistence::computeEpsilonSamplingEdgeVectorMultiScaleND(std::vector<std::vector<SimplexId>>& sampling,
//                                                     std::vector<std::vector<float>>& fitCoordinates,
//                                                     std::vector<std::vector<float>>& mls3DCoordinates,
//                                                     std::vector<float>& pointCloudBarycenter,
//                                                     std::function<double(const float*, const float*, const float*, const float*, const float*)> distanceFunction,
//                                                     std::vector<EpsilonSampling<pointCloud::MyPoint5D::Scalar, pointCloud::MyPoint5D>>& epsilonSamplers,
//                                                     int dim,
//                                                     double maxEdgeLength){

//   Timer tm{};
//   printMsg("Computing exact list of edges in "+std::to_string(dim)+"D, WITH Epsilon SUBSAMPLING, all scales mixed, below "+std::to_string(maxEdgeLength));

//   SimplexId nbScales = sampling.size();

//   edges_.clear();
//   edges_.reserve(10000000);

//   if(maxEdgeLength<0){
//     maxEdgeLength = std::numeric_limits<double>::max();
//   }

//   SimplexId numberOfEdges = 0;


//   std::vector<SimplexId> pointsAtScaleStartIdx(nbScales,0);
//   // std::vector<EpsilonSampling<pointCloud::MyPoint5D::Scalar, pointCloud::MyPoint5D>> epsilonSamplers(nbScales);
//   epsilonSamplers.resize(nbScales);
//   std::vector<ScaleSpaceSampling<pointCloud::MyPoint5D::Scalar, pointCloud::MyPoint5D>> sampler5D(nbScales);
//   double diag_aabb_5D = 2;

//   std::vector<std::vector<SimplexId>> epsilon_sampling(nbScales);

//   SimplexId totalNbPoints=0;
//   SimplexId epsSamplingNbPoints=0;
//   for (int iScale=0; iScale<nbScales; iScale++){
//     totalNbPoints+=sampling[iScale].size();
//     SimplexId nbPointsAtScale = sampling[iScale].size();

//     auto& epsSampler = epsilonSamplers[iScale];
//     epsSampler.setEpsilon(1e-5);
//     epsSampler.buildSampling(nbPointsAtScale, fitCoordinates[iScale].data(), sampling[iScale]);
//     epsilon_sampling[iScale] = epsSampler.samples();  // copy

//     epsSamplingNbPoints+=epsilon_sampling[iScale].size();
//   }

//   std::cout<<" Total nb of points : "<<totalNbPoints<<std::endl;
//   std::cout<<" Eps sampling nb of points : "<<epsSamplingNbPoints<<std::endl;

// // #ifdef TTK_ENABLE_OPENMP
// // #pragma omp parallel for num_threads(threadNumber_)
// // #endif // TTK_ENABLE_OPENMP
//   for (int iScale=0; iScale<nbScales; iScale++){
//     auto& sampler = sampler5D[iScale];
//     sampler.setDebugLevel(debugLevel_);
//     sampler.setSamplingFactor(0.1);
//     sampler.sampleScaleLog(0.0001*diag_aabb_5D, 1*diag_aabb_5D, 1);
//     sampler.poissonDiskSamplingBucketIndexing(fitCoordinates, epsilon_sampling, iScale, iScale+1);

//     SimplexId nbPointsAtScale = epsilon_sampling[iScale].size();
//     pointsAtScaleStartIdx[iScale+1] = pointsAtScaleStartIdx[iScale] + nbPointsAtScale;

//   }


//   // { // last level
//   //   auto& sampler = sampler5D[nbScales-1];
//   //   sampler.setDebugLevel(debugLevel_);
//   //   sampler.setSamplingFactor(0.1);
//   //   sampler.sampleScaleLog(0.001*diag_aabb_5D, 1*diag_aabb_5D, 5);
//   //   sampler.poissonDiskSamplingBucketIndexing(fitCoordinates, sampling, nbScales-1);
//   // }

//   for (int iScale=0; iScale<nbScales-1; iScale++){
//     auto& sampler = sampler5D[iScale];

//     int nbLevels = sampler.scale_count();

//     for(int level=nbLevels-1; level>=nbLevels-1; level--){
//       // double distMax = level == nbLevels-1 ? std::numeric_limits<double>::infinity() : sampler.scale(level)*sampler.samplingFactor();
//       float distMax = std::numeric_limits<float>::max();
//       distMax = std::min(distMax, (float)maxEdgeLength);
//       SimplexId addedEdges=0;
//       double max_added_dist=0;
//       std::cout<<" Level "<<level<<", distmax = "<<distMax<<std::endl;
//       auto& tree = sampler.kdTree(level);
//       LongSimplexId nbPointsAtLevel = tree.samples().size();

//       std::cout<<" tree samples "<<nbPointsAtLevel<< " (potentially "<<nbPointsAtLevel*(nbPointsAtLevel-1)<<" edges )"<<std::endl;
//       for(LongSimplexId id=0; id<nbPointsAtLevel; id++){
//         // std::cout<<" id "<<id<<std::endl;
//         SimplexId index_biscale_i = tree.samples()[id];
//         SimplexId index_monoscale_i, I; 
//         int idScale=-1;
//         const float* coordinates_i{}; 
//         const float* mlsCoordinates_i;
//         if(index_biscale_i < epsilon_sampling[iScale].size()){ // point in iScale
//           index_monoscale_i = index_biscale_i;
//           idScale=iScale;
//           // std::cout<<"there 0 id: "<<id<<" "<<" biscale : "<<index_biscale_i<<" monoscale :"<<index_monoscale_i<<std::endl;
//           // std::cout<<"sampling size "<<sampling[iScale].size()<<std::endl;
//         } else{ // point in iScale+1
//           index_monoscale_i = index_biscale_i - epsilon_sampling[iScale].size(); // verifier
//           idScale = iScale+1;
//         }
//         I = epsilon_sampling[idScale][index_monoscale_i];
//         coordinates_i = fitCoordinates[idScale].data();
//         mlsCoordinates_i = mls3DCoordinates[idScale].data();
//         const float* coordsi = &coordinates_i[dim*I];


//         // if(index_monoscale_i+pointsAtScaleStartIdx[idScale]==176765){
//         //   std::cout<<" found 176765 in last level with indices "<<index_biscale_i<<" "<<index_monoscale_i<<std::endl;
//         // }

//         // std::cout<<"go for jd"<<std::endl;
//         // if(index_biscale_i==0 or index_biscale_i==11091){
//           // std::cout<<"  looking at point "<<index_biscale_i<<std::endl;
//         // }
//         auto tree_query = tree.range_neighbors(index_biscale_i,distMax);

//         SimplexId nMatches=0;
//         for(LongSimplexId jd : tree_query){
//           nMatches++;
//           // std::cout<<" "<<jd; //<<std::endl;
//           // std::cout<<" tree samples "<<tree.samples().size()<<std::endl;

//           LongSimplexId index_biscale_j = jd;
//           SimplexId index_monoscale_j, J; 
//           const float* coordinates_j;
//           const float* mlsCoordinates_j;
//           int jdScale=-1;

//           // std::cout<<" jd here 0 "<<std::endl;
//           if(index_biscale_j < epsilon_sampling[iScale].size()){ // point in iScale
//           // std::cout<<" there 0 0 "<<std::endl;
//             index_monoscale_j = index_biscale_j;
//             jdScale=iScale;
//             // std::cout<<" jd: "<<jd<<" "<<" biscale : "<<index_biscale_j<<" monoscale :"<<index_monoscale_j<<std::endl;
//             // std::cout<<"sampling size "<<sampling[iScale].size()<<std::endl;
//           // std::cout<<" here 0 0 "<<std::endl;
//           } else{ // point in iScale+1
//           // std::cout<<" there 0 1 "<<std::endl;
//             index_monoscale_j = index_biscale_j - epsilon_sampling[iScale].size(); // verifier
//             jdScale=iScale+1;
//             // std::cout<<"sampling size "<<sampling[iScale+1].size()<<std::endl;
//           // std::cout<<" here 0 1 "<<std::endl;
//           }


//         // if(index_monoscale_j+pointsAtScaleStartIdx[jdScale]==176765){
//         //   std::cout<<" found 176765 in last level query with indices "<<index_biscale_j<<" "<<index_monoscale_j<<std::endl;
//         // }


//           // std::cout<<" jd: "<<jd<<" "<<" biscale : "<<index_biscale_j<<" monoscale :"<<index_monoscale_j<<std::endl;
//           J = epsilon_sampling[jdScale][index_monoscale_j];
//           coordinates_j = fitCoordinates[jdScale].data();
//           mlsCoordinates_j = mls3DCoordinates[jdScale].data();
//           // std::cout<<" jd here 1 "<<std::endl;
//           const float* coordsj = &coordinates_j[dim*J];

//             // if(index_biscale_j==11091 or index_biscale_j==0 or index_biscale_i==0 or index_biscale_i==11091){
//             // std::cout<<"                   original point "<<index_biscale_i<<std::endl;
//             // std::cout<<"                   neighbor point "<<index_biscale_j<<std::endl;
//             // }

//           // std::cout<<" jd here 2"<<std::endl;
//           double dist = distanceFunction(coordsi, coordsj, &mlsCoordinates_i[3*I], &mlsCoordinates_j[3*J], pointCloudBarycenter.data());

//           // std::cout<<"                     potential edge "<<index_monoscale_i+pointsAtScaleStartIdx[idScale]<<" - "<<index_monoscale_j+pointsAtScaleStartIdx[jdScale]<<" :"<<dist<<" and maxedgelength "<<maxEdgeLength<<std::endl;
//           if(dist <= maxEdgeLength){
//             // if(index_biscale_j==11091 or index_biscale_j==0 or index_biscale_i==0 or index_biscale_i==11091){
//             // }
//             edges_.emplace_back(index_monoscale_i+pointsAtScaleStartIdx[idScale], index_monoscale_j+pointsAtScaleStartIdx[jdScale], dist);
//           // std::cout<<"                    added edge "<<index_monoscale_i+pointsAtScaleStartIdx[idScale]<<" - "<<index_monoscale_j+pointsAtScaleStartIdx[jdScale]<<" :"<<dist<<std::endl;
//             numberOfEdges++;
//             addedEdges++;
//             if(dist>max_added_dist){
//               max_added_dist=dist;
//             }
//           }

//         }
//         // std::cout<<" found "<<nMatches<<" matches w/ dist "<<distMax<<std::endl;
//         // std::cout<<" done jd "<<std::endl;
//       }
//       std::cout<<" added "<<addedEdges<<" edges with max dist "<<max_added_dist<<std::endl;
//     }

//     // now add edges of the lowest level
//     // if(nbLevels>1){
//     //   int level = 0;
//     for(int level=nbLevels-2; level>=0; level--){
//       float distMax = level < nbLevels-1 ? sampler.scale(level)*sampler.samplingFactor()*1.1 : std::numeric_limits<float>::max();
//       SimplexId addedEdges=0;
//       double max_added_dist=0;
//       double min_added_dist=std::numeric_limits<double>::infinity();
//       std::cout<<" Level "<<level<<", distmax = "<<distMax<<std::endl;
//       auto& tree1 = sampler.kdTree(level+1);
//       auto& tree0 = sampler.kdTree(level);
//       LongSimplexId nbPointsAtLevel1 = tree1.samples().size();
//       LongSimplexId nbPointsAtLevel0 = tree0.samples().size();
//       // std::map<SimplexId,int> isMatched{};
//       // std::cout<<" points: "<<nbPointsAtLevel1<<" and "<<nbPointsAtLevel0<<std::endl;
//       // std::cout<<" tree1: ";
//       // for(LongSimplexId id=0; id<nbPointsAtLevel1; id++){
//       //   SimplexId index_biscale_i = tree1.samples()[id];
//       //   std::cout<<" "<<index_biscale_i;
//       // }
//       // std::cout<<std::endl;
//       // std::cout<<" tree0: ";
//       // for(LongSimplexId id=0; id<nbPointsAtLevel0; id++){
//       //   SimplexId index_biscale_i = tree0.samples()[id];
//       //   std::cout<<" "<<index_biscale_i;
//       //   isMatched[index_biscale_i]=0;
//       // }
//       // std::cout<<std::endl;


//       for(LongSimplexId id=0; id<nbPointsAtLevel1; id++){
//         // std::cout<<" id "<<id<<std::endl;
//         // std::cout<<" tree samples "<<tree.samples().size()<<std::endl;
//         SimplexId index_biscale_i = tree1.samples()[id];
//         SimplexId index_monoscale_i, I; 

//         // isMatched[index_biscale_i]+=1;

//         int idScale=-1;
//         const float* coordinates_i{}; 
//         const float* mlsCoordinates_i;
//         if(index_biscale_i < sampling[iScale].size()){ // point in iScale
//           index_monoscale_i = index_biscale_i;
//           idScale=iScale;
//           // std::cout<<"there 0 id: "<<id<<" "<<" biscale : "<<index_biscale_i<<" monoscale :"<<index_monoscale_i<<std::endl;
//           // std::cout<<"sampling size "<<sampling[iScale].size()<<std::endl;
//         } else{ // point in iScale+1
//           index_monoscale_i = index_biscale_i - sampling[iScale].size(); // verifier
//           idScale = iScale+1;
//         }
//         I = sampling[idScale][index_monoscale_i];


//         // if(index_monoscale_i+pointsAtScaleStartIdx[idScale]==176765){
//         //   std::cout<<" found 176765 in tree 1 with indices "<<index_biscale_i<<" "<<index_monoscale_i<<std::endl;
//         // }

//         // if(index_biscale_i==0 or index_biscale_i==11091){
//         //   std::cout<<"  looking at point "<<index_biscale_i<<std::endl;
//         // }
//         coordinates_i = fitCoordinates[idScale].data();
//         mlsCoordinates_i = mls3DCoordinates[idScale].data();
//         const float* coordsi = &coordinates_i[dim*I];


//         // query in the lowest level

//         // std::cout<<" query for index "<<index_biscale_i<<" : "<<std::endl;
//         auto tree_query =tree0.range_neighbors(index_biscale_i,distMax);
//         SimplexId nMatches=0;
//         for(LongSimplexId jd : tree_query){
//           nMatches++;
//           // std::cout<<"jd "<<jd<<std::endl;
//           // std::cout<<" tree samples "<<tree.samples().size()<<std::endl;
//           LongSimplexId index_biscale_j = jd;

//           // isMatched[index_biscale_j]+=1;

//           // std::cout<<" "<<index_biscale_j;
//           SimplexId index_monoscale_j, J; 
//           const float* coordinates_j;
//           const float* mlsCoordinates_j;
//           int jdScale=-1;

//           // std::cout<<" jd here 0 "<<std::endl;
//           if(index_biscale_j < sampling[iScale].size()){ // point in iScale
//           // std::cout<<" there 0 0 "<<std::endl;
//             index_monoscale_j = index_biscale_j;
//             jdScale=iScale;
//             // std::cout<<" jd: "<<jd<<" "<<" biscale : "<<index_biscale_j<<" monoscale :"<<index_monoscale_j<<std::endl;
//             // std::cout<<"sampling size "<<sampling[iScale].size()<<std::endl;
//           // std::cout<<" here 0 0 "<<std::endl;
//           } else{ // point in iScale+1
//           // std::cout<<" there 0 1 "<<std::endl;
//             index_monoscale_j = index_biscale_j - sampling[iScale].size(); // verifier
//             jdScale=iScale+1;
//             // std::cout<<" jd: "<<jd<<" "<<" biscale : "<<index_biscale_j<<" monoscale :"<<index_monoscale_j<<std::endl;
//             // std::cout<<"sampling size "<<sampling[iScale+1].size()<<std::endl;
//           // std::cout<<" here 0 1 "<<std::endl;
//           }
//           J = sampling[jdScale][index_monoscale_j];
//             // if(index_biscale_j==11091 or index_biscale_j==0 or index_biscale_i==0 or index_biscale_i==11091){
//             // std::cout<<"                   original point "<<index_biscale_i<<std::endl;
//             // std::cout<<"                   neighbor point "<<index_biscale_j<<std::endl;
//             // }
//             // if(J+pointsAtScaleStartIdx[jdScale]==11091){
//           // std::cout<<"                   neighbor point "<<J+pointsAtScaleStartIdx[jdScale]<<std::endl;
//             // }
//             // if(J+pointsAtScaleStartIdx[jdScale]==0){
//           // std::cout<<"                   neighbor point "<<J+pointsAtScaleStartIdx[jdScale]<<std::endl;
//             // }
//         // if(I+pointsAtScaleStartIdx[idScale]==11091){
//           // std::cout<<"                   neighbor point "<<J+pointsAtScaleStartIdx[jdScale]<<std::endl;
//         // }
//         // if(I+pointsAtScaleStartIdx[idScale]==0){
//           // std::cout<<"                   neighbor point "<<J+pointsAtScaleStartIdx[jdScale]<<std::endl;
//         // }
//         // if(index_monoscale_j+pointsAtScaleStartIdx[jdScale]==176765){
//         //   std::cout<<" found 176765 in tree 0 with indices "<<index_biscale_j<<" "<<index_monoscale_j<<std::endl;
//         // }
//           coordinates_j = fitCoordinates[jdScale].data();
//           mlsCoordinates_j = mls3DCoordinates[jdScale].data();
//           // std::cout<<" jd here 1 "<<std::endl;
//           const float* coordsj = &coordinates_j[dim*J];


//           // std::cout<<" jd here 2"<<std::endl;
//           double dist = distanceFunction(coordsi, coordsj, &mlsCoordinates_i[3*I], &mlsCoordinates_j[3*J], pointCloudBarycenter.data());

//           if(dist <= maxEdgeLength){
//             // if(index_biscale_j==11091 or index_biscale_j==0 or index_biscale_i==0 or index_biscale_i==11091){
//             //   std::cout<<"                      edge "<<index_monoscale_i+pointsAtScaleStartIdx[idScale]<<" - "<<index_monoscale_j+pointsAtScaleStartIdx[jdScale]<<" :"<<dist<<std::endl;
//             // }
//             // if(I+pointsAtScaleStartIdx[idScale] ==0 or J+pointsAtScaleStartIdx[jdScale]==0){
//             // std::cout<<"                      edge "<<I+pointsAtScaleStartIdx[idScale]<<" - "<<J+pointsAtScaleStartIdx[jdScale]<<" :"<<dist<<std::endl;
//             // }
//             // if(I+pointsAtScaleStartIdx[idScale] ==11091 or J+pointsAtScaleStartIdx[jdScale]==11091){
//             // std::cout<<"                      edge "<<I+pointsAtScaleStartIdx[idScale]<<" - "<<J+pointsAtScaleStartIdx[jdScale]<<" :"<<dist<<std::endl;
//             // }
//             edges_.emplace_back(index_monoscale_i+pointsAtScaleStartIdx[idScale], index_monoscale_j+pointsAtScaleStartIdx[jdScale], dist);
//             numberOfEdges++;
//             addedEdges++;
//             if(dist>max_added_dist){
//               max_added_dist=dist;
//             }
//             if(dist<min_added_dist){
//               min_added_dist=dist;
//             }
//           }else{
//             // std::cout<<"edge "<<dist<<" discarded"<<std::endl;
//           }

//           // std::cout<<" done jd "<<std::endl;
//         }
//         // std::cout<<std::endl;
//         // std::cout<<" found "<<nMatches<<" matches w/ dist "<<distMax<<std::endl;
//       }
//       // for(LongSimplexId id=0; id<nbPointsAtLevel0; id++){
//       //   SimplexId index_biscale_i = tree0.samples()[id];
//       //   if(isMatched[index_biscale_i]==0){
//       //     std::cout<<"  ERROR: POINT "<<index_biscale_i<<" was not matched to any point in tree1"<<std::endl;
//       //   }
//       // }
//       std::cout<<" added "<<addedEdges<<" edges with max dist "<<max_added_dist<<" and min dist "<<min_added_dist<<std::endl;
//     }
//   }



//   std::cout<<" Found: "<<numberOfEdges<<" edges in "<<tm.getElapsedTime()<<" s."<<std::endl;
//   // std::cout<<" compared to : "<<numberOfEdges<<" effective edges"<<std::endl;

//   const auto cmp = [=](const edgeTuple &e0, const edgeTuple &e1){
//     return compareEdgeTuple(e0, e1);
//   };

//   TTK_PSORT(this->threadNumber_, edges_.begin(), edges_.end(), cmp);
//   printMsg("Done, found ");
// }

// template<class PointND>
// void RipsPersistence::addEpsilonPairs(std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSamplers,
//                                       std::vector<std::vector<SimplexId>>& filteredIdxToEpsSampling, 
//                                       std::vector<std::vector<SimplexId>>& filteredToLocalId,
//                                       std::vector<SimplexId>& pointsAtScaleStartIdx){

//   std::cout<<" Diagram size "<<diagram_.size()<<std::endl;
  
//   auto nbScales = epsilonSamplers.size();
//   if(filteredIdxToEpsSampling.size()!=nbScales or filteredToLocalId.size()!=nbScales){
//     printErr("Size mismatch");
//     return;
//   }

//   for(SimplexId iScale=0; iScale<nbScales; iScale++){
    
//     auto& epsReps = epsilonSamplers[iScale].representant();
//     float epsilon = epsilonSamplers[iScale].epsilon();

//     for(SimplexId id_filtered=0; id_filtered<epsReps.size(); id_filtered++){

//       SimplexId rep_filtered = epsReps[id_filtered]; // monoscale full filtered index
//                                            // if(rep_filtered>filteredIdxToEpsSampling[iScale].size()){
//                                            //   std::cout<<"size mistake "<<rep_filtered<<" "<<filteredToLocalIdVec[iScale].size()<<std::endl;
//                                            // }
//                                            // SimplexId rep_local = filteredToLocalIdVec[iScale][rep_filtered];
//       SimplexId birthId_ms = filteredIdxToEpsSampling[iScale][id_filtered] + pointsAtScaleStartIdx[iScale]; 
//       SimplexId rep_ms = filteredIdxToEpsSampling[iScale][rep_filtered] + pointsAtScaleStartIdx[iScale];

//       diagram_.emplace_back(PersistencePair{
//           CriticalVertex{birthId_ms, CriticalType::Local_minimum, 0,{}},
//           CriticalVertex{0, CriticalType::Local_maximum, epsilon,{}},
//           0, true, 1, rep_ms,{}});

//     }
//   }
//   std::cout<<" now is "<<diagram_.size()<<std::endl;
// }
