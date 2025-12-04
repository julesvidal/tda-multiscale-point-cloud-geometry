#include "DataTypes.h"
#include "Debug.h"
#include <RipsTree.h>
#include <algorithm>
#include <utility>
#include <vector>
#include <iterator>
#include <numeric>
#include <sstream>

using namespace ttk;

ttk::RipsTree::RipsTree() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("RipsTree");
  this->setDebugLevel(3);
}



// void RipsTree::buildArcMT(std::vector<SimplexId>& nodes, std::vector<SimplexId>& arcs){
//   for(auto& arc : ripsArcs_)  {
//     nodes.push_back(arc.birthId);
//     arcs.push_back()
//   }
// }

void RipsTree::buildArcDiagram(){

    for(size_t iArc=0; iArc<ripsArcs_.size(); iArc++){
      auto& arc = ripsArcs_[iArc];
      SimplexId birthId = this->getPointId(arc.ccId);
      double birth = arc.birth;
      double death = arc.death;
      int isFinite = (iArc!=0);
      arcDiagram_.emplace_back(PersistencePair{
          CriticalVertex{birthId, CriticalType::Local_minimum, birth,{}},
          CriticalVertex{-1, CriticalType::Local_maximum, death,{}},
          0, isFinite, 0, 0,{}});
    }
}

void RipsTree::compute5DAABBDiagonals(SimplexId iCC,
                             const std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                             const std::vector<std::vector<SimplexId>>& filteredToLocalId,
                             const std::vector<SimplexId>& pointsAtScaleStartIdx,
                             const std::vector<std::vector<double>>& fitCoordsVec,
                             const int dim,
                             std::vector<double>& output_aabb){
    
  auto& cc = this->getComponent(iCC);

  int d=dim;

  output_aabb.resize(2*d);


  auto& children = cc.getChildren();
  SimplexId nbChildren = children.size();

  SimplexId ms_id = cc.birthId();
  int iScale = findScaleId(ms_id, pointsAtScaleStartIdx);
  SimplexId id_mono_eps = ms_id - pointsAtScaleStartIdx[iScale];
  SimplexId local_id = filteredToLocalId[iScale][epsilonIdToFilteredId[iScale][id_mono_eps]];

  // output_aabb[0]   = fitCoordsVec[iScale][dim*local_id+1];
  // output_aabb[1] = fitCoordsVec[iScale][dim*local_id+1];

  for(int i=0; i<d; i++){
    output_aabb[2*i]   = fitCoordsVec[iScale][dim*local_id+dim-d+i];
    output_aabb[2*i+1] = fitCoordsVec[iScale][dim*local_id+dim-d+i];
  }

  // termination condition
  if(nbChildren==0){
    return;
  }

  // for(int i=0; i<dim; i++){
  //   output_aabb[2*i]=std::numeric_limits<double>::max();
  //   output_aabb[2*i+1]=std::numeric_limits<double>::lowest();
  // }
  
  // recursion 
  std::vector<double> diag_aabb_5D(nbChildren);

  for(int iChild=nbChildren-1; iChild>=0; iChild--){
    SimplexId child_ccId = children[iChild];
    std::vector<double> child_aabb(2*d);

    compute5DAABBDiagonals(child_ccId,
                           epsilonIdToFilteredId,
                           filteredToLocalId,
                           pointsAtScaleStartIdx,
                           fitCoordsVec,
                           dim, child_aabb);
    
    for(int i=0; i<d; i++){ // add bounding boxes
      if(child_aabb[2*i] <  output_aabb[2*i]){
        output_aabb[2*i]= child_aabb[2*i];
      }
      if(child_aabb[2*i+1] >  output_aabb[2*i+1]){
        output_aabb[2*i+1]= child_aabb[2*i+1];
      }
    }
    // emplace diag
    double diag = AABBToDiag(output_aabb, d);
    diag_aabb_5D[iChild] = diag;
  }
  cc.set_diag_aabb_5D(diag_aabb_5D);
  return;
}


void RipsTree::buildRipsArcsWithBirth( const std::vector<double>& births){


 ripsArcs_.clear();
 ripsArcs2_.clear();

 std::vector<RipsArc> ripsArcs2{};

 SimplexId nbCC = this->getNumberOfPoints();

 double min_persistence = 1e-7;
 int max_ccId = 1000;



  for(int iRoot=0; iRoot<tree_roots_.size(); iRoot++){
    
    SimplexId r = tree_roots_[iRoot];

    std::stack<std::tuple<SimplexId, SimplexId, SimplexId, double, SimplexId, double>> stack; // stack of pairs <iCC, i_current_split>


    stack.emplace(r, -1, this->getComponent(r).size(), this->getComponent(r).death(), r+nbCC, -1);

    while(!stack.empty()){
      auto iCC = std::get<0>(stack.top());
      // std::cout<<"icc "<<iCC<<std::endl;
      auto i_currentThreshold = std::get<1>(stack.top());
      auto size = std::get<2>(stack.top());
      double last_split_value = std::get<3>(stack.top());
      SimplexId last_split_id = std::get<4>(stack.top());
      stack.pop();


      auto& cc = this->getComponent(iCC);
      double birth = cc.birth();
      double persistence = cc.death()-cc.birth();

      // std::vector<SimplexId> ms_ids;


      if(persistence<min_persistence or iCC > max_ccId){
        continue;
      }

      double currentThreshold{};



      if(i_currentThreshold == -1){
        currentThreshold = cc.death();
      }else{
        currentThreshold = this->getComponent(cc.getChildren()[i_currentThreshold]).death();
      }


        bool create_leave_arc = true;
        if(i_currentThreshold+1 < cc.getChildren().size()){



          auto iChild = cc.getChildren()[i_currentThreshold+1];
          SimplexId child_size = this->getComponent(iChild).size();
          double child_death = this->getComponent(iChild).death();
          double child_birth = this->getComponent(iChild).birth();
          double child_persistence = child_death-child_birth;


          bool push_child = ( iChild < max_ccId and child_persistence>min_persistence);
          bool push_childless = true;
          create_leave_arc = !push_child and !push_childless;

          int arc_death_id = last_split_id;
          double arc_death = last_split_value;


          // process arcs for merge tree vis
          if(push_child and push_childless){
            ripsArcs2.push_back({iCC, child_death, last_split_value, iChild+nbCC, last_split_id, 2, iRoot, -1});
            arc_death = child_death;
            arc_death_id = iChild + nbCC;
          }

          if(push_childless){
            stack.emplace(iCC, i_currentThreshold+1, size-child_size, arc_death, arc_death_id, -1);
          }
          if(push_child){
            auto& child = this->getComponent(iChild);

            stack.emplace(iChild, -1, child_size, arc_death, arc_death_id, -1);
          }

        }

        if(create_leave_arc){ 
          int type = 0; 
          ripsArcs_.push_back({iCC, birth, last_split_value, iCC, last_split_id, type, iRoot, -1});
        }



    }
  }

 ripsArcs2_= std::move(ripsArcs2);

 // SimplexId nbNodes = ripsArcs_.size() + ripsArcs2_.size();
 SimplexId newNbCC = ripsArcs_.size();

 // std::cout<<" old nb cc "<<nbCC<<" to new nb cc "<<newNbCC<<std::endl;

 // std::cout<<"Leaves arc : "<<std::endl;
 for(SimplexId iArc = 0; iArc<ripsArcs_.size(); iArc++){
   auto& arc = ripsArcs_[iArc];
      // std::cout<<"   arc "<<iArc<<": "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<" type "<<arc.type<<std::endl;
 }
 // std::cout<<"Other arc : "<<std::endl;
 for(SimplexId iArc = 0; iArc<ripsArcs2_.size(); iArc++){
   auto& arc = ripsArcs2_[iArc];
      // std::cout<<"   arc "<<iArc<<": "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<" type "<<arc.type<<std::endl;
      // std::cout<<"                 : "<<arc.ccId<<" "<<arc.birthId-nbCC<<"-"<<arc.deathId-nbCC<<"("<<arc.birth<<"-"<<arc.death<<")"<<std::endl;
 }

 int test=0;
 // std::cout<<"test "<<test++<<std::endl;

 const auto cmp_pers = [=](const RipsArc &a0, const RipsArc &a1){
   return (a0.ccId<a1.ccId);
 };

 TTK_PSORT(this->threadNumber_, ripsArcs_.begin(), ripsArcs_.end(), cmp_pers);

 std::vector<SimplexId> oldIdToNewId(nbCC,-1);

 // new contiguous ids for components (leaf arcs)
 for(SimplexId i = 0; i<ripsArcs_.size(); i++){
  oldIdToNewId[ripsArcs_[i].birthId] = i;
 }

 // std::cout<<"test "<<test++<<std::endl;
 int extra_nodes=0;
 for(SimplexId i = 0; i<ripsArcs2_.size(); i++){
  SimplexId bid = ripsArcs2_[i].birthId-nbCC;
  SimplexId did = ripsArcs2_[i].deathId-nbCC;
  if(oldIdToNewId[bid] == -1){
    oldIdToNewId[bid] = extra_nodes+2*newNbCC;
    extra_nodes++;
    // std::cout<<" new node bid "<<bid<<std::endl;
  }
  if(oldIdToNewId[did] == -1){
    oldIdToNewId[did] = extra_nodes+2*newNbCC;
    extra_nodes++;
    // std::cout<<" new node did "<<did<<std::endl;;
  }
 }
 // std::cout<<" Extra nodes : "<<extra_nodes<<std::endl;
 // std::cout<<"test "<<test++<<std::endl;

 for(SimplexId i = 0; i<ripsArcs_.size(); i++){
  ripsArcs_[i].birthId = oldIdToNewId[ripsArcs_[i].birthId];
  ripsArcs_[i].deathId = oldIdToNewId[ripsArcs_[i].deathId-nbCC]+newNbCC;
 }

 for(SimplexId i = 0; i<ripsArcs2_.size(); i++){
  ripsArcs2_[i].birthId = oldIdToNewId[ripsArcs2_[i].birthId-nbCC] + newNbCC;
  ripsArcs2_[i].deathId = oldIdToNewId[ripsArcs2_[i].deathId-nbCC] + newNbCC;
 }



 const auto cmp = [=](const RipsArc &a0, const RipsArc &a1){
   return ((a0.death-a0.birth) > (a1.death-a1.birth)) || (((a0.death-a0.birth) == (a1.death-a1.birth)) && a0.ccId<a1.ccId);
 };

 TTK_PSORT(this->threadNumber_, ripsArcs_.begin(), ripsArcs_.end(), cmp);

 std::cout<<" found "<<ripsArcs_.size()<<" leaves arcs"<<std::endl;
 std::cout<<" and found "<<ripsArcs2_.size()<<" branching arcs"<<std::endl;


ripsArcs_.insert(
      ripsArcs_.end(),
      std::make_move_iterator(ripsArcs2_.begin()),
      std::make_move_iterator(ripsArcs2_.end()));



 TTK_PSORT(this->threadNumber_, ripsArcs_.begin(), ripsArcs_.end(), cmp);

 std::cout<<" Found "<<ripsArcs_.size()<<" rips arcs"<<std::endl;
}
void RipsTree::buildAABBDiagonals(const std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                  const std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                  const std::vector<SimplexId>& pointsAtScaleStartIdx,
                                  const std::vector<std::vector<double>>& fitCoordsVec,
                                  const int dim){

  for(auto r : tree_roots_){
    std::vector<double> aabb{};
    compute5DAABBDiagonals(r,
        epsilonIdToFilteredId,
        filteredToLocalId,
        pointsAtScaleStartIdx,
        fitCoordsVec,
        dim,
        aabb);
  }
}

void RipsTree::buildRipsArcsWithBirth(SimplexId sizeMin,
                             double diagMax,
                             const std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                             const std::vector<std::vector<SimplexId>>& filteredToLocalId,
                             const std::vector<SimplexId>& pointsAtScaleStartIdx,
                             const std::vector<std::vector<double>>& fitCoordsVec,
                             const std::vector<double>& births,
                             const int dim){

 std::cout<<" build rips arcs with sizemin="<<sizeMin<<" and diagmax="<<diagMax<<std::endl;

 ripsArcs_.clear();
 ripsArcs2_.clear();

 std::vector<RipsArc> ripsArcs2{};

 SimplexId nbCC = this->getNumberOfPoints();

 double min_persistence = 1e-7;
 int max_ccId = 1000;

  // std::cout<<" nbCC "<<nbCC<<std::endl;


  for(int iRoot=0; iRoot<tree_roots_.size(); iRoot++){
    
    SimplexId r = tree_roots_[iRoot];

    std::stack<std::tuple<SimplexId, SimplexId, SimplexId, double, SimplexId, double>> stack; // stack of pairs <iCC, i_current_split>

    std::vector<double> aabb{};
    // std::cout<<"compute 5d diags for "<<r<<std::endl;
    compute5DAABBDiagonals(r,
        epsilonIdToFilteredId,
        filteredToLocalId,
        pointsAtScaleStartIdx,
        fitCoordsVec,
        dim,
        aabb);
    // std::cout<<"done"<<std::endl;

    stack.emplace(r, -1, this->getComponent(r).size(), this->getComponent(r).death(), r+nbCC, -1);

    while(!stack.empty()){
      auto iCC = std::get<0>(stack.top());
      // std::cout<<"icc "<<iCC<<std::endl;
      auto i_currentThreshold = std::get<1>(stack.top());
      auto size = std::get<2>(stack.top());
      double last_split_value = std::get<3>(stack.top());
      SimplexId last_split_id = std::get<4>(stack.top());
      double last_split_diag = std::get<5>(stack.top());
      stack.pop();


      auto& cc = this->getComponent(iCC);
      double birth = cc.birth();
      double persistence = cc.death()-cc.birth();

      // std::vector<SimplexId> ms_ids;


      if(size<sizeMin or persistence<min_persistence or iCC > max_ccId){
        continue;
      }

      double currentThreshold{};



      if(i_currentThreshold == -1){
        currentThreshold = cc.death();
      }else{
        // if(i_currentThreshold >= cc.getChildren().size()){
        //   std::cout<<" ERROR children: "<<i_currentThreshold<<" outta "<<cc.getChildren().size()<<std::endl;
        // }
        currentThreshold = this->getComponent(cc.getChildren()[i_currentThreshold]).death();
      }

      // if(i_currentThreshold+1 >= cc.diag_aabb_5D().size()){
      //   std::cout<<" ERROR diag: "<<i_currentThreshold+1<<" outta "<<cc.diag_aabb_5D().size()<<std::endl;
      // }

      double diag = i_currentThreshold+1 < cc.diag_aabb_5D().size() ? cc.diag_aabb_5D()[i_currentThreshold+1] : 0.0;

      // std::cout<<" pop "<<iCC<<" at "<<currentThreshold<<", i_currentThreshold="<<i_currentThreshold<<", children_size="<<cc.getChildren().size()<<", diag_size="<<cc.diag_aabb_5D().size()<<std::endl;
      // std::cout<<"                                          s="<<size<<", d="<<diag<<std::endl;

      if(diag<diagMax){

        if(i_currentThreshold>-1 and cc.getChildren()[i_currentThreshold]+nbCC != last_split_id){
          int dummy_split_child = cc.getChildren()[i_currentThreshold];

          ripsArcs_.push_back({iCC, birth, currentThreshold, iCC, dummy_split_child+nbCC, 0, iRoot, diag});
          // std::cout<<"   -> create arc "<<iCC<<" "<<birth<<"-"<<currentThreshold<<" of type "<<0<<std::endl;
          // std::cout<<"          bid "<<iCC<<"  did"<<dummy_split_child+nbCC<<std::endl;
          //   if(iCC==dummy_split_child){
          //     std::cout<<"ERROR HERE"<<std::endl;
          //   }

          double diag_child = this->getComponent(dummy_split_child).getDiag();
          double birth_child = this->getComponent(dummy_split_child).birth();

          ripsArcs_.push_back({dummy_split_child, birth_child, currentThreshold, dummy_split_child, dummy_split_child+nbCC, 1, iRoot, diag_child});
          // std::cout<<"   -> create arc "<<dummy_split_child<<" "<<birth_child<<"-"<<currentThreshold<<" w/ diag "<<diag_child<<" of type "<<1<<std::endl;
          // std::cout<<"          bid "<<dummy_split_child<<"  did"<<dummy_split_child+nbCC<<std::endl;
          //   if(dummy_split_child==dummy_split_child+nbCC){
          //     std::cout<<"ERROR HERE"<<std::endl;
          //   }

          double diag_trunk = last_split_diag;
          ripsArcs2.push_back({iCC, currentThreshold, last_split_value, dummy_split_child+nbCC, last_split_id, 2, iRoot, diag_trunk});
          // std::cout<<"   -> create arc "<<iCC<<" "<<currentThreshold<<"-"<<last_split_value<<" w/ diag "<<diag_trunk<<" of type "<<2<<std::endl;
          // std::cout<<"          bid "<<dummy_split_child+nbCC<<"  did "<<last_split_id<<std::endl;
          //   if(dummy_split_child+nbCC==last_split_id){
          //     std::cout<<"ERROR HERE"<<std::endl;
          //   }

        }else{
          ripsArcs_.push_back({iCC, birth, currentThreshold, iCC, last_split_id, 0, iRoot, diag});
          // std::cout<<"   -> create arc "<<iCC<<" "<<birth<<"-"<<currentThreshold<<" of size "<<size<<" w/ diag "<<diag<<" and type "<<0<<std::endl;
          // std::cout<<"          bid "<<iCC<<"  did"<<last_split_id<<std::endl;
          //   if(iCC==last_split_id){
          //     std::cout<<"ERROR HERE"<<std::endl;
          //   }
        }

      }else{


        bool create_leave_arc = true;
        if(i_currentThreshold+1 < cc.getChildren().size()){

          // split: add split arc
          // arc iChild+nbArcs -> last_split_id
          // child.death -> last_split_value


          auto iChild = cc.getChildren()[i_currentThreshold+1];
          SimplexId child_size = this->getComponent(iChild).size();
          double child_death = this->getComponent(iChild).death();
          double child_birth = this->getComponent(iChild).birth();
          double child_persistence = child_death-child_birth;


          bool push_child = (child_size >= sizeMin and iChild < max_ccId and child_persistence>min_persistence);
          bool push_childless = (size-child_size >= sizeMin);
          create_leave_arc = !push_child and !push_childless;

          int arc_death_id = last_split_id;
          double arc_death = last_split_value;


          // process arcs for merge tree vis
          if(push_child and push_childless){
            // real split: add "split arc"
            ripsArcs2.push_back({iCC, child_death, last_split_value, iChild+nbCC, last_split_id, 2, iRoot, diag});
            // std::cout<<"   -> create arc "<<iCC<<" "<<child_death<<"-"<<last_split_value<<" w/ diag "<<diag<<" and type "<<2<<std::endl;
            // std::cout<<"          bid "<<iChild+nbCC<<"  did"<<last_split_id<<std::endl;
            // if(iChild+nbCC==last_split_id){
            //   std::cout<<"ERROR HERE"<<std::endl;
            // }
            arc_death = child_death;
            arc_death_id = iChild + nbCC;
          }

          // std::cout<<" childless "<<iCC <<" size "<<size-child_size<<" vs sizemin "<<sizeMin<<std::endl;
          if(push_childless){
            double diag_childless = i_currentThreshold+2 < cc.diag_aabb_5D().size() ? cc.diag_aabb_5D()[i_currentThreshold+2] : 0.0;
            stack.emplace(iCC, i_currentThreshold+1, size-child_size, arc_death, arc_death_id, diag);
            // std::cout<<"   -> split into "<<iCC<<" at "<<child_death<<", s="<<size-child_size<<std::endl;
          }
          // std::cout<<" child "<<iChild <<" size "<<child_size<<" vs sizemin "<<sizeMin<<std::endl;
          if(push_child){
            auto& child = this->getComponent(iChild);
            double diag_child = 0 < child.diag_aabb_5D().size() ? child.diag_aabb_5D()[0] : 0.0;

            stack.emplace(iChild, -1, child_size, arc_death, arc_death_id, diag);
            // std::cout<<"   -> split into child "<<iChild<<" at "<<child_death<<", s="<<child_size<<std::endl;
          }

        }

        if(create_leave_arc){ // add leave arc with diag > diagmax
          int type = diagMax < 0 ? 0 : 1; // treat this arc as accepted arc whenever diagmax is negative
          ripsArcs_.push_back({iCC, birth, last_split_value, iCC, last_split_id, type, iRoot, diag});
          // std::cout<<"   -> create arc "<<iCC<<" "<<birth<<"-"<<currentThreshold<<" of size "<<size<<" w/ diag "<<diag<<" and type "<<type<<std::endl;
          // std::cout<<"          bid "<<iCC<<"  did"<<last_split_id<<std::endl;
          // if(birth>last_split_value){
            // std::cout<<" /!\ ERROR persistence negative"<<std::endl;
          // }
          //   if(iCC==last_split_id){
          //     std::cout<<"ERROR HERE"<<std::endl;
          //   }
        }
      }
    }
    // std::cout<<" stack is empty ! "<<std::endl;
  }

 // clean intermediate branches
 // SimplexId current_icc = -1;
 // RipsArc currentArc = ripsArcs2[0];
 // std::cout<<"   arc "<<currentArc.ccId<<" "<<currentArc.birthId<<"-"<<currentArc.deathId<<"("<<currentArc.birth<<"-"<<currentArc.death<<")"<<std::endl;
 // for(int iArc = 1; iArc<ripsArcs2.size(); iArc++){
 //    auto arc = ripsArcs2[iArc];
 //    if(arc.ccId == currentArc.ccId){
 //      currentArc.birthId = arc.birthId;
 //      currentArc.birth = arc.birth;
 //      std::cout<<"   arc "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<std::endl;
 //    }else{
 //      ripsArcs2_.push_back(currentArc);
 //      std::cout<<" sum up by arc "<<currentArc.ccId<<" "<<currentArc.birthId<<"-"<<currentArc.deathId<<"("<<currentArc.birth<<"-"<<currentArc.death<<")"<<std::endl;
 //      currentArc = arc;
 //      std::cout<<" current arc "<<currentArc.ccId<<" "<<currentArc.birthId<<"-"<<currentArc.deathId<<"("<<currentArc.birth<<"-"<<currentArc.death<<")"<<std::endl;
 //      std::cout<<"   new arc "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<std::endl;
 //    }
 // }
 //  ripsArcs2_.push_back(currentArc);
 //  std::cout<<" sum up by arc "<<currentArc.ccId<<" "<<currentArc.birthId<<"-"<<currentArc.deathId<<"("<<currentArc.birth<<"-"<<currentArc.death<<")"<<std::endl;
 ripsArcs2_= std::move(ripsArcs2);

 // std::cout<<"\nLeaves arc : "<<std::endl;
 // for(SimplexId iArc = 0; iArc<ripsArcs_.size(); iArc++){
 //   auto& arc = ripsArcs_[iArc];
 //      std::cout<<"   arc "<<iArc<<": "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<" type "<<arc.type<<std::endl;
 // }
 // std::cout<<"Other arc : "<<std::endl;
 // for(SimplexId iArc = 0; iArc<ripsArcs2_.size(); iArc++){
 //   auto& arc = ripsArcs2_[iArc];
 //      std::cout<<"   arc "<<iArc<<": "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<" type "<<arc.type<<std::endl;
 //      std::cout<<"                 : "<<arc.ccId<<" "<<arc.birthId-nbCC<<"-"<<arc.deathId-nbCC<<"("<<arc.birth<<"-"<<arc.death<<")"<<std::endl;
 // }
 // std::cout<<std::endl;

  // ================= clean dangling branches ===============

 // std::cout<<"\ncleaning dangling branches"<<std::endl;
 // const auto cmp_inc_birth = [=](const RipsArc &a0, const RipsArc &a1){
 //   return ((a0.birth) < (a1.birth)) || (((a0.birth) == (a1.birth)) && a0.ccId<a1.ccId);
 // };

 //  TTK_PSORT(this->threadNumber_, ripsArcs2_.begin(), ripsArcs2_.end(), cmp_inc_birth);
 //  std::vector<RipsArc> cleaned_ripsArcs2{};
 //  std::vector<int> match_count(nbCC,0);
 //  for(auto& arc : ripsArcs_){
 //    auto dId = arc.deathId - nbCC;
 //    match_count[dId]++;
 //  }
 //  for(auto arc : ripsArcs2_){
 //    auto bId = arc.birthId - nbCC;
 //    auto dId = arc.deathId - nbCC;
 //    if(match_count[bId]){ // a leave connects this arc by the bottom -> keepit
 //      match_count[dId]++;
 //    }
 //  }
 //  for(auto arc : ripsArcs2_){
 //    auto bId = arc.birthId - nbCC;
 //    if(match_count[bId]){ // a leave connects this arc by the bottom -> keepit
 //      cleaned_ripsArcs2.push_back(arc);
 //    }
 //  }
 //  ripsArcs2_ = std::move(cleaned_ripsArcs2);

 // std::cout<<"Other arcs after cleaning : "<<std::endl;
 // for(SimplexId iArc = 0; iArc<ripsArcs2_.size(); iArc++){
 //   auto& arc = ripsArcs2_[iArc];
 //      std::cout<<"   arc "<<iArc<<": "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<" type "<<arc.type<<std::endl;
 //      std::cout<<"                 : "<<arc.ccId<<" "<<arc.birthId-nbCC<<"-"<<arc.deathId-nbCC<<"("<<arc.birth<<"-"<<arc.death<<")"<<std::endl;
 // }
 // std::cout<<std::endl;
  // =========================================================

 // SimplexId nbNodes = ripsArcs_.size() + ripsArcs2_.size();
 SimplexId newNbCC = ripsArcs_.size();

 // std::cout<<" old nb cc "<<nbCC<<" to new nb cc "<<newNbCC<<std::endl;

 // std::cout<<"Leaves arc : "<<std::endl;
 for(SimplexId iArc = 0; iArc<ripsArcs_.size(); iArc++){
   auto& arc = ripsArcs_[iArc];
      // std::cout<<"   arc "<<iArc<<": "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<" type "<<arc.type<<std::endl;
 }
 // std::cout<<"Other arc : "<<std::endl;
 for(SimplexId iArc = 0; iArc<ripsArcs2_.size(); iArc++){
   auto& arc = ripsArcs2_[iArc];
      // std::cout<<"   arc "<<iArc<<": "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<" type "<<arc.type<<std::endl;
      // std::cout<<"                 : "<<arc.ccId<<" "<<arc.birthId-nbCC<<"-"<<arc.deathId-nbCC<<"("<<arc.birth<<"-"<<arc.death<<")"<<std::endl;
 }

 int test=0;
 // std::cout<<"test "<<test++<<std::endl;

 const auto cmp_pers = [=](const RipsArc &a0, const RipsArc &a1){
   return (a0.ccId<a1.ccId);
 };

 TTK_PSORT(this->threadNumber_, ripsArcs_.begin(), ripsArcs_.end(), cmp_pers);

 std::vector<SimplexId> oldIdToNewId(nbCC,-1);

 // new contiguous ids for components (leaf arcs)
 for(SimplexId i = 0; i<ripsArcs_.size(); i++){
  oldIdToNewId[ripsArcs_[i].birthId] = i;
 }

 // std::cout<<"test "<<test++<<std::endl;
 int extra_nodes=0;
 for(SimplexId i = 0; i<ripsArcs2_.size(); i++){
  SimplexId bid = ripsArcs2_[i].birthId-nbCC;
  SimplexId did = ripsArcs2_[i].deathId-nbCC;
  if(oldIdToNewId[bid] == -1){
    oldIdToNewId[bid] = extra_nodes+2*newNbCC;
    extra_nodes++;
    // std::cout<<" new node bid "<<bid<<std::endl;
  }
  if(oldIdToNewId[did] == -1){
    oldIdToNewId[did] = extra_nodes+2*newNbCC;
    extra_nodes++;
    // std::cout<<" new node did "<<did<<std::endl;;
  }
 }
 // std::cout<<" Extra nodes : "<<extra_nodes<<std::endl;
 // std::cout<<"test "<<test++<<std::endl;

 for(SimplexId i = 0; i<ripsArcs_.size(); i++){
   // if(oldIdToNewId[ripsArcs_[i].birthId]==-1 or oldIdToNewId[ripsArcs_[i].deathId-nbCC]==-1){
     // std::cout<<"shit here"<<oldIdToNewId[ripsArcs_[i].birthId]<<" "<<oldIdToNewId[ripsArcs_[i].deathId-nbCC]<<std::endl;
   // }
  ripsArcs_[i].birthId = oldIdToNewId[ripsArcs_[i].birthId];
  ripsArcs_[i].deathId = oldIdToNewId[ripsArcs_[i].deathId-nbCC]+newNbCC;
 }
 // std::cout<<"test "<<test++<<std::endl;

 for(SimplexId i = 0; i<ripsArcs2_.size(); i++){
  ripsArcs2_[i].birthId = oldIdToNewId[ripsArcs2_[i].birthId-nbCC] + newNbCC;
  ripsArcs2_[i].deathId = oldIdToNewId[ripsArcs2_[i].deathId-nbCC] + newNbCC;
 }



 const auto cmp = [=](const RipsArc &a0, const RipsArc &a1){
   return ((a0.death-a0.birth) > (a1.death-a1.birth)) || (((a0.death-a0.birth) == (a1.death-a1.birth)) && a0.ccId<a1.ccId);
 };

 TTK_PSORT(this->threadNumber_, ripsArcs_.begin(), ripsArcs_.end(), cmp);

 std::cout<<" found "<<ripsArcs_.size()<<" leaves arcs"<<std::endl;
 std::cout<<" and found "<<ripsArcs2_.size()<<" branching arcs"<<std::endl;

 // std::cout<<"Leaves arc : "<<std::endl;
 // for(SimplexId iArc = 0; iArc<ripsArcs_.size(); iArc++){
 //   auto& arc = ripsArcs_[iArc];
 //      std::cout<<"   arc "<<iArc<<": "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<" type "<<arc.type<<std::endl;
 // }
 // std::cout<<"Other arc : "<<std::endl;
 // for(SimplexId iArc = 0; iArc<ripsArcs2_.size(); iArc++){
 //   auto& arc = ripsArcs2_[iArc];
 //      std::cout<<"   arc "<<iArc<<": "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<" type "<<arc.type<<std::endl;
 // }
ripsArcs_.insert(
      ripsArcs_.end(),
      std::make_move_iterator(ripsArcs2_.begin()),
      std::make_move_iterator(ripsArcs2_.end()));
 // std::cout<<" now size is "<<ripsArcs_.size()<<" full arcs"<<std::endl;

 // for(auto& arc : ripsArcs_){
 //   std::cout<<" arc "<<arc.ccId<<" "<<arc.death-arc.birth<<std::endl;
 // }
 // for(SimplexId iArc = 0; iArc< ripsArcs_.size(); iArc++){

 //       auto ccId = ripsArcs_[iArc].ccId;
 //       auto threshold = ripsArcs_[iArc].death;
 //       auto birth = ripsArcs_[iArc].birth;
 //       auto persistence = ripsArcs_[iArc].death-ripsArcs_[iArc].birth;
 //       auto type = ripsArcs_[iArc].type;

 //     std::cout<<"    arc  "<<iArc<<": "<<ccId<<" "<<birth<<"-"<<threshold<<" of size "<<sizeAtThreshold(ccId, threshold)<<std::endl;
 // }
  
 // decrease the infinite arc a bit
 // ripsArcs_[0].death = ripsArcs_[1].death*2;


 TTK_PSORT(this->threadNumber_, ripsArcs_.begin(), ripsArcs_.end(), cmp);

 std::cout<<" Found "<<ripsArcs_.size()<<" rips arcs"<<std::endl;
}
void RipsTree::buildRipsArcs(SimplexId sizeMin,
                             double diagMax,
                             const std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                             const std::vector<std::vector<SimplexId>>& filteredToLocalId,
                             const std::vector<SimplexId>& pointsAtScaleStartIdx,
                             const std::vector<std::vector<double>>& fitCoordsVec,
                             const int dim){

 std::cout<<" build rips arcs with sizemin="<<sizeMin<<" and diagmax="<<diagMax<<std::endl;

 ripsArcs_.clear();
 ripsArcs2_.clear();

 std::vector<RipsArc> ripsArcs2{};

 SimplexId nbCC = this->getNumberOfPoints();

  // std::cout<<" nbCC "<<nbCC<<std::endl;


  for(int iRoot=0; iRoot<tree_roots_.size(); iRoot++){
    
    SimplexId r = tree_roots_[iRoot];

    std::stack<std::tuple<SimplexId, SimplexId, SimplexId, double, SimplexId, double>> stack; // stack of pairs <iCC, i_current_split>

    std::vector<double> aabb{};
    std::cout<<"compute 5d diags for "<<r<<std::endl;
    compute5DAABBDiagonals(r,
        epsilonIdToFilteredId,
        filteredToLocalId,
        pointsAtScaleStartIdx,
        fitCoordsVec,
        dim,
        aabb);
    std::cout<<"done"<<std::endl;

    stack.emplace(r, -1, this->getComponent(r).size(), this->getComponent(r).death(), r+nbCC, -1);

    while(!stack.empty()){
      auto iCC = std::get<0>(stack.top());
      // std::cout<<"icc "<<iCC<<std::endl;
      auto i_currentThreshold = std::get<1>(stack.top());
      auto size = std::get<2>(stack.top());
      double last_split_value = std::get<3>(stack.top());
      SimplexId last_split_id = std::get<4>(stack.top());
      double last_split_diag = std::get<5>(stack.top());
      stack.pop();


      auto& cc = this->getComponent(iCC);

      // std::vector<SimplexId> ms_ids;


      if(size<sizeMin){
        continue;
      }

      double currentThreshold{};



      if(i_currentThreshold == -1){
        currentThreshold = cc.death();
      }else{
        // if(i_currentThreshold >= cc.getChildren().size()){
        //   std::cout<<" ERROR children: "<<i_currentThreshold<<" outta "<<cc.getChildren().size()<<std::endl;
        // }
        currentThreshold = this->getComponent(cc.getChildren()[i_currentThreshold]).death();
      }

      // if(i_currentThreshold+1 >= cc.diag_aabb_5D().size()){
      //   std::cout<<" ERROR diag: "<<i_currentThreshold+1<<" outta "<<cc.diag_aabb_5D().size()<<std::endl;
      // }

      double diag = i_currentThreshold+1 < cc.diag_aabb_5D().size() ? cc.diag_aabb_5D()[i_currentThreshold+1] : 0.0;

      // std::cout<<" pop "<<iCC<<" at "<<currentThreshold<<", i_currentThreshold="<<i_currentThreshold<<", children_size="<<cc.getChildren().size()<<", diag_size="<<cc.diag_aabb_5D().size()<<std::endl;
      // std::cout<<"                                          s="<<size<<", d="<<diag<<std::endl;

      if(diag<diagMax){

        if(i_currentThreshold>-1 and cc.getChildren()[i_currentThreshold]+nbCC != last_split_id){
          int dummy_split_child = cc.getChildren()[i_currentThreshold];

          ripsArcs_.push_back({iCC, 0, currentThreshold, iCC, dummy_split_child+nbCC, 0, iRoot, diag});
          // std::cout<<"   -> create arc "<<iCC<<" "<<0<<"-"<<currentThreshold<<" of type "<<0<<std::endl;
          // std::cout<<"          bid "<<iCC<<"  did"<<dummy_split_child+nbCC<<std::endl;
          //   if(iCC==dummy_split_child){
          //     std::cout<<"ERROR HERE"<<std::endl;
          //   }

          double diag_child = this->getComponent(dummy_split_child).getDiag();

          ripsArcs_.push_back({dummy_split_child, 0, currentThreshold, dummy_split_child, dummy_split_child+nbCC, 1, iRoot, diag_child});
          // std::cout<<"   -> create arc "<<dummy_split_child<<" "<<0<<"-"<<currentThreshold<<" w/ diag "<<diag_child<<" of type "<<1<<std::endl;
          // std::cout<<"          bid "<<dummy_split_child<<"  did"<<dummy_split_child+nbCC<<std::endl;
          //   if(dummy_split_child==dummy_split_child+nbCC){
          //     std::cout<<"ERROR HERE"<<std::endl;
          //   }

          double diag_trunk = last_split_diag;
          ripsArcs2.push_back({iCC, currentThreshold, last_split_value, dummy_split_child+nbCC, last_split_id, 2, iRoot, diag_trunk});
          // std::cout<<"   -> create arc "<<iCC<<" "<<currentThreshold<<"-"<<last_split_value<<" w/ diag "<<diag_trunk<<" of type "<<2<<std::endl;
          // std::cout<<"          bid "<<dummy_split_child+nbCC<<"  did "<<last_split_id<<std::endl;
          //   if(dummy_split_child+nbCC==last_split_id){
          //     std::cout<<"ERROR HERE"<<std::endl;
          //   }

        }else{
          ripsArcs_.push_back({iCC, 0, currentThreshold, iCC, last_split_id, 0, iRoot, diag});
          // std::cout<<"   -> create arc "<<iCC<<" "<<0<<"-"<<currentThreshold<<" of size "<<size<<" w/ diag "<<diag<<" and type "<<0<<std::endl;
          // std::cout<<"          bid "<<iCC<<"  did"<<last_split_id<<std::endl;
          //   if(iCC==last_split_id){
          //     std::cout<<"ERROR HERE"<<std::endl;
          //   }
        }

      }else{


        bool create_leave_arc = true;
        if(i_currentThreshold+1 < cc.getChildren().size()){

          // split: add split arc
          // arc iChild+nbArcs -> last_split_id
          // child.death -> last_split_value


          auto iChild = cc.getChildren()[i_currentThreshold+1];
          SimplexId child_size = this->getComponent(iChild).size();
          double child_death = this->getComponent(iChild).death();


          bool push_child = (child_size >= sizeMin);
          bool push_childless = (size-child_size >= sizeMin);
          create_leave_arc = !push_child and !push_childless;

          int arc_death_id = last_split_id;
          double arc_death = last_split_value;


          // process arcs for merge tree vis
          if(push_child and push_childless){
            // real split: add "split arc"
            ripsArcs2.push_back({iCC, child_death, last_split_value, iChild+nbCC, last_split_id, 2, iRoot, diag});
            // std::cout<<"   -> create arc "<<iCC<<" "<<child_death<<"-"<<last_split_value<<" w/ diag "<<diag<<" and type "<<2<<std::endl;
            // std::cout<<"          bid "<<iChild+nbCC<<"  did"<<last_split_id<<std::endl;
            // if(iChild+nbCC==last_split_id){
            //   std::cout<<"ERROR HERE"<<std::endl;
            // }
            arc_death = child_death;
            arc_death_id = iChild + nbCC;
          }

          if(push_childless){
            double diag_childless = i_currentThreshold+2 < cc.diag_aabb_5D().size() ? cc.diag_aabb_5D()[i_currentThreshold+2] : 0.0;
            stack.emplace(iCC, i_currentThreshold+1, size-child_size, arc_death, arc_death_id, diag);
            // std::cout<<"   -> split into "<<iCC<<" at "<<child_death<<", s="<<size-child_size<<std::endl;
          }
          if(push_child){
            auto& child = this->getComponent(iChild);
            double diag_child = 0 < child.diag_aabb_5D().size() ? child.diag_aabb_5D()[0] : 0.0;

            stack.emplace(iChild, -1, child_size, arc_death, arc_death_id, diag);
            // std::cout<<"   -> split into child "<<iChild<<" at "<<child_death<<", s="<<child_size<<std::endl;
          }

        }

        if(create_leave_arc){ // add leave arc with diag > diagmax
          ripsArcs_.push_back({iCC, 0, last_split_value, iCC, last_split_id, 1, iRoot, diag});
          // std::cout<<"   -> create arc "<<iCC<<" "<<0<<"-"<<currentThreshold<<" of size "<<size<<" w/ diag "<<diag<<" and type "<<1<<std::endl;
          // std::cout<<"          bid "<<iCC<<"  did"<<last_split_id<<std::endl;
          //   if(iCC==last_split_id){
          //     std::cout<<"ERROR HERE"<<std::endl;
          //   }
        }
      }
    }
    // std::cout<<" stack is empty ! "<<std::endl;
  }

 // clean intermediate branches
 // SimplexId current_icc = -1;
 // RipsArc currentArc = ripsArcs2[0];
 // std::cout<<"   arc "<<currentArc.ccId<<" "<<currentArc.birthId<<"-"<<currentArc.deathId<<"("<<currentArc.birth<<"-"<<currentArc.death<<")"<<std::endl;
 // for(int iArc = 1; iArc<ripsArcs2.size(); iArc++){
 //    auto arc = ripsArcs2[iArc];
 //    if(arc.ccId == currentArc.ccId){
 //      currentArc.birthId = arc.birthId;
 //      currentArc.birth = arc.birth;
 //      std::cout<<"   arc "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<std::endl;
 //    }else{
 //      ripsArcs2_.push_back(currentArc);
 //      std::cout<<" sum up by arc "<<currentArc.ccId<<" "<<currentArc.birthId<<"-"<<currentArc.deathId<<"("<<currentArc.birth<<"-"<<currentArc.death<<")"<<std::endl;
 //      currentArc = arc;
 //      std::cout<<" current arc "<<currentArc.ccId<<" "<<currentArc.birthId<<"-"<<currentArc.deathId<<"("<<currentArc.birth<<"-"<<currentArc.death<<")"<<std::endl;
 //      std::cout<<"   new arc "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<std::endl;
 //    }
 // }
 //  ripsArcs2_.push_back(currentArc);
 //  std::cout<<" sum up by arc "<<currentArc.ccId<<" "<<currentArc.birthId<<"-"<<currentArc.deathId<<"("<<currentArc.birth<<"-"<<currentArc.death<<")"<<std::endl;
 ripsArcs2_= std::move(ripsArcs2);

 // std::cout<<"\nLeaves arc : "<<std::endl;
 // for(SimplexId iArc = 0; iArc<ripsArcs_.size(); iArc++){
 //   auto& arc = ripsArcs_[iArc];
 //      std::cout<<"   arc "<<iArc<<": "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<" type "<<arc.type<<std::endl;
 // }
 // std::cout<<"Other arc : "<<std::endl;
 // for(SimplexId iArc = 0; iArc<ripsArcs2_.size(); iArc++){
 //   auto& arc = ripsArcs2_[iArc];
 //      std::cout<<"   arc "<<iArc<<": "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<" type "<<arc.type<<std::endl;
 //      std::cout<<"                 : "<<arc.ccId<<" "<<arc.birthId-nbCC<<"-"<<arc.deathId-nbCC<<"("<<arc.birth<<"-"<<arc.death<<")"<<std::endl;
 // }
 // std::cout<<std::endl;

  // ================= clean dangling branches ===============

 // std::cout<<"\ncleaning dangling branches"<<std::endl;
 // const auto cmp_inc_birth = [=](const RipsArc &a0, const RipsArc &a1){
 //   return ((a0.birth) < (a1.birth)) || (((a0.birth) == (a1.birth)) && a0.ccId<a1.ccId);
 // };

 //  TTK_PSORT(this->threadNumber_, ripsArcs2_.begin(), ripsArcs2_.end(), cmp_inc_birth);
 //  std::vector<RipsArc> cleaned_ripsArcs2{};
 //  std::vector<int> match_count(nbCC,0);
 //  for(auto& arc : ripsArcs_){
 //    auto dId = arc.deathId - nbCC;
 //    match_count[dId]++;
 //  }
 //  for(auto arc : ripsArcs2_){
 //    auto bId = arc.birthId - nbCC;
 //    auto dId = arc.deathId - nbCC;
 //    if(match_count[bId]){ // a leave connects this arc by the bottom -> keepit
 //      match_count[dId]++;
 //    }
 //  }
 //  for(auto arc : ripsArcs2_){
 //    auto bId = arc.birthId - nbCC;
 //    if(match_count[bId]){ // a leave connects this arc by the bottom -> keepit
 //      cleaned_ripsArcs2.push_back(arc);
 //    }
 //  }
 //  ripsArcs2_ = std::move(cleaned_ripsArcs2);

 // std::cout<<"Other arcs after cleaning : "<<std::endl;
 // for(SimplexId iArc = 0; iArc<ripsArcs2_.size(); iArc++){
 //   auto& arc = ripsArcs2_[iArc];
 //      std::cout<<"   arc "<<iArc<<": "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<" type "<<arc.type<<std::endl;
 //      std::cout<<"                 : "<<arc.ccId<<" "<<arc.birthId-nbCC<<"-"<<arc.deathId-nbCC<<"("<<arc.birth<<"-"<<arc.death<<")"<<std::endl;
 // }
 // std::cout<<std::endl;
  // =========================================================

 // SimplexId nbNodes = ripsArcs_.size() + ripsArcs2_.size();
 SimplexId newNbCC = ripsArcs_.size();

 // std::cout<<"Leaves arc : "<<std::endl;
 // for(SimplexId iArc = 0; iArc<ripsArcs_.size(); iArc++){
 //   auto& arc = ripsArcs_[iArc];
 //      std::cout<<"   arc "<<iArc<<": "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<" type "<<arc.type<<std::endl;
 // }
 // std::cout<<"Other arc : "<<std::endl;
 // for(SimplexId iArc = 0; iArc<ripsArcs2_.size(); iArc++){
 //   auto& arc = ripsArcs2_[iArc];
 //      std::cout<<"   arc "<<iArc<<": "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<" type "<<arc.type<<std::endl;
 //      std::cout<<"                 : "<<arc.ccId<<" "<<arc.birthId-nbCC<<"-"<<arc.deathId-nbCC<<"("<<arc.birth<<"-"<<arc.death<<")"<<std::endl;
 // }

 int test=0;
 // std::cout<<"test "<<test++<<std::endl;

 const auto cmp_pers = [=](const RipsArc &a0, const RipsArc &a1){
   return (a0.ccId<a1.ccId);
 };

 TTK_PSORT(this->threadNumber_, ripsArcs_.begin(), ripsArcs_.end(), cmp_pers);

 std::vector<SimplexId> oldIdToNewId(nbCC,-1);
 for(SimplexId i = 0; i<ripsArcs_.size(); i++){
  // if(ripsArcs_[i].birthId > nbCC){
  //   std::cout<<" error "<<ripsArcs_[i].birthId<<" outta "<<nbCC<<std::endl;
  // }
  oldIdToNewId[ripsArcs_[i].birthId] = i;
  // oldIdToNewId[ripsArcs_[i].deathId] = i+newNbCC;
 }
 // std::cout<<"test "<<test++<<std::endl;
 int extra_nodes=0;
 for(SimplexId i = 0; i<ripsArcs2_.size(); i++){
  SimplexId bid = ripsArcs2_[i].birthId-nbCC;
  SimplexId did = ripsArcs2_[i].deathId-nbCC;
  if(oldIdToNewId[bid] == -1){
    oldIdToNewId[bid] = extra_nodes+2*newNbCC;
    extra_nodes++;
    // std::cout<<" new node bid "<<bid<<std::endl;
  }
  if(oldIdToNewId[did] == -1){
    oldIdToNewId[did] = extra_nodes+2*newNbCC;
    extra_nodes++;
    // std::cout<<" new node did "<<did<<std::endl;;
  }
 }
 // std::cout<<" Extra nodes : "<<extra_nodes<<std::endl;
 // std::cout<<"test "<<test++<<std::endl;

 for(SimplexId i = 0; i<ripsArcs_.size(); i++){
   // if(oldIdToNewId[ripsArcs_[i].birthId]==-1 or oldIdToNewId[ripsArcs_[i].deathId-nbCC]==-1){
     // std::cout<<"shit here"<<oldIdToNewId[ripsArcs_[i].birthId]<<" "<<oldIdToNewId[ripsArcs_[i].deathId-nbCC]<<std::endl;
   // }
  ripsArcs_[i].birthId = oldIdToNewId[ripsArcs_[i].birthId];
  ripsArcs_[i].deathId = oldIdToNewId[ripsArcs_[i].deathId-nbCC]+newNbCC;
 }
 // std::cout<<"test "<<test++<<std::endl;

 for(SimplexId i = 0; i<ripsArcs2_.size(); i++){
  ripsArcs2_[i].birthId = oldIdToNewId[ripsArcs2_[i].birthId-nbCC] + newNbCC;
  ripsArcs2_[i].deathId = oldIdToNewId[ripsArcs2_[i].deathId-nbCC] + newNbCC;
 }



 // std::cout<<" sort rips arcs"<<std::endl;
 const auto cmp = [=](const RipsArc &a0, const RipsArc &a1){
   return ((a0.death-a0.birth) > (a1.death-a1.birth)) || (((a0.death-a0.birth) == (a1.death-a1.birth)) && a0.ccId<a1.ccId);
 };

 TTK_PSORT(this->threadNumber_, ripsArcs_.begin(), ripsArcs_.end(), cmp);
 // std::cout<<" sorting done!"<<std::endl;

 std::cout<<" found "<<ripsArcs_.size()<<" leaves arcs"<<std::endl;
 std::cout<<" and found "<<ripsArcs2_.size()<<" branching arcs"<<std::endl;

 // std::cout<<"Leaves arc : "<<std::endl;
 // for(SimplexId iArc = 0; iArc<ripsArcs_.size(); iArc++){
 //   auto& arc = ripsArcs_[iArc];
 //      std::cout<<"   arc "<<iArc<<": "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<" type "<<arc.type<<std::endl;
 // }
 // std::cout<<"Other arc : "<<std::endl;
 // for(SimplexId iArc = 0; iArc<ripsArcs2_.size(); iArc++){
 //   auto& arc = ripsArcs2_[iArc];
 //      std::cout<<"   arc "<<iArc<<": "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<"("<<arc.birth<<"-"<<arc.death<<")"<<" type "<<arc.type<<std::endl;
 // }
ripsArcs_.insert(
      ripsArcs_.end(),
      std::make_move_iterator(ripsArcs2_.begin()),
      std::make_move_iterator(ripsArcs2_.end()));
 // std::cout<<" now size is "<<ripsArcs_.size()<<" full arcs"<<std::endl;

 // for(auto& arc : ripsArcs_){
 //   std::cout<<" arc "<<arc.ccId<<" "<<arc.death-arc.birth<<std::endl;
 // }
 // for(SimplexId iArc = 0; iArc< ripsArcs_.size(); iArc++){

 //       auto ccId = ripsArcs_[iArc].ccId;
 //       auto threshold = ripsArcs_[iArc].death;
 //       auto birth = ripsArcs_[iArc].birth;
 //       auto persistence = ripsArcs_[iArc].death-ripsArcs_[iArc].birth;
 //       auto type = ripsArcs_[iArc].type;

 //     std::cout<<"    arc  "<<iArc<<": "<<ccId<<" "<<birth<<"-"<<threshold<<" of size "<<sizeAtThreshold(ccId, threshold)<<std::endl;
 // }
  
 // decrease the infinite arc a bit
 // ripsArcs_[0].death = ripsArcs_[1].death*2;

 TTK_PSORT(this->threadNumber_, ripsArcs_.begin(), ripsArcs_.end(), cmp);

 std::cout<<" Found "<<ripsArcs_.size()<<" rips arcs"<<std::endl;
}


void RipsTree::buildRipsArcs(SimplexId sizeMin){
  std::cout<<" build rips arcs "<<std::endl;
    // SimplexId nbCC = this->getNumberOfPoints();


    std::vector<std::vector<RipsArc>> ripsArcsPerThread(1);

    SimplexId nbCC = this->getNumberOfFilteredComponents(1);

    SimplexId nbArc=0;
    for(SimplexId filteredId=0; filteredId<nbCC; filteredId++){

      SimplexId iCC = filteredIdToNodeId_[filteredId];
      auto& cc = this->getComponent(iCC);

// #ifdef TTK_ENABLE_OPENMP
//       const auto tid = omp_get_thread_num();
// #else
      const auto tid = 0;
// #endif // TTK_ENABLE_OPENMP

      // skip small components
      if(cc.size()<sizeMin){
        // printWrn("this should not happen: filtered cc < sizemin");
        continue;
      }

      // bool addedArc=false;

      double currentThreshold = cc.death();
      SimplexId currentDeathId = nbCC+filteredId;
      SimplexId currentSize = cc.size();

      for(SimplexId iChild : cc.getChildren()){
        auto& child = this->getComponent(iChild);

        SimplexId nextSize = currentSize-child.size();
        if(nextSize<sizeMin){
          break;
        }
        if(child.size() >= sizeMin){

          int criticalType = 2; //(currentDeathId==cc.deathId()) ? 3 : 1;
          SimplexId filteredChildId = child.filteredId();

          ripsArcsPerThread[tid].push_back({iCC, child.death(), currentThreshold, nbCC+filteredChildId, currentDeathId, criticalType});
           // std::cout<<"    arc  "<<nbArc<<": "<<iCC<<" "<<child.death()<<"-"<<currentThreshold<<" of size "<<sizeAtThreshold(iCC, currentThreshold)<<std::endl;
           nbArc++;
           currentDeathId = nbCC+filteredChildId;
        }

        currentThreshold = child.death();
        currentSize = nextSize;
      }

      // add last arc
      ripsArcsPerThread[tid].push_back({iCC, 0, currentThreshold, filteredId, currentDeathId, 0});
      // std::cout<<"    arc  "<<nbArc<<": "<<iCC<<" "<<0<<"-"<<currentThreshold<<" of size "<<sizeAtThreshold(iCC, currentThreshold)<<std::endl;
      nbArc++;
    }

  // aggregate per-threads vectors in global vector
  size_t ripsArcsSize=0;
  for(auto& vec : ripsArcsPerThread){
    ripsArcsSize+=vec.size();
  }
  ripsArcs_.clear();
  ripsArcs_.reserve(ripsArcsSize);
  for(auto& vec : ripsArcsPerThread){
    for(auto arc : vec){
      ripsArcs_.push_back(arc);
    }
  }


  const auto cmp = [=](const RipsArc &a0, const RipsArc &a1){
    return ((a0.death-a0.birth) > (a1.death-a1.birth)) || (((a0.death-a0.birth) == (a1.death-a1.birth)) && a0.ccId<a1.ccId);
  };

  TTK_PSORT(this->threadNumber_, ripsArcs_.begin(), ripsArcs_.end(), cmp);

  // for(auto& arc : ripsArcs_){
  //   std::cout<<" arc "<<arc.ccId<<" "<<arc.death-arc.birth<<std::endl;
  // }
  for(SimplexId iArc = 0; iArc< ripsArcs_.size(); iArc++){

        auto ccId = ripsArcs_[iArc].ccId;
        auto threshold = ripsArcs_[iArc].death;
        auto birth = ripsArcs_[iArc].birth;
        auto persistence = ripsArcs_[iArc].death-ripsArcs_[iArc].birth;
        auto type = ripsArcs_[iArc].type;

      std::cout<<"    arc  "<<iArc<<": "<<ccId<<" "<<birth<<"-"<<threshold<<" of size "<<sizeAtThreshold(ccId, threshold)<<std::endl;
  }
  
  // decrease the infinite arc a bit
  // ripsArcs_[0].death = ripsArcs_[1].death*2;

  std::cout<<" Found "<<ripsArcs_.size()<<" rips arcs"<<std::endl;
}

void RipsTree::buildRipsTree(DiagramType& diagram){

  // printMsg("Build Rips Tree", 0, debug::LineMode::REPLACE);
  std::cout<<"build rips tree"<<std::endl;
  SimplexId nbPair = diagram.size();

  this->clear();
  // std::cout<<"\n\n\nBuild Rips Tree with "<<nbPair<<" nodes"<<std::endl;

  components_.resize(nbPair);
  nodeIds_.resize(nbPair);
  std::iota(nodeIds_.begin(), nodeIds_.end(), 0);

  std::vector<SimplexId>& roots = tree_roots_;
  roots.clear();

  for(SimplexId iPair=0; iPair<nbPair; iPair++){
    auto& pp = diagram[iPair];
    auto& cc = components_[iPair];
    cc.copy(pp, iPair);
    cc.clearChildren();
    cc.setParent(pp.cc_parent);
    if(iPair!=pp.cc_parent){
      components_[pp.cc_parent].addChild(iPair);
    }else{
      roots.push_back(iPair);
      // std::cout<<"found root "<<iPair<<std::endl;
      // std::cout<<" size is "<<pp.cc_size<<std::endl;
      // std::cout<<" size is also "<<components_[iPair].size()<<std::endl;
    }
    // std::cout<<" done node "<<iPair<<std::endl;
  }
  // printNodeIds();
  maxDeath_=components_[0].death();
  SimplexId currentIndex=0;
  for(auto root : roots){
    currentIndex = fillIndices(root, currentIndex);
  }
  // std::cout<<"filled indices"<<std::endl;
  // buildPointIdToNodeId();
  // printNodeIds();
  printMsg("Built Rips Tree", 1);
}

SimplexId RipsTree::fillIndices(SimplexId ccId, SimplexId currentId){

  auto& cc = components_[ccId];

  cc.setStartNodeId(currentId);
  nodeIds_[currentId] = ccId;
  currentId+=1;
  for(const auto& child : cc.getChildren()){
    currentId = fillIndices(child, currentId);
  }
  cc.setEndNodeId(currentId-1);
  return currentId;
}

SimplexId RipsTree::sizeAtThreshold(SimplexId ccId, double threshold) const{

  auto& cc = components_[ccId];
  if(cc.death()<threshold){
    return 0;
  }
  SimplexId size = 1;

  for(auto child : cc.getChildren()){
    if(components_[child].death()<threshold){
      size+=components_[child].size();
    }
  }
  return size;
}
SimplexId RipsTree::augmentedSizeAtThreshold(SimplexId ccId, double threshold, const SimplexId sizeMin) const{

  auto& cc = components_[ccId];
  if(cc.death()<threshold){
    return 0;
  }
  SimplexId size = 1;

  for(auto child : cc.getChildren()){
    if(components_[child].death()<threshold or components_[child].size() < sizeMin){
      size+=components_[child].size();
    }
  }
  return size;
}

// // returns 1 if component exists at given threshold
// int RipsTree::getComponentNodeIdsAtThreshold(SimplexId ccId, double threshold, std::vector<SimplexId>& points)const{
//     auto& cc = components_[ccId];
//     points.clear();
//     if(cc.death()<threshold){
//       return 0;
//     }
//     points.push_back(ccId);
//     SimplexId start = cc.startNodeId();
//     SimplexId end = cc.endNodeId();

//     SimplexId i = start+1;
//     while(i<=end){
//       if(components_[nodeIds_[i]].death()<threshold){ //component is dead, its points must be added
//         points.push_back(nodeIds_[i]);
//         i++;
//       }else{ // component is not dead, must skip until the end of its points
//         i=components_[nodeIds_[i]].endNodeId()+1;
//       }
//     }
//     return 1;
// }

// returns 1 if component exists at given threshold
int RipsTree::getComponentNodeIdsAtThreshold(SimplexId ccId, double threshold, std::vector<SimplexId>& points, const SimplexId sizeMin)const{
  if(threshold>maxDeath_){
    threshold=maxDeath_;
  }
    points.clear();
    if(ccId>=getNumberOfPoints()){
      return -1;
    }
    auto& cc = components_[ccId];
    if(cc.death()<threshold){
      return 0;
    }
    points.push_back(ccId);
    SimplexId start = cc.startNodeId();
    SimplexId end = cc.endNodeId();

    SimplexId i = start+1;
    while(i<=end){
      if(components_[nodeIds_[i]].death()<threshold /*or fullSize(nodeIds_[i]) < sizeMin*/){ //component is dead, its points must be added
        SimplexId last = components_[nodeIds_[i]].endNodeId();
        for(int j=i; j<=last; j++){
          points.push_back(nodeIds_[j]);
        }
        i=last+1;
      }else{ // component is not dead, must skip until the end of its points
        i=components_[nodeIds_[i]].endNodeId()+1;
      }
    }
    return 1;
}
// returns 1 if component exists at given threshold
int RipsTree::getComponentPointIdsAtThreshold(SimplexId ccId, double threshold, std::vector<SimplexId>& points, const SimplexId sizeMin)const{
  if(threshold>maxDeath_){
    threshold=maxDeath_;
  }
    points.clear();
    if(ccId>=getNumberOfPoints()){
      return -1;
    }
    auto& cc = components_[ccId];
    if(cc.death()<threshold){
      return 0;
    }
    points.push_back(getPointId(ccId));
    SimplexId start = cc.startNodeId();
    SimplexId end = cc.endNodeId();

    SimplexId i = start+1;
    while(i<=end){
      if(components_[nodeIds_[i]].death()<threshold /*or fullSize(nodeIds_[i]) < sizeMin*/){ //component is dead, its points must be added
        SimplexId last = components_[nodeIds_[i]].endNodeId();
        for(int j=i; j<=last; j++){
          points.push_back(getPointId(nodeIds_[j]));
        }
        i=last+1;
      }else{ // component is not dead, must skip until the end of its points
        i=components_[nodeIds_[i]].endNodeId()+1;
      }
    }
    return 1;
}

// returns 1 if component exists at given threshold
void RipsTree::getFullComponentNodeIds(SimplexId ccId, std::vector<SimplexId>& points) const{

    auto& cc = components_[ccId];
    SimplexId start = cc.startNodeId();
    SimplexId end = cc.endNodeId();
    points.resize(end-start+1);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i=0; i<points.size(); i++){
      points[i] = nodeIds_[start+i];
    }
}

void RipsTree::printNodeIds() const{

  std::stringstream msg;
  std::cout<<"Node Ids: ("<<nodeIds_.size()<<" nodes)";
  printMsg(msg.str());

  for(auto id : nodeIds_){
    std::cout<<" "<<id;
  }
  std::cout<<std::endl;

  for(auto i=0; i<components_.size(); i++){
    auto& cc = components_[i];
    std::cout<<"component "<<i<<": "<<cc.string()<<std::endl;
  }
}

void RipsTree::printSizes(double threshold) const{
  std::stringstream msg;
  msg<<"Node sizes at thresold: "<<threshold;
  printMsg(msg.str());
  for(auto i=0; i<components_.size(); i++){
    // auto& cc = components_[i];
    std::cout<<"cc "<< i <<": "<<sizeAtThreshold(i, threshold)<<std::endl;
  }

}

// sorted vectors as input
// double RipsTree::computeJaccardIndex(const std::vector<SimplexId>& cc1, const std::vector<SimplexId>& cc2) const{
//     std::vector<SimplexId> interVec;
//     std::vector<SimplexId> unionVec;
//     std::set_intersection(cc1.cbegin(), cc1.cend(), cc2.cbegin(), cc2.cend(), std::back_inserter(unionVec));
//     std::set_union(cc1.cbegin(), cc1.cend(), cc2.cbegin(), cc2.cend(), std::back_inserter(unionVec));

//     return (double)interVec.size()/(double)unionVec.size();
// }

// double RipsTree::computeJaccardIndex(SimplexId cc1, SimplexId cc2) const{
//   std::vector<SimplexId> points1, points2;
//   getFullComponentPoints(cc1, points1);
//   getFullComponentPoints(cc2, points2);
//   double jaccard=computeJaccardIndex(points1, points2);
//   return jaccard;
// }

void RipsTree::findContainingComponents(std::vector<SimplexId>& pointsIds, double threshold, std::vector<SimplexId>& ccIds) const{

  ccIds.clear();
  for(SimplexId p : pointsIds){
    ccIds.push_back(pointIdToNodeId(p, threshold));
  }
  std::sort(ccIds.begin(), ccIds.end());
  auto last=std::unique(ccIds.begin(), ccIds.end());
  ccIds.erase(last, ccIds.end());
}

void RipsTree::buildPointIdToNodeId(){
  if(pointIdToNodeId_.size()==getNumberOfPoints()){
    return;
  }
  pointIdToNodeId_.resize(getNumberOfPoints());
  for(SimplexId i=0; i<getNumberOfPoints(); i++){
    pointIdToNodeId_[components_[i].birthId()] = i;
  }
  return;
}

SimplexId RipsTree::pointIdToNodeId(SimplexId pointId, double threshold)const{
  // if(pointIdToNodeId_.empty()){
  //   buildPointIdToNodeId();
  // }
  if(threshold>maxDeath_){
    threshold=maxDeath_;
  }
  SimplexId ccId = pointIdToNodeId_[pointId];
  while(components_[ccId].death()<threshold){
    ccId = components_[ccId].parent();
  }
  return ccId;
}

SimplexId RipsTree::getNumberOfFilteredComponents(SimplexId sizeMin){

  std::cout<<"FILTERING with sizemin "<<sizeMin<<std::endl;
  if(1 /*sizeMin!=sizeMin_ or filteredNumberOfComponents_==-1*/){
    sizeMin_=sizeMin;
    filteredIdToNodeId_.clear();
    SimplexId nbCC = this->getNumberOfPoints();
    for(SimplexId iCC = 0; iCC<nbCC; iCC++){
      auto& cc = this->getComponent(iCC);
      // std::cout<<" cc "<<iCC<<" of size "<<cc.size()<<std::endl;
      // skip small components
      if(cc.size()>=sizeMin){
        cc.setFilteredId(filteredIdToNodeId_.size());
        filteredIdToNodeId_.push_back(iCC);
      }
    }
    filteredNumberOfComponents_=filteredIdToNodeId_.size();
  }
  // std::cout<<" final filtering: "<<std::endl;
  // for(SimplexId id : filteredIdToNodeId_){
  //   std::cout<<" "<<id;
  // }
  // std::cout<<"\n and cc have filtered indices: "<<std::endl;
  // for(SimplexId iCC = 0; iCC<getNumberOfPoints(); iCC++){
  //     auto& cc = this->getComponent(iCC);
  //     std::cout<<" ,  "<<iCC<<" : "<<cc.filteredId();
  // }
  // std::cout<<" "<<std::endl;
  return filteredNumberOfComponents_;
}

// void RipsTree::getClusteringAtThreshold(double threshold, SimplexId sizeMin, std::vector<SimplexId>& sampling_indexes){
//   SimplexId nbPoints = sampling_indexes.size();

//   std::vector<SimplexId> clustering(nbPoints,-1);

//   for(SimplexId iCC=0; iCC<nbPoints; iCC++){
//     auto& c = components_[iCC];

//     if(c.death()<threshold){
//       break;
//     }

//     std::vector<SimplexId> children;
//     getComponentNodeIdsAtThreshold(iCC, threshold, children);
//     if(children.size()>=sizeMin){
//       for(auto& iChild : children){
//         clustering[iChild] = iCC;
//       }
//     }
//   }
// }
//

void RipsTree::computeArcAutoSimilarity(const std::vector<std::vector<double>>& distMat, const std::vector<SimplexId>& sampling){

  SimplexId nbArcs = ripsArcs_.size();
  // std::vector<double> autosimilarity(nbArcs);

  for(SimplexId iArc=0; iArc<nbArcs; iArc++){
    auto& arc = ripsArcs_[iArc];

    std::vector<SimplexId> pointIds;
    getComponentPointIdsAtThreshold(arc.ccId, arc.death, pointIds,1);
    SimplexId nPoints = pointIds.size();
    // std::cout<<"iArc "<<iArc<<" and size "<<nPoints<<std::endl;

    double sum_dist=0;
    for(SimplexId i=0; i<nPoints; i++){
      for(SimplexId j=i+1; j<nPoints; j++){
        // std::cout<<i<<" "<<j<<" "<<pointIds[i]<<" "<<pointIds[j]<<" "<<sampling[pointIds[i]]<<" "<<sampling[pointIds[j]]<<" "<<distMat[sampling[pointIds[i]]][sampling[pointIds[j]]]<<std::endl;
        sum_dist += distMat[sampling[pointIds[i]]][sampling[pointIds[j]]];
      }
    }
    // std::cout<<"found autoSim "<<sum_dist<<std::endl;
    arc.autosimilarity = sum_dist/nPoints;

  }
  // std::sort(ripsArcs_.begin(), ripsArcs_.end(), 
      // [=](const auto& i, const auto& j) { return i.autosimilarity< j.autosimilarity; } );

  // std::sort(autosimilarity.begin(), autosimilarity.end());
}

void ttk::compute5DAABBMultiscale(const std::vector<std::vector<SimplexId>>& local_ids, const std::vector<std::vector<double>>& fitCoordsVec, const int dim, std::vector<double>& outputAABB){
  outputAABB.resize(dim*2);

  outputAABB[0]=std::numeric_limits<double>::min();
  for(int i=0; i<dim; i++){
    outputAABB[2*i]=std::numeric_limits<double>::max();
    outputAABB[2*i+1]=std::numeric_limits<double>::lowest();
  }

  int nbScales = local_ids.size();
  if(fitCoordsVec.size()!=nbScales){
    std::cout<<"ERROR SIZE MISMATCH 5d aabb computation"<<std::endl;
    return;
  }

  for(int iScale=0; iScale<nbScales; iScale++){
    auto& fitCoords = fitCoordsVec[iScale];

    for(auto& localId : local_ids[iScale]){
      for(int i=0; i<dim; i++){
        double v_i = fitCoords[localId*dim + i];
        if(v_i <  outputAABB[2*i]){
          outputAABB[2*i]=v_i;
        }
        if(v_i >  outputAABB[2*i+1]){
          outputAABB[2*i+1]=v_i;
        }
      }
    }
  }
}
void ttk::compute5DAABB(const std::vector<SimplexId>& local_ids, const std::vector<double>& fitCoords, const int dim, std::vector<double>& outputAABB){
  outputAABB.resize(dim*2);

  for(int i=0; i<dim; i++){
    outputAABB[2*i]=std::numeric_limits<double>::max();
    outputAABB[2*i+1]=std::numeric_limits<double>::lowest();
  }

  for(auto& localId : local_ids){
    for(int i=0; i<dim; i++){
      double v_i = fitCoords[localId*dim + i];
      if(v_i <  outputAABB[2*i]){
        outputAABB[2*i]=v_i;
      }
      if(v_i >  outputAABB[2*i+1]){
        outputAABB[2*i+1]=v_i;
      }
    }
  }
}

double ttk::AABBToDiag(std::vector<double>& aabb, int dim){
    if(aabb.size()!=2*dim){
      std::cout<<"ERROR wrong dims aabb"<<std::endl;
    }
  double diag = 0;
  for(int i=0; i<dim; i++){
    double l = aabb[2*i+1]-aabb[2*i];
    diag+= l*l;
  }
  return std::sqrt(diag);
}

double ttk::compute5DAABBDiagonalMultiscale(const std::vector<std::vector<SimplexId>>& local_ids, const std::vector<std::vector<double>>& fitCoordsVec, const int dim){
  std::vector<double> aabb;
  compute5DAABBMultiscale(local_ids, fitCoordsVec, dim, aabb);
  return AABBToDiag(aabb, dim);
}

double ttk::compute5DAABBDiagonal(const std::vector<SimplexId>& local_ids, const std::vector<double>& fitCoords, const int dim){
  std::vector<double> aabb;
  compute5DAABB(local_ids, fitCoords, dim, aabb);
  return AABBToDiag(aabb, dim);
}

int ttk::findScaleId(SimplexId ripsBirthId, const std::vector<SimplexId>& pointsAtScaleStartIdx){
  for(int iScale = pointsAtScaleStartIdx.size()-1; iScale>=0; iScale--){
    if(ripsBirthId>=pointsAtScaleStartIdx[iScale]){
      return iScale;
    }
  }
  std::cout<<"ERROR Found no scale id from pointId"<<std::endl;
  return -1;
}

double ttk::compute5DAABBDiagonalMultiscale(const std::vector<SimplexId>& multiscale_ids,
                                            const std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                            const std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                            const std::vector<SimplexId>& pointsAtScaleStartIdx,
                                            const std::vector<std::vector<double>>& fitCoordsVec,
                                            const int dim){

  int nbScales = pointsAtScaleStartIdx.size();
  std::vector<std::vector<SimplexId>> local_ids(nbScales);
  for(auto id : multiscale_ids){
    int iScale = findScaleId(id, pointsAtScaleStartIdx);
    SimplexId id_mono_eps = id - pointsAtScaleStartIdx[iScale];
    local_ids[iScale].push_back(filteredToLocalId[iScale][epsilonIdToFilteredId[iScale][id_mono_eps]]);
  }

  std::vector<double> aabb;
  compute5DAABBMultiscale(local_ids, fitCoordsVec, dim, aabb);
  return AABBToDiag(aabb, dim);
}
