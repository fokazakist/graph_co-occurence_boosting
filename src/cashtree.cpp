#include "gspan.h"
#include <sstream>
#include <set>
#include <strstream>

void Ctree::print(){
  std::cout << g2tracers.size()<<":"<<pattern <<std::endl;
  if(children.size()!=0){
    for(list<Ctree*>::iterator it=children.begin();it!=children.end();++it){
      for(unsigned int i=1;i<pattern.size();i++){
	std::cout <<" ";
      }
      std::cout << "-";
      (*it)->print();
    }
  }
}
void CRoot::print(){
    for(list<Ctree*>::iterator it=one_edge_graphs.begin();it!=one_edge_graphs.end();++it){
      (*it)->print();
    }
}

void Tdelete(Ctree* tree){
  if(tree->children.size()!=0){
    for(list<Ctree*>::iterator it=tree->children.begin();it!=tree->children.end();++it){
      Tdelete(*it);
    }
  }
  delete tree;
}

void Gspan::Crun(){
  if(first_flag==true){
    first_tree_make();
    first_flag=false;
    std::cout << TNnum <<std::endl;
    return;
  }

  DFSCode  dcode;
  for(list<Ctree*>::iterator it=croot->one_edge_graphs.begin();it!=croot->one_edge_graphs.end();++it){
    dcode = (*it)->pattern[(*it)->pattern.size()-1];
    cash_tsearch(croot->heap[dcode.labels],**it);
  }
  std::cout << TNnum <<std::endl;
}

void Gspan::cash_tsearch(GraphToTracers& g2tracers,Ctree& node){

  double gain=0.0;
  double upos=0.0;
  double uneg=0.0;
  gain=-wbias;
  upos=-wbias;
  uneg=wbias;

  for(GraphToTracers::iterator it=g2tracers.begin();it!=g2tracers.end();++it){
    int gid = it->first;
    gain += 2 * corlab[gid] * weight[gid];
    if(corlab[gid]>0){
      upos += 2 * weight[gid];
	}else{
      uneg += 2 * weight[gid];
    }
  }

  if(fabs(opt_pat.gain)-std::max(upos,uneg)>=-1e-10) {
    return;
  }

  double gain_abs = fabs(gain);
  if(gain_abs > fabs(opt_pat.gain) || (fabs(gain_abs - fabs(opt_pat.gain)) < 1e-10 && node.pattern.size() < opt_pat.size)){
    opt_pat.gain = gain;
    opt_pat.size = node.pattern.size();
    opt_pat.locsup.clear();
    for(GraphToTracers::iterator it=g2tracers.begin();it!=g2tracers.end();++it){
      opt_pat.locsup.push_back(it->first);
    }
    std::ostrstream ostrs;
    ostrs <<node.pattern;
    ostrs << std::ends;
    opt_pat.dfscode = ostrs.str();
  }
  if(node.children.size()==0){
    edge_grow(node);
    return;
  }
  DFSCode  dcode;
  Pair pkey;
  for(list<Ctree*>::iterator it=node.children.begin();it!=node.children.end();++it){
    dcode = (*it)->pattern[(*it)->pattern.size()-1];
    if(dcode.labels.z == -1){
      pkey.set(dcode.time.b,dcode.labels.y);
      cash_tsearch(node.b_heap[pkey],**it);
    }else{
      pkey.set(dcode.labels.y,dcode.labels.z);
      cash_tsearch(node.f_heap[dcode.time.a][pkey],**it);
    }
  }
}

void Gspan::first_tree_make(){
  /***   init CRoot         ***/
  croot = new CRoot;
  croot->one_edge_graphs.resize(0);

  /****  construct CRoot   ****/
  map<Triplet,GraphToTracers>& heap = croot->heap;
    p_count = 0;
  for(unsigned int gid = 0; gid < gdata.size(); ++gid){
    EdgeTracer cursor;
    Triplet wild_edge;
    Graph& g = gdata[gid];

    for(unsigned int v=0; v<g.size(); ++v){
      for(vector<Edge>::iterator e = g[v].begin(); e != g[v].end(); ++e){

	if (e->labels.x <= e->labels.z){
	  cursor.set(v,e->to,e->id,0);
	  heap[e->labels][gid].push_back(cursor);

	  if(wildcard_r>0){
	    wild_edge = e->labels;
	    wild_edge.z =999;
	    heap[wild_edge][gid].push_back(cursor);
	    wild_edge = e->labels.reverse();
	    wild_edge.z =999;
	    cursor.set(e->to,v,e->id,0);
	    heap[wild_edge][gid].push_back(cursor);
	  }
	}
      }
    }
  }
  pattern.resize(1);
  for(map<Triplet,GraphToTracers>::iterator it = croot->heap.begin(); it != croot->heap.end(); ++it){		
    pattern[0].labels = it->first;
    pattern[0].time.set(0,1);
    edge_grow(it->second,*croot);
    pattern.resize(1);
  }
  std::cout << p_count << std::endl;
}

void Gspan::coocsearch(){
  //this function is start position to search cooc pattern
  

}
