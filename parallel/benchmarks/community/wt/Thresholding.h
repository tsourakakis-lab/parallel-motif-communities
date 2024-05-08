// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#pragma once

#include <algorithm>
#include <fstream>

#include "gbbs/gbbs.h"

#include "benchmarks/Connectivity/WorkEfficientSDB14/Connectivity.h"
#include "benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h"

namespace gbbs
{

  // using Clustering = std::vector<std::vector<uintE>>;

  double ComputeModularity(
    std::vector<std::vector<gbbs::uintE>>& initial_clustering,
    gbbs::symmetric_graph<gbbs::symmetric_vertex, gbbs::empty>& graph, parlay::sequence<gbbs::uintE>& cluster_ids,
    double resolution){
    double total_edge_weight = 0;
    double modularity = 0;
    for (std::size_t i = 0; i < graph.n; i++) {
      auto vtx = graph.get_vertex(i);
      auto nbhrs = vtx.out_neighbors();
      double deg_i = vtx.out_degree();
      for (std::size_t j = 0; j < deg_i; j++) {
        total_edge_weight++;
        auto nbhr = nbhrs.get_neighbor(j);
        // auto nbhr = std::get<0>(nbhrs[static_cast<int>(j)]);
        if (cluster_ids[i] == cluster_ids[nbhr]) {
          modularity++;
        }
      }
    }
    for (std::size_t i = 0; i < initial_clustering.size(); i++) {
      double degree = 0;
      for (std::size_t j = 0; j < initial_clustering[i].size(); j++) {
        auto vtx_id = initial_clustering[i][j];
        auto vtx = graph.get_vertex(vtx_id);
        degree += vtx.out_degree();
      }
      modularity -= (resolution * degree * degree) / (total_edge_weight);
    }
    modularity = modularity / (total_edge_weight);
    return modularity;
  }

  double ComputeConductance(
    std::vector<std::vector<gbbs::uintE>>& initial_clustering,
    gbbs::symmetric_graph<gbbs::symmetric_vertex, gbbs::empty>& graph, 
    parlay::sequence<gbbs::uintE>& cluster_ids){
    double conductance = 0;
    for (std::size_t i = 0; i < initial_clustering.size(); i++) {
      double cut_size = 0;
      double vol = 0;
      for (std::size_t j = 0; j < initial_clustering[i].size(); j++) {
        auto vtx_id = initial_clustering[i][j];
        auto vtx = graph.get_vertex(vtx_id);
        double deg_i = vtx.out_degree();
        vol += deg_i;
        auto nbhrs = vtx.out_neighbors();
        for (std::size_t j = 0; j < deg_i; j++) {
          auto nbhr = nbhrs.get_neighbor(j);
          if (cluster_ids[vtx_id] != cluster_ids[nbhr]) {
            cut_size++;
          }
        }
      }
      conductance += cut_size/vol;
    }
    return conductance/initial_clustering.size();
  }

  void split(const std::string &s, char delim, std::vector<uintE> &elems)
  {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim))
    {
      elems.push_back(std::stoi(item));
    }
  }

  void ReadCommunities(const char* filename,
    std::vector<std::vector<gbbs::uintE>>& communities, std::vector<std::vector<int>>& node_community_map,
    std::size_t n) {
    std::ifstream infile(filename);
    for (gbbs::uintE i=0; i<n; i++){
      std::vector<int> node_community_belonging;
      node_community_map.push_back(node_community_belonging);
    }
    // if (!infile.is_open()) {
    //   return absl::NotFoundError("Unable to open file.");
    // }
    std::string line;
    int community_idx = 0;
    while (std::getline(infile, line)) {
      std::vector<gbbs::uintE> row_values;
      split(line, '\t', row_values);
      if (row_values.size()<3)
        continue;
      std::sort(row_values.begin(), row_values.end());
      communities.push_back(row_values);
      for (auto node : row_values) {
        node_community_map[static_cast<int>(node)].push_back(community_idx);
      }
      community_idx++;
    }
  }

  std::unordered_map<int, int> counter(const std::vector<int>& vec) {
    std::unordered_map<int, int> counts;
    for (int num : vec) {
        counts[num]++;
    }
    return counts;
  }

  int CompareCommunities(const char* filename, std::vector<std::vector<uintE>>& clustering, std::size_t n) {
    std::cout << "community init" << std::endl;
  std::vector<std::vector<gbbs::uintE>> communities;
  std::vector<std::vector<int>> node_community_map;
  ReadCommunities(filename, communities, node_community_map, n);
  auto precision_vec = sequence<double>::from_function(clustering.size(), [&](size_t i) { return 0.0; });
  auto recall_vec = sequence<double>::from_function(clustering.size(), [&](size_t i) { return 0.0; });
  auto valid_size_vec = sequence<int>::from_function(clustering.size(), [&](size_t i) { return 0; });
    std::cout << "read community" << std::endl;
  // gbbs::sequence<double> precision_vec(clustering.size(), [](std::size_t i){return 0.;});
  // gbbs::sequence<double> recall_vec(clustering.size(), [](std::size_t i){return 0.;});
  // gbbs::sequence<int> valid_size_vec(clustering.size(), [](std::size_t i){return 0;});
  parallel_for(0, clustering.size(), [&](std::size_t i) {
    auto cluster = clustering[i];
    std::sort(cluster.begin(), cluster.end());
  });

  parallel_for(0, clustering.size(), [&](int j) {
    auto partition = clustering[j];
    std::vector<int> candidates;
    int valid_size = 0;
    for (auto v : partition) {
      if (node_community_map[v].size()>0) valid_size++;
      for (auto community_of_v:node_community_map[v])
        candidates.push_back(community_of_v);

    }
    std::unordered_map<int, int> candidate_count = counter(candidates);
    double best_jac = 0.0;
    int intersect_size = 0;
    int best_candidate=0;
    double best_gt_size=0.0, cluster_size=0.0;
    for (const auto& pair : candidate_count) {
        int candidate = pair.first;
        int cur_count = pair.second;
        int gt_size = communities[candidate].size();
        int cluster_size = partition.size();
        if ((double)(cur_count)/(gt_size+cluster_size-cur_count)>best_jac){
          best_jac = (double)(cur_count)/(gt_size+cluster_size-cur_count);
          best_candidate = candidate;
          intersect_size = cur_count;
        }
    }
    if (intersect_size!=0){
        best_gt_size = (double) communities[best_candidate].size();
        cluster_size = (double) partition.size();
        precision_vec[j] = (double) intersect_size / cluster_size * valid_size;
        recall_vec[j] = (double) intersect_size / best_gt_size * valid_size;
        valid_size_vec[j] = valid_size;
      }
  });

  double avg_precision = parlay::reduce(precision_vec);
  double avg_recall = parlay::reduce(recall_vec);
  int valid_size_sum = parlay::reduce(valid_size_vec);
  avg_precision /= valid_size_sum;
  avg_recall /= valid_size_sum;
  std::cout << "Avg precision recall: [" << std::setprecision(17) << avg_precision << "," << std::setprecision(17) << avg_recall << "," << clustering.size() << "]" << std::endl;
  return 0;
  }

  // int CompareCommunities(const char *filename, Clustering clustering)
  // {
  //   // std::cout << filename << "\n";
  //   std::vector<std::vector<uintE>> communities;
  //   ReadCommunities(filename, communities);

  //   auto precision_vec = sequence<double>::from_function(communities.size(), [&](size_t i) { return 0.0; });
  //   auto recall_vec = sequence<double>::from_function(communities.size(), [&](size_t i) { return 0.0; });
  //   auto f1_vec = sequence<double>::from_function(communities.size(), [&](size_t i) { return 0.0; });
  //   parallel_for(0, clustering.size(), [&](std::size_t i)
  //                {
  //   auto cluster = clustering[i];
  //   std::sort(cluster.begin(), cluster.end()); });
  //   parallel_for(0, communities.size(), [&](std::size_t j)
  //                {
  //   auto community = communities[j];
  //   std::vector<uintE> intersect(community.size());
  //   std::size_t max_intersect = 0;
  //   std::size_t max_idx = 0;

  //   for (std::size_t i = 0; i < clustering.size(); i++) {
  //     auto cluster = clustering[i];
  //     auto it = std::set_intersection(cluster.begin(), cluster.end(), 
  //       community.begin(), community.end(), intersect.begin());
  //     std::size_t it_size = it - intersect.begin();
  //     if (it_size > max_intersect) {
  //       max_intersect = it_size;
  //       max_idx = i;
  //     }
  //   }
  //   precision_vec[j] = (double) max_intersect / (double) clustering[max_idx].size();
  //   recall_vec[j] = (communities[j].size() == 0) ? 0 : 
  //     (double) max_intersect / (double) communities[j].size(); 
  //   f1_vec[j] = 2.0 * precision_vec[j] * recall_vec[j] / (precision_vec[j] + recall_vec[j]);
  //     });
  //   double avg_precision = parlay::reduce(precision_vec);
  //   double avg_recall = parlay::reduce(recall_vec);
  //   avg_precision /= communities.size();
  //   avg_recall /= communities.size();
  //   double avg_f1 = parlay::reduce(f1_vec);
  //   avg_f1 /= communities.size();
  //   std::cout << "Avg precision: " << avg_precision << std::endl;
  //   std::cout << "Avg recall: " << avg_recall << std::endl;
  //   std::cout << "Avg F1 score: " << avg_f1 << std::endl;
  //   return 0;
  // }

  template <class Graph, class F>
  inline void Thresholding(Graph &G, const F &f, commandLine &P)
  {
    std::string community_file = P.getOptionValue("-community", "");
    double threshold = std::stod(P.getOptionValue("-threshold", "8.0"));
    std::string method = P.getOptionValue("-method", "WT");
    using W = typename Graph::weight_type;
    timer t;
    t.start();
    auto pack_predicate = [&](const uintE &u, const uintE &v, const W &wgh)
    {
      auto their_neighbors = G.get_vertex(v).out_neighbors();
      auto triangle_count = static_cast<int>(G.get_vertex(u).out_neighbors().intersect_f_par(&their_neighbors, f));
      if (method.compare("TW")==0){
        // std::cout << "run WT" << std::endl;
        auto wedge_count = static_cast<int>(G.get_vertex(u).out_degree() + G.get_vertex(v).out_degree()) - 2 * triangle_count - 2;
        return triangle_count - wedge_count > threshold;
      } else if (method.compare("WTnorm")==0){
        // std::cout << "run WT" << std::endl;
        auto wedge_count = static_cast<int>(G.get_vertex(u).out_degree() + G.get_vertex(v).out_degree()) - 2 * triangle_count - 2;
        return (wedge_count - triangle_count)/static_cast<double>(G.get_vertex(u).out_degree() + G.get_vertex(v).out_degree()) < threshold;
      }
      else if (method.compare("TECTONIC")==0)
      {
        return triangle_count/static_cast<double>(G.get_vertex(u).out_degree() + G.get_vertex(v).out_degree()) > threshold;
      }
      else if (method.compare("JACCARD")==0)
      {
        return triangle_count/(static_cast<double>(G.get_vertex(u).out_degree() + G.get_vertex(v).out_degree())-triangle_count) > threshold;
      } 
      else if (method.compare("K3")==0)
      {
        return triangle_count > threshold;
      }
      else{
        return false;
      } 
    };
    auto DG = filterGraph(G, pack_predicate);
    double tt = t.stop();
    if (community_file.size()==0){
      std::cout << "### sparse running Time: " << tt << std::endl;
      auto K3count1 = Triangle_degree_ordering<Graph, F>(G, f);
      auto K3count2 = Triangle_degree_ordering<Graph, F>(DG, f);
      
      std::cout << "#num of edges left: " << float(DG.m)/G.m << " triangle " << float(K3count2)/K3count1 << std::endl;

    }
    timer t2;
    t2.start();
    auto components =
        workefficient_cc::CC(DG, 0.2, P.getOption("-pack"), P.getOption("-permute"));
    contract::RelabelIds(components);
    double tt2 = t2.stop();
    if (community_file.size()==0){
      std::cout << "### CC Running Time: " << tt2 << std::endl;
    }
    auto cc_f = [&](size_t i)
    { return components[i]; };
    auto cc_im = parlay::delayed_seq<uintE>(DG.n, cc_f);
    size_t numCC = workefficient_cc::num_cc(cc_im);
    std::cout << "#num of CC normalized: " << float(numCC)/G.n  << std::endl;
    std::vector<std::vector<uintE>> cluster(numCC);
    // std::vector<gbbs::uintE> cluster_ids(DG.n);
    for(gbbs::uintE  i=0;i<components.size(); i++) {
      cluster[components[i]].push_back(i);
    }

    size_t largest_cc = 0;
    for (std::size_t i = 0; i < cluster.size(); i++) {
      if (largest_cc<cluster[i].size())
        largest_cc = cluster[i].size();
    }
    std::cout << "size of the largest CC normalized: " << float(largest_cc)/G.n  << std::endl;
    
    double modularity = ComputeModularity(cluster, G, components, 1.0);
    double avg_conductance = ComputeConductance(cluster, G, components);
    std::cout << "finished partition with modularity " << modularity << ", conductance " << avg_conductance << std::endl;
    if (community_file.size()>0){
          std::cout << "threshold: " << threshold << std::endl;
          CompareCommunities(community_file.c_str(), cluster, DG.n);
        }
    return;
  }

} // namespace gbbs
