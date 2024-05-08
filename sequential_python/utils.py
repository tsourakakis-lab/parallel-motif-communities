import networkx as nx
import pandas as pd
import numpy as np
# import graph_tool.all as gt
import pyintergraph 
import matplotlib.pyplot as plt
import time
import pickle
from itertools import combinations
from collections import defaultdict
import networkx.algorithms.community as nx_comm

from collections import Counter

def evaluate2(relabel_G, re_G, communities, node_community_map, verbose=True):
    n = relabel_G.number_of_nodes()
    total_precision = 0.0
    total_recall = 0.0
    total_f1 = 0
    ttl_valid_size = 0
    if type(re_G)==nx.Graph:
        re_G = nx.connected_components(re_G)
    for cluster in re_G:
        candidate_gt = []
        valid_size = 0
        for i in cluster:
            candidate_gt += list(node_community_map[i])
            if len(node_community_map[i])!=0:
                valid_size += 1
        if len(candidate_gt)==0:
            continue
        counter = Counter(candidate_gt)
        best_gt = max(counter, key=lambda x:counter[x]/(len(cluster)+len(communities[x])-counter[x]))
        intersect_size = counter[best_gt]
        gt_size = len(communities[best_gt])
        precision = intersect_size / float(len(cluster))
        recall = intersect_size / float(gt_size)
        total_precision += precision * valid_size
        total_recall += recall * valid_size
        total_f1 += (2*precision*recall/(precision+recall)) * valid_size
        ttl_valid_size += valid_size
    return total_precision / ttl_valid_size, total_recall / ttl_valid_size, total_f1/ttl_valid_size

def thresholding_by_partition(relabel_G, thetas, edge_func, communities=None, node_community_map=None, additional_metrics=True):
    precisions = []
    recalls = []
    f1s = []
    total_t = []
    largest_size = []
    cond, num_comp, edges_remove, modulars = [], [], [], []
    percent_k3, percent_wedge = [], []
    metrics = {"precision": precisions, "recall":recalls,
              "f1":f1s, "time":total_t, "conductance":cond,
              "num_components":num_comp, "edges_removed": edges_remove,
              "modularity":modulars, "largest_size":largest_size,
              "num_k3":percent_k3, "num_wedge":percent_wedge}
#     ori_num_k3 = sum(nx.triangles(relabel_G).values()) / 3
#     ori_num_wedge = 0
#     for e in relabel_G.edges():
#         ori_num_wedge += relabel_G[e[0]][e[1]]['wedge']
    
    for theta in thetas:
        start_time = time.time()
        re_G = nx.Graph()
        re_G.add_nodes_from([i for i in range(relabel_G.number_of_nodes())])
        c = 0.0
        for e in relabel_G.edges():
            if edge_func(relabel_G, e, theta):
                re_G.add_edge(e[0], e[1])
                c += 1
        components = list(nx.connected_components(re_G))
        total_t.append(time.time() - start_time)
        print('theta=', theta)
        if communities is not None:
            if node_community_map is not None:
                p,r,f1 = evaluate2(relabel_G, components, communities, node_community_map, False)
            else:
                p,r,f1 = evaluate1(relabel_G, components, communities, False)
            precisions.append(p)
            recalls.append(r)
            f1s.append(f1)
            largest_size.append(max([len(c) for c in components]))
            modulars.append(nx_comm.modularity(relabel_G,components))
            print('precision, recall and F1:', round(p, 2), round(r, 2), round(f1, 2))
        if additional_metrics:
            num_comp.append(len(components))
#             percent_k3.append(sum(nx.triangles(re_G).values()) / 3 / ori_num_k3)
#             num_wedge = 0
#             for e in re_G.edges():
#                 num_wedge += len(set(re_G[e[0]])|set(re_G[e[1]]))
#             num_wedge = num_wedge-sum(nx.triangles(re_G).values())
#             percent_wedge.append(num_wedge/ori_num_wedge)
            edges_remove.append(c/relabel_G.number_of_edges())
#             tmp_cond, tmp_cond_big = [], []
#             for comp in components:
#                 tmp_cond.append(float(nx.cut_size(relabel_G, comp))/nx.volume(relabel_G, comp))
#                 if len(comp)>2:
#                     tmp_cond_big.append(tmp_cond[-1])
#             cond.append(np.average(tmp_cond))
#             big_comp_cond.append(np.average(tmp_cond_big))
#             print("modularity:", modulars[-1], "Num of components:",len(components), "percent of k3:", round(percent_k3[-1],4), "percent of edges:", round(c/relabel_G.number_of_edges(),4), )
    print('avg running time', np.average(total_t))
    return metrics

def evaluate1(relabel_G, re_G, communities, verbose=True):
    n = relabel_G.number_of_nodes()
    cluster_no = [-1] * n
    cluster_size = []
    if type(re_G)==nx.Graph:
        re_G = nx.connected_components(re_G)
    for cluster in re_G:
    #     cluster = map(int, line.split())
        for i in cluster:
            cluster_no[i] = len(cluster_size)
        cluster_size.append(len(cluster))

    inside = [0] * len(cluster_size)

    total_precision = 0.0
    total_recall = 0.0
    total_f1 = 0
    exact = 0
    cnt = 0

    for ty in communities:
        match = -1
        community = communities[ty]
        for i in community:
            if cluster_no[i] == -1:
                continue
            inside[cluster_no[i]] += 1
        hits = 0
        csize = 1
        for i in community:
            if cluster_no[i] == -1:
                continue
            if inside[cluster_no[i]] > hits:
                hits = inside[cluster_no[i]]
                csize = cluster_size[cluster_no[i]]
                match = cluster_no[i]
        if verbose:
            print(inside[match], cluster_size[match])
        for i in community:
            if cluster_no[i] == -1:
                continue
            inside[cluster_no[i]] = 0
        precision = hits / float(csize)
        recall = hits / float(len(community))
        total_precision += precision
        total_recall += recall
        total_f1 += (2*precision*recall/(precision+recall))
        cnt += 1
        if precision == 1.0 and recall == 1.0:
            exact += 1
    return total_precision / float(cnt), total_recall / float(cnt), total_f1/float(cnt)

def thresholding(relabel_G, thetas, edge_func, communities=None, additional_metrics=True):
    precisions = []
    recalls = []
    f1s = []
    total_t = []
    big_comp_cond = []
    cond, num_comp, edges_remove, modulars = [], [], [], []
    percent_k3, percent_wedge = [], []
    metrics = {"precision": precisions, "recall":recalls,
              "f1":f1s, "time":total_t, "conductance":cond,
              "num_components":num_comp, "edges_removed": edges_remove,
              "modularity":modulars, "big_comp_cond":big_comp_cond,
              "num_k3":percent_k3, "num_wedge":percent_wedge}
    ori_num_k3 = sum(nx.triangles(relabel_G).values()) / 3
#     ori_num_wedge = 0
#     for e in relabel_G.edges():
#         ori_num_wedge += relabel_G[e[0]][e[1]]['wedge']
    
    for theta in thetas:
        start_time = time.time()
        re_G = nx.Graph()
        re_G.add_nodes_from([i for i in range(relabel_G.number_of_nodes())])
        c = 0.0
        for e in relabel_G.edges():
            if edge_func(relabel_G, e, theta):
                re_G.add_edge(e[0], e[1])
                c += 1
        components = list(nx.connected_components(re_G))
        total_t.append(time.time() - start_time)
        print('theta=', theta)
        if communities is not None:
            p,r,f1 = evaluate1(relabel_G, components, communities, False)
            precisions.append(p)
            recalls.append(r)
            f1s.append(f1)
            print('precision, recall and F1:', round(p, 2), round(r, 2), round(f1, 2))
        if additional_metrics:
            modulars.append(nx_comm.modularity(relabel_G,components))
            num_comp.append(len(components))
            percent_k3.append(sum(nx.triangles(re_G).values()) / 3 / ori_num_k3)
#             num_wedge = 0
#             for e in re_G.edges():
#                 num_wedge += len(set(re_G[e[0]])|set(re_G[e[1]]))
#             num_wedge = num_wedge-sum(nx.triangles(re_G).values())
#             percent_wedge.append(num_wedge/ori_num_wedge)
            edges_remove.append(c/relabel_G.number_of_edges())
            tmp_cond, tmp_cond_big = [], []
            for comp in components:
                tmp_cond.append(float(nx.cut_size(relabel_G, comp))/nx.volume(relabel_G, comp))
                if len(comp)>2:
                    tmp_cond_big.append(tmp_cond[-1])
            cond.append(np.average(tmp_cond))
            big_comp_cond.append(np.average(tmp_cond_big))
            print("modularity:", modulars[-1], "Num of components:",len(components), "percent of k3:", round(percent_k3[-1],4), "percent of edges:", round(c/relabel_G.number_of_edges(),4), )
    print('avg running time', np.average(total_t))
    return metrics
    
    
def plot(metrics, xticks=None):
    plt.plot(np.array(metrics["conductance"]), label='conductance', marker='d')
    plt.plot(np.array(metrics["big_comp_cond"]), label='conductance(big)', marker='d')
    plt.plot(np.array(metrics["num_components"])/max(metrics["num_components"]), 
             label='normalized num_comp', marker='d')
    plt.plot(np.array(metrics["edges_removed"]), label='percent_edges_left', marker='d')
    plt.plot(np.array(metrics["num_k3"]), label='percent_k3_left', marker='d')
    plt.plot(np.array(metrics["modularity"]), label='modulars', marker='d')

    plt.plot(metrics["precision"], label='precision', marker='x')
    plt.plot(metrics["recall"], label='recall', marker='x')
    plt.plot(metrics["f1"], label='F1', marker='x')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    if xticks is not None:
        plt.xticks([i for i in range(len(metrics["f1"]))], xticks, rotation=30)
    plt.show()
    
def WT_edge_func(relabel_G, e, theta):
    return relabel_G[e[0]][e[1]]['wedge']-relabel_G[e[0]][e[1]]['triangle'] < theta

def TW_edge_func(relabel_G, e, theta):
    return relabel_G[e[0]][e[1]]['triangle']-relabel_G[e[0]][e[1]]['wedge'] > theta

def tec_edge_func(relabel_G, e, theta):
    return relabel_G[e[0]][e[1]]['triangle']/float(relabel_G.degree[e[0]]+relabel_G.degree[e[1]]) > theta

def jac_edge_func(relabel_G, e, theta):
    return relabel_G[e[0]][e[1]]['triangle'] / \
        float(relabel_G.degree[e[0]]+relabel_G.degree[e[1]]-relabel_G[e[0]][e[1]]['triangle']) > theta

def jac_var_edge_func(relabel_G, e, theta):
    return (1+relabel_G[e[0]][e[1]]['triangle']) / \
        float(relabel_G.degree[e[0]]+relabel_G.degree[e[1]]-relabel_G[e[0]][e[1]]['triangle']-1) > theta

def k3_edge_func(relabel_G, e, theta):
    return relabel_G[e[0]][e[1]]['triangle'] >= theta