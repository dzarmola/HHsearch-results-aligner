#!/usr/bin/env python2

import time
import sys,glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from matplotlib.patches import Rectangle as rect
from matplotlib.lines import Line2D as line
import argparse

from my_little_hhpred_reader import *

EVAL_CUTOFF = 0.001

def load_data(files):
    """Creates a internally consistent dictionary of format
    {'PFX':{'PFY':(float evalue,[(x_profile_position,y_profile_position),...]),...}...}
    For each profile pair alignment with better evalue is taken.
    Also - filters e-values above cutoff"""
    #TODO: maybe other criterion (e.g. length) is better?
    matrix = {}
    hhrs = []
    lens = {}
    for file in files:
        h = HHpredOutput(file)
        hhrs.append(h)
        lens[h.query] = h.len
        for _h in h.hits:
            if _h.eval<EVAL_CUTOFF:
                matrix[h.query] = matrix.get(h.query,{})
                matrix[_h.target] = matrix.get(_h.target,{})
                if _h.eval < matrix[_h.target].get(h.query,[1000,0,0])[0]:
                    matrix[_h.target][h.query] = (_h.eval,map(lambda x:(x[1],x[0]), _h.match))
                    matrix[h.query][_h.target] = (_h.eval,_h.match)
    return matrix,hhrs,lens

def sort_fams_by_dist(mat):
    """Sorts families by average matching e-value.
    Currently mostly to get a list of families"""
    fams = sorted(mat.keys(), key=lambda x: sum(_[0] for _ in mat[x].values())/len(mat[x].values()) if mat[x] else 1000 )
    for key in mat.values():
        for k in sorted(key,key=lambda x:x[0]):
            if k not in fams:
                fams.append(k)
    return fams

def read_in_clusters(file):
    """Reads MCL dump format"""
    cluster = {}
    with open(file) as input:
        for i,line in enumerate(input):
            for fam in line.split():
                cluster[fam] = i
    return cluster

def save_point(data,labels,clusters,lens):
    global savepoint_name
    with open(savepoint_name,"w",0) as output:
        output.write(str(clusters)+"\n")
        output.write(str(labels)+"\n")
        output.write(str(data)+"\n")
        output.write(str(lens))

def from_save(filename):
    with open(filename) as save:
        clusters,labels,data,lens = map(eval,save.readlines())
    for k,v in lens.items():
        lens[k] = int(v)
    pretty_plot(data,labels,clusters,lens)

def prepare_data(columns,clusters):
    """Transforms from columnwise to rowwise alignment,
    assigns family names as labels"""
    labels = [x[0] for x in columns[0]]
    data = [[] for x in columns[0]]
    for c in columns:
        for i,x in enumerate(c):
            data[i].append(x[1])
    labels,data = zip(*sorted(zip(labels,data), key=lambda x:clusters.get(x[0],len(clusters)) ))
    return data,labels

def pretty_plot(data,labels,clusters,lens):
    global plot_name
    colors = cm.rainbow(np.linspace(0, 1, len(data)))

    fig = plt.figure()
    ax = plt.axes()

    width = 2
    height = 6
    ax.set_xlim(0, len(data[0])*width)
    ax.set_ylim(0, len(data)*height)

    label_colors = cm.gist_ncar(np.linspace(0, 1, len(set(clusters.values()))+1 ))

    for r,row in enumerate(data):
        colors = cm.rainbow(np.linspace(0, 1, lens.get(labels[r],max(row))-2))
        colors = np.vstack((np.array([0.,0.,0.,1.]),colors,np.array([0.,0.,0.,1.])))
        for i,pos in enumerate(row):
            if pos is not None:
                ax.add_patch(rect((i*width,r*height),width,height-2,color=colors[pos-1] ))
                if labels[r] in clusters:
                    ax.add_line(line([i*width,(i+1)*width],[r*height,r*height] , lw=1,color = label_colors[clusters[labels[r]]]))
                    ax.add_line(line([i*width,(i+1)*width],[(r+1)*height-1,(r+1)*height-1] , lw =1, color = label_colors[clusters[labels[r]]]))
    plt.yticks( map(lambda x:x*height+height/2.,range(len(labels))), labels)
    ax.tick_params(axis = 'y', which = 'major', labelsize = 8)
    plt.tight_layout()
    plt.savefig(plot_name,dpi=300)

def comparable(c1,c2):
    """Checks if two columns are comparable - have any common profile"""
    return any(c1[x][1] is not None and c2[x][1] is not None for x in xrange(len(c1)))

def half(it):
    """Deprecated"""
    lit = list(it)
    return sum(lit)*1./len(lit)>=.5

def first_lower(c1,c2):
    """should c1 come first in the alignment"""
    return all(c1[x][1] < c2[x][1]  for x in xrange(len(c1)) if c1[x][1] is not None and c2[x][1] is not None)


def sorter3(columns):
    """Brute force final sorting of columns"""
    columns = sorted(columns, key=lambda x: len([_ for _ in x if _[1] is not None]))
    csorted = []
    wstawiony=0

    while columns:
        c = columns.pop()
        csorted.append(c)
        no_change = 0
        while columns:
            col = columns.pop()
            found_comp = None
            for i,s in enumerate(csorted):
                if comparable(col,s):
                    found_comp = i
                    no_change = 0
                    if first_lower(col,s):
                        csorted.insert(i,col)
                        wstawiony+=1
                        break
            else:
                if found_comp is not None: 
                        csorted.insert(found_comp+1,col)
                        wstawiony+=1
                else:
                    columns.append(col)
                    no_change += 1
                if len(columns)+1 == no_change:
                    break
    return csorted

def fix_columns(columns, fams):
    """Adds missing families to the column"""
    for c,col in enumerate(columns):
        present = [_[0] for _ in col]
        for f in fams:
            if f not in present:
                col.append((f,None))
        columns[c] = sorted(col)

def make_a_graph_from_matrix(matrix):
    """Creates match-length sorted list of graph edges of format:
    (('PFX',x_profile_position),('PFY',y_profile_position))"""
    edges = []
    popular = {}
    evals = []
    for q,ts in matrix.items():
        for t in ts:
            #eval already checked
            eval,match = ts[t]
            for qc,tc in match:
                if qc!=None and tc!=None:
                    edges.append(((q,qc),(t,tc)))
                    evals.append(len(match))
                    popular[(q,qc)] = popular.get((q,qc),0) + 1
                    popular[(t,tc)] = popular.get((t,tc),0) + 1
    evals,edges = map(list,zip(*sorted(zip(evals,edges), reverse=True))) #reverse, cause by length
    # no fixing, since both ways same matching
    return edges,sorted(popular.keys(),key=lambda x: popular[x], reverse=True)[0]

def neigh(edges,vert):
    """Gets neighbours of a vertex, removes travelled edges"""
    new = []
    delete = []
    for i,e in enumerate(edges):
        if vert in e:
            nv = e[1] if e[0]==vert else e[0]
            new.append(nv)
            delete.append(i)
    while delete:
        edges.pop(delete.pop())
    return new

def in_column(column,what):
    """Checks if this profile is already in a given list.
    To fight against graph inconsistencies"""
    return what[0] in [_[0] for _ in column]

def bfs(edges,start,column,columns):
    """Breadth-first search of our graph.
    Only allows adding vertices (positions to alignment columns)
    if they will not create inconsistencies in the alignment,such as:
    >PFX
    14 15
    >PFY
    27 26
    """
    queue = [start]
    pos = None
    while queue:
        cur = queue.pop(0)
        pos = position_in_columns(columns, column+[cur])
        if pos is None:
            continue
        column.append(cur)
        ns = neigh(edges,cur)
        for n in ns:
            if not in_column(column,n) and not in_column(queue,n):
                    pos = position_in_columns(columns, column+queue+[n])
                    if pos is not None:
                        queue.append(n)
    return position_in_columns(columns, column)

def c2d(c):
    return {k:v for k,v in c}

def smaller(c1,c2):
    """Should column c1 be in alignment before c2"""
    d1 = c2d(c1)
    d2 = c2d(c2)
    for k,v in d1.items():
        if v>=d2.get(k,v+1):
            return False
    if not set(d1.keys()).intersection(set(d2.keys())):
        return avg_pos(c1)<=avg_pos(c2)
    return True

def avg_pos(l):
    """Average position index"""
    return sum(_[1] for _ in l)*1./len(l)

def position_in_columns(columns,ncol):
    """Checks if ncol can be positioned in alignment
    (is consistent). If not returns None, else returns
    proper index, where it can be inserted."""
    if not columns:
        return 0
    if niespojnosc(columns,ncol):
        return None
    for i, c in enumerate(columns):
        if (i==0 or smaller(columns[i-1],ncol)) and smaller(ncol,c):
            return i
    else:
        if smaller(columns[-1],ncol):
            return len(columns)
        return None

def cc_nsp(c1,c2):
    """Are two columns inconsistent (some positions
    indicate c1 should be first, some - that c2 should be first"""
    d1 = c2d(c1)
    d2 = c2d(c2)
    comm = set(d1.keys()).intersection(set(d2.keys()))
    return all([
        any(d1[k]<d2[k] for k in comm),
        any(d1[k]>d2[k] for k in comm)])

def niespojnosc(columns,ncol):
    return any(cc_nsp(c,ncol) for c in columns)

def go_through_graph(edges,starter):
    """Runs BFS until all edges have been consumed.
    Each run produces one column of the alignment.
    Columns are partially sorted"""
    columns = []
    while edges:
        column = []
        pos = bfs(edges,starter,column,columns)
        columns.insert(pos,column)
        if edges:
            starter = edges[0][0]
    return columns

def main(files,cluster_file):
    matrix,hhrs,lens = load_data(files)
    clusters = read_in_clusters(cluster_file) if cluster_file else {}
    fqueue = sort_fams_by_dist(matrix)
    edges,starter = make_a_graph_from_matrix(matrix)
    columns = go_through_graph(edges,starter)
    fix_columns(columns,fqueue)
    columns = sorter3(columns)
    data,labels = prepare_data(columns,clusters)
    save_point(data,labels,clusters,lens)
    pretty_plot(data,labels,clusters,lens)


if __name__ == "__main__":
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    parser = argparse.ArgumentParser(description='Merge hhpred alignments')
    parser.add_argument('hhrs', metavar='HHR', type=str, nargs='*', help='.hhr files to be merged')
    parser.add_argument('--hhr_dir', type=str, default="", help='Directory from which all .hhr will be taken')
    parser.add_argument('--eval', type=float, default=0.001, help='Maximum e-value of used hhsearch results')
    parser.add_argument('--plot_from_save', type=str, default="", help='Savepoint from a previous run')
    parser.add_argument('--cluster_file', type=str, default="", help='Cluster defining file: each cluster in one line')
    parser.add_argument('--save_name', type=str, default="savepoint_{}.txt".format(timestamp),
                         help='Custom savepoint name')
    parser.add_argument('--plot_name', type=str, default="merged_alignment_{}.png".format(timestamp),
                         help='Custom plot name')
    args = parser.parse_args()

    savepoint_name = args.save_name

    plot_name = args.plot_name

    just_plot = args.plot_from_save

    cluster_file = args.cluster_file

    EVAL_CUTOFF = args.eval

    if just_plot:
        from_save(just_plot)
    else:
        if not args.hhrs and not args.hhr_dir:
            exit("No .hhr files indicated")
        hhrs = args.hhrs
        d_hhrs = glob.glob("{}/*.hhr".format(args.hhr_dir))
        hhrs = list(set(hhrs + d_hhrs))
        main(hhrs,cluster_file)

