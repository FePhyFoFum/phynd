#!/usr/bin/env python3
import sys,os,math
import tree_reader
import scipy.stats
import subprocess
import node
import seq
import tree_utils
import gzip
import argparse
from sink_tree import sink_tree

remove_intermediate_files = True

class Bipart:
    def __init__ (self,lf,rt):
        self.left = lf
        self.right = rt

    def __str__(self):
        x = ",".join(list(self.left))
        y = ",".join(list(self.right))
        return x+" | "+y
    
    def conflict(self, inbp):
        if len(inbp.right.intersection(self.right)) > 0 and len(inbp.right.intersection(self.left)) > 0:
            if len(inbp.left.intersection(self.right)) > 0 and len(inbp.left.intersection(self.left)) > 0 :
                return True
        if len(inbp.left.intersection(self.left)) > 0 and len(inbp.left.intersection(self.right)) > 0:
            if len(inbp.right.intersection(self.left)) > 0 and len(inbp.right.intersection(self.right)) > 0:
                return True
        return False

    def same(self, inbp):
        if len(inbp.right) != len(self.right) and len(inbp.right) != len(self.left):
            return False
        if inbp.right == self.right and inbp.left == self.left:
            return True
        if inbp.right == self.left and inbp.left == self.right:
            return True
        return False

"""
return one bipart for the node
"""
def get_bipart(nd,rt):
    out = set(rt.lvsnms())
    right = set(nd.lvsnms())
    left = set(out-right)
    bp = Bipart(right,left)
    return bp

"""
get all the biparts for a tree
"""
def get_biparts(rt):
    bps = []
    out = set(rt.lvsnms())
    for i in rt.iternodes():
        if len(i.children) == 0:
            continue
        if i == rt:
            continue
        right = set(i.lvsnms())
        left = set(out-right)
        bp = Bipart(right,left)
        bps.append(bp)
    return bps

"""
basic iqtree run with alrt 1000 bb 10000
"""
def run_iqtree(fn, numthreads):
    #ORIGINAL
    #cmd = "iqtree -s "+fn+" -bb 1000 > log"
    cmd = "iqtree -nt "+str(numthreads)+" -m GTR+G -s "+fn+" -alrt 1000 -nm 10000 -bb 1000 -redo > log"
    os.system(cmd)
    sink_tree(fn+".treefile",fn+".treefile.sunk")
    otfs = [fn+".ckp.gz",fn+".bionj",fn+".log",fn+".iqtree",fn+".mldist",fn+".model.gz",fn+".contree",fn+".splits.nex"]
    for otf in otfs:
        if os.path.exists(otf):
            os.remove(otf)
    return fn+".treefile.sunk"

"""
make the iqtrees for an interval
"""
def make_trees(seqs,fn,a_len,sw,step,numthreads):
    #print mmcutoff
    num_trees = [] # range steps are the list and then sets are the values
    segs = []
    count = 0
    for i in range(0,a_len,step):
        j = a_len if i+sw > a_len else i+sw
        if (j-i)/float(sw) < .9:
            continue
        tfn = open(fn+"."+str(count),"w")
        print(" making tree ( range:",i,"-",j,")",file=sys.stderr)
        seqnames = [] #have to have 4 seqs
        for k in seqs:
            ts = k.seq[i:j]
            if len(ts)-(ts.count("-")+ts.upper().count("N")) < 300:
                continue
            else:
                tfn.write(">"+k.name+"\n"+ts+"\n")
                seqnames.append(k.name)
        tfn.close()
        if len(seqnames) < 4:
            ttf = open(fn+"."+str(count)+".treefile.sunk","w")
            ttf.write("("+",".join(seqnames)+");")
            ttf.close()
        else:
            run_iqtree(fn+"."+str(count),numthreads)
            if remove_intermediate_files:
                otf = fn+"."+str(count)
                if os.path.exists(otf):
                    os.remove(otf)
                otf = fn+"."+str(count)+".treefile"
                if os.path.exists(otf):
                    os.remove(otf)
        tfn = fn+"."+str(count)+".treefile.sunk"
        num_trees.append(tfn)
        segs.append((i,j))
        count += 1
        if j == a_len:
            break
    return num_trees,segs

"""
get the index in a list of biparts for a specific bipart
"""
def get_bp_ind(mlbps,bp):
    for i in range(len(mlbps)):
        if mlbps[i].same(bp):
            return i

"""
add a bipart to a set of biparts
  ignores same biparts
"""
def add_bp(inset,inbp):
    add = True
    for i in inset:
        if inbp.same(i):
            add = False
            break
    if add:
        inset.add(inbp)
    return inset

"""
write the rscript
"""
def write_r():
    s = """args <- commandArgs(trailingOnly = TRUE)
a = read.table(args[1],sep=" ",h=T)
library(gplots)
png(args[2])
heatmap.2(as.matrix(a),dendrogram = 'none',sepwidth=c(0.05,0.05), sepcolor="white",
          lhei = c(0.2,1),margins=c(9,10),main=args[3],key=F,colsep=c(1:10000),rowsep=(1:10000), 
          Rowv=FALSE,Colv=FALSE,trace='none', cexRow=1.0,cexCol = 1,labRow="",
          breaks=seq(-1,1,length=4),col=c("red","white","blue"))
dev.off()
"""
    out = open("rplot.r","w")
    out.write(s)
    out.close()

"""
process the interval trees for whether there are conflicts
plots to r
"""
def run_bp_window(infn,tsegfiles,mltr,segc,outf):
    write_r()
    mlto = tree_reader.read_tree_file_iter(mltr).__next__()
    mlbps = get_biparts(mlto)
    segs = {}
    segslong = {}
    segstree = {}
    plotsegs = [] # each row is a seg, each column is a node
    conflictsegscount = []
    count = 0
    for i in range(len(tsegfiles)):
        conflictcount = 0
        segs[count] = []
        segslong[count] = set()
        plotsegs.append([0] * len(mlbps))
        cmd = "bp -c "+mltr+" -t "+tsegfiles[i]+" -tv"
        segstree[count] = open(tsegfiles[i],"r").readline()
        o = subprocess.check_output(cmd.split(" "),stderr=subprocess.STDOUT)
        keepo = str(o).split("\\n")
        cf = keepo[-8]
        cc = keepo[-7]
        cft = tree_reader.read_tree_string(cf)
        cct = tree_reader.read_tree_string(cc)
        for j,k in zip(cft.iternodes(),cct.iternodes()):
            if len(j.children) > 1:
                sbp = None 
                inde = None
                if j.label != "" or k.label != "":
                    sbp =get_bipart(j,cft)
                    inde = get_bp_ind(mlbps,sbp) 
                # conflict one
                if j.label != "":
                    if int(j.label) > 0:
                        conflictcount += 1
                        segs[count].append(sbp)
                        plotsegs[i][inde] = -1
                        # process the bp out from above to record the actual split that conflicts
                        start = False
                        for l in keepo:
                            if start:
                                if "  (" == l[0:3]:
                                    tttt = tree_reader.read_tree_string(l.strip().split(" ")[-1])
                                    segslong[count] = add_bp(segslong[count],get_biparts(tttt)[0])
                            if "read " == l[0:5]:
                                start = True
                            if "TREES " == l[0:6]:
                                break
                # concordant one
                if k.label != "":
                    if int(k.label) > 0:
                        plotsegs[i][inde] = 1
        if remove_intermediate_files:
            otf = tsegfiles[i]
            if os.path.exists(otf):
                os.remove(otf)
        conflictsegscount.append(conflictcount) #just a running tally of the conflicts per segment 
        count += 1

    # print the number of conflicts per segment
    print(infn+" "+" ".join([str(k) for k in conflictsegscount]))

    # print the verbose stuff to the gzip
    detfile = gzip.open(outf+".details.gz","wt")
    for i,sc,sl,t in zip(segs,segc,segslong,segstree):
        detfile.write (str(i)+" "+"-".join([str(k) for k in list(sc)])+"\n")
        for j in segs[i]:
            detfile.write  (" conflicts with:"+str(j)+"\n")
        for j in segslong[i]:
            detfile.write  (" prefers:"+str(j)+"\n")
        detfile.write (" tree:"+segstree[t]+"\n")
    detfile.close()

    # write the plotting information
    ouf = open(outf,"w")
    first = True
    for sc in segc:
        if first == True:
            first = False
        else:
            ouf.write(" ")
        ouf.write("\""+"-".join([str(k) for k in list(sc)])+"\"")
    ouf.write("\n")
    for i in range(len(mlbps)):
        s = []
        for j in range(len(plotsegs)):
            s.append(str(plotsegs[j][i]))
        ouf.write(" ".join(s)+"\n")
    ouf.close()
    cmd = "Rscript rplot.r "+outf+" "+outf+".png "+infn[0:min(15,len(infn))]+" > rlog 2>&1"
    os.system(cmd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-s", "--seqfile", help="sequence file",required=True)
    parser.add_argument("-w", "--sw", help="sliding window (e.g., length of segment)", type=int,required=True)
    parser.add_argument("-i", "--increment", help="increment (e.g., every i basepair)",type=int, required=True)
    parser.add_argument("-p", "--plot", help="should we plot or just make the trees (default = False)",
                        action="store_true", required=False, default=False)
    parser.add_argument("-t","--threads",help="how many threads for iqtree?",required=False,type=int,default=2)
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    args = parser.parse_args()

    seqfile = args.seqfile
    print("running "+seqfile,file=sys.stderr)
    seqs = seq.read_fasta_file(seqfile)
    sw = args.sw #sliding window
    step = args.increment #steps
    numthreads = args.threads
    a_len = len(seqs[0].seq)
    if a_len < sw:
        print("sequence is less than "+str(args.sw),file=sys.stderr)
        sys.exit(0)
    if a_len < sw+step:
        print("wouldn't make more than one segment: "+str(a_len)+"<"+str(sw+step),file=sys.stderr)
        sys.exit(0)
    x,segc = make_trees(seqs,seqfile,a_len,sw,step,numthreads)
    y = run_iqtree(seqfile,numthreads)
    
    if args.plot == True:
        print(" writing outputs and figures",file=sys.stderr)
        run_bp_window(seqfile,x,y,segc,seqfile+".interval_plotdata")
