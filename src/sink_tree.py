import sys
import tree_reader

def sink_tree(filen,outfn):
    t = tree_reader.read_tree_file_iter(filen).__next__()
    cutoffa = 80
    cutoffb = 95
    nds = set()
    for i in t.iternodes():
        if len(t.children) < 2:
            continue
        l = i.label
        if "/" in l:
            s = l.split("/")
            a = float(s[0])
            b = float(s[1])
            if a < cutoffa or b < cutoffb:
                nds.add(i)
            i.label = ""
    for j in nds:
        chs = j.children
        par = j.parent
        par.remove_child(j)
        for k in chs:
            k.parent = par
            par.add_child(k)
    outf = open(outfn,"w")
    outf.write(t.get_newick_repr(False)+";")
    outf.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print ("python3 "+sys.argv[0]+" filen outn")
        sys.exit(0)
    sink_tree(sys.argv[1],sys.argv[2])