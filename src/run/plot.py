#!/bin/env python3
from hapunit.lib import hapunit, hapgrp, block
from hapunit.zbio import io, fa, tools, plot
import time
import tabix
import numpy as np
from scipy.cluster import hierarchy
#import matplotlib.pylab as plt
from ete3 import Tree
import sys

def set_parser(parser):
  parser.add_argument("-i", type=str, dest="input", required=True, help="Input haplotype bed file from the output of group module")
  parser.add_argument("-r", type=str, dest="region", required=True, help="Input haplotype block region (format: chr:start-stop)")
  parser.add_argument("-o", type=str, dest="output", required=True, help="Output figure file (pdf)")
  #parser.add_argument("-p", type=int, dest="numproc", default=1, help="Number of processes (1)")
  #parser.add_argument("--minvars", type=int, default=3, help="Minimum haplotype variants filter, default: 3")
  #parser.add_argument("--chr", type=str, default=None, help="Chr to process, default: auto")
  parser.add_argument('--use', type=str, help='Use only the listed variant sites')
  parser.add_argument("--filt", type=str, default=None, help="Ignore the listed variant sites")
  parser.add_argument("--minaf", type=float, default=0.05, help="Minor allele frequency threshold, default: 0.05")
  parser.add_argument("--figsize", type=str, help="Figure size (default: 10,8)")
  parser.add_argument('--title', type=str, help='Figure title')
  #parser.add_argument("--chrmap", type=str, help="Input chromosome id mapping table file if chr ids are not the same in different files")
  #parser.add_argument("--weightmaf", action="store_true", help="Weight distance by MAF")
  #parser.add_argument("--weightmaf05", action="store_true", help="Weight distance by MAF + 0.5")
  #parser.add_argument("--tsweight", type=float, default=1, help="Weight for transition SNP distances, default: 1")
  #parser.add_argument("--tvweight", type=float, default=1, help="Weight for transversion SNP distances, default: 1")
  #parser.add_argument("--indelweight", type=float, default=1, help="Weight for indel distances, default: 1")
  #parser.add_argument("--cutoff", type=float, default=0.7, help="Cutoff for grouping, range 0-1, lower cutoffs result in more branches, default: 0.7")
  parser.add_argument("--pop", type=str, help="Use only samples with selected population labels, comma separated")
  parser.add_argument("--popskip", type=str, help="Ignore certain pops in total counts, comma separated")
  parser.add_argument('--cmpop', type=str, help='Plot two figures to compare distribution between two populations, format: grp1pop1;grp1pop2,grp2pop1..')
  #parser.add_argument("--euclidean", action="store_true", help="Using euclidean distance for grouping")
  #parser.add_argument("--dynamic", action="store_true", help="Using dynamic grouping cutoffs")
  #parser.add_argument('--tr', action='store_true', help='Use TR length (longer, long outlier, short outlier) as alleles')
  parser.add_argument("--nosv", action="store_true", help="Skip SV variants")
  parser.add_argument("--plot", action="store_true", help="Plot tree dendrogram for each hapunit")
  parser.add_argument("-p", type=int, dest="numProc", default=1, help="Number of processes")
  #parser.add_argument("--gid", type=int, help="Column of gene ID in block file, 0-based")
  #parser.add_argument("--symbol", type=int, help="Column of gene symbol in block file, 0-based")
  #parser.add_argument("--addchr", action="store_true", help="Add chr in chr names")
  #parser.add_argument("--pair", action="store_true", help="Paired end")
  parser.add_argument('--markgrp', type=str, help='Mark groups in the tree, comma separated')

tsd = {'A>G': 1, 'G>A': 1, 'C>T': 1, 'T>C': 1}

vcftpl = ['', '', '', '', '<HUG>', '.', '.', '', 'GT:HUG']

maxsel = 200

def var_type(ref, alt):
  if len(ref) > 1 or len(alt) > 1: return 'indel'
  key = '{}>{}'.format(ref, alt)
  if key in tsd: return 'ts'
  else: return 'tv'

def in_iter(di, use, filt):
  for lst in di:
    chr, pos = lst[0], int(lst[1])
    if pos in filt[chr]: continue
    #if args.nosv and lst[4] in ('<DUP>', '<DEL>', '<INS>', '<INV>'): continue
    if use is not None: # chr in use:
      if not use[chr].find(pos):
        #print('pos {} not found in use[{}], {} {}'.format(pos, chr, use[chr].current(), use[chr].i))
        continue # pos not in use[chr]: continue
    yield lst

def getvar(lst, args):
  #lst = l.rstrip('\n').split('\t')
  if args.chr is not None and lst[0] != args.chr: return None ##
  if args.nosv and lst[4] in ('<DUP>', '<DEL>', '<INS>', '<INV>'): return None ##
  #if len(lst[3]) > 1 or len(lst[4]) > 1: return None # continue
  var = hapunit.Var(lst)
  if var.is_multi(): return None
  if var.maf() < args.minaf: return None # continue
  return var

def shared(s1, s2):
  l1, l2 = len(s1), len(s2)
  l = min(l1, l2)
  for i in range(l, 0, -1):
    s = s1[0:i]
    if s2.startswith(s): return s
  else: return ''

def run(args):
  #var_weight = {'ts': args.tsweight, 'tv': args.tvweight, 'indel': args.indelweight}
  #if args.cutoff >= 1 or args.cutoff <= 0: args.cutoff = 0.7 # default

  chrmap = {}
  filt = {}
  if args.filt is not None:
    print('{} Loading filt data {}...'.format(time.ctime(), args.filt))
    for lst in io.splitIter(args.filt):
      if lst[0] not in filt:
        filt[lst[0]] = {} #chr
        if lst[0] in chrmap: filt[chrmap[lst[0]]] = filt[lst[0]]
      filt[lst[0]][int(lst[1])] = 1

  use = None # {}
  if args.use is not None:
    print('{} Loading use data {}...'.format(time.ctime(), args.use))
    use = {}
    pc, pp = 0, 1
    if args.use.split('.')[-1].lower() == 'bim': pp = 3
    for lst in io.splitIter(args.use):
      if lst[pc] not in use:
        use[lst[pc]] = io.OrderedList() # {} #chr
        if lst[pc] in chrmap: use[chrmap[lst[pc]]] = use[lst[pc]]
      use[lst[pc]].data.append(int(lst[pp])) # = 1
    for chr in use: use[chr].check()

  samples = []

  chr = ''
  pop = {}
  if args.cmpop is None:
    if args.pop is not None:
      for s in args.pop.split(','):
        pop[s] = 0
  else:
    lst = args.cmpop.split(',')
    for j, sg in enumerate(lst):
      for s in sg.split(';'): pop[s] = j
    #for lst in io.splitIter(args.pop, sep=None):
    #  pops[lst[0]] = lst[1:]
    #hapunit.BlockData.pops = hapunit.Var.pops = pops

  popskip = {}
  if args.popskip is not None:
    for s in args.popskip.split(','):
      popskip[s] = 1
    #hapunit.BlockData.popskip = hapunit.Var.popskip = popskip

  inhaps = tabix.open(args.input)
  #outvcf = open('{}_gt.vcf'.format(args.output), 'w')
  #outree = open('{}_unit.bed'.format(args.output), 'w')
  #outhaps = open('{}_haps.bed'.format(args.output), 'w')

  lst = args.region.split(':')
  lst2 = lst[1].split('-')
  chr, start, stop = lst[0], lst2[0], lst2[1]
  data = {}
  dgroup = {}
  for lst in inhaps.query(chr, int(start), int(stop)):
    if lst[1] != start or lst[2] != stop: continue
    if lst[3] == 'variants': vars = lst[4].split(',')
    elif lst[3] == 'hapunit_tree': nwk = lst[4] # t = Tree(nwk)
    else:
      lst2 = lst[3].split('|')
      if len(lst2) < 2: continue
      grp, cnt = lst2[1], int(lst2[0].split('_')[1])
      cnts = [cnt, 0]
      #if grp not in data: data[grp] = []
      #data[grp].append([cnt, lst[4]])
      if args.pop is not None or args.cmpop is not None:
        cnts = [0, 0]
        for s in lst[5].split(','):
          lstp = s.split(':')
          if lstp[0] in pop: cnts[pop[lstp[0]]] += int(lstp[1])
        #if cnt == 0: continue
      data[lst2[0]] = [cnts, lst[4]]
      dgroup[lst2[0]] = lst2[1]

  vari = [i for i, s in enumerate(vars) if float(s.split(':')[3]) >= args.minaf]
  print(','.join([vars[i] for i in vari]))
  t = Tree(nwk)
  hapdata, dist = [[], []], {}
  t.name = 'root'
  for node in t.traverse(strategy='postorder'):
    if node.name in dgroup:
      node.group = dgroup[node.name]
      cnts, hap = data[node.name]
      hap1 = ''.join([hap[i] for i in vari])
      dis = t.get_distance('root', node.name)
      for j in range(2):
        if cnts[j] > 0: hapdata[j].append([node.group, cnts[j], hap1, dis])
    else:
      chd = node.get_children()
      if chd[0].group == chd[1].group: node.group = chd[0].group
      else:
        gshare = shared(chd[0].group, chd[1].group)
        node.group = gshare
        #print(gshare, chd[0].group, chd[1].group)
        for i in range(2):
          if chd[i].group != gshare:
            chd[i].name = chd[i].group
            dist[chd[i].group] = t.get_distance('root', chd[i].group)


  #hapdata = []
  #for grp in data:
  #  for d in data[grp]:
  #    cnt, hap = d
  #    hapdata.append([grp, cnt, hap])
  #hapdata.sort()

  #print(args.markgrp)
  figsize = 10, 8
  if args.figsize is not None: figsize = tuple(eval(args.figsize))
  fig = plot.figure(figsize=figsize)
  if args.cmpop is None:
    title = args.title
    if title is None:
      title = args.region
      if args.pop is not None:
        title += ' ' + args.pop
    #fig = plot.figure(figsize=figsize)
    #fig.set_size_inches(10, 8)
    gs = plot.GridSpec(1, 4, wspace=0.0)
    tree_ax = plot.subplot(gs[0,0])
    main_ax = plot.subplot(gs[0,1:4], sharey=tree_ax)
    #tree_ax = plot.subplot(gs[0,0],sharey=main_ax)
    plot.hapShow(hapdata[0], main_ax, tree_ax, dist=dist, markgrp=args.markgrp, title = title)
    plot.savefig(args.output)
  else:
    gs = plot.GridSpec(1, 10, wspace=0.0)
    if args.title is not None:
      titles = args.title.split(',')
    else:
      titles = [args.region + ' ' + s for s in args.cmpop.split(',')]
    for j in range(2):
      tree_ax = plot.subplot(gs[0,5*j])
      main_ax = plot.subplot(gs[0,1+5*j:4+5*j], sharey=tree_ax)
      #tree_ax = plot.subplot(gs[0,5*j],sharey=main_ax)
      plot.hapShow(hapdata[j], main_ax, tree_ax, dist=dist, markgrp=args.markgrp, title = titles[j])
    plot.savefig(args.output)

'''
  for lst in io.splitIter(args.block):
    print("{} Processing {}:{}-{}".format(time.ctime(), lst[0], lst[1], lst[2]))
    chr = lst[0]
    chr1 = hapunit.changechr(chr)
    if chr not in filt:
      #chr1 = hapunit.changechr(chr)
      if chr1 in filt: filt[chr] = filt[chr1]
      else: filt[chr] = filt[chr1] = {}
    else:
      #chr1 = hapunit.changechr(chr)
      filt[chr1] = filt[chr]
    if use is not None:
      if chr not in use:
      #chr1 = hapunit.changechr(chr)
        if chr1 in use: use[chr] = use[chr1]
        #else: use[chr] = io.OrderedList() # {}
      #else: use[chr1] = use[chr]
      if chr1 not in use:
        if chr in use: use[chr1] = use[chr]
      #print(use)
    start, stop = int(lst[1])-1, int(lst[2])
    skip = []
    if len(lst) > 4 and lst[4] != '':
      for s in lst[4].split(','):
        try: skip.append(int(s.split(':')[1]))
        except: break
    bdata = hapunit.BlockData(afth = 0.05, skip = [1,1])
    try: di = vcf.query(chr, start, stop)
    except: di = vcf.query(chr1, start, stop)
    fdi = in_iter(di, use, filt)
    for var in io.multiProcIter(fdi, getvar, args, args.numProc-1):
    #for data in di: # vcf.query(chr, start, stop):
      #if var.pos in filt[chr]: continue # data[1]
      #if len(lst[3]) > 1 or len(lst[4]) > 1: continue
      #var = hapunit.Var(data)
      print('var pos: {}, af: {}'.format(var.pos, var.af))
      #if var.maf() < args.minaf: continue
      if args.weightmaf: var.weight = var.maf()
      elif args.weightmaf05: var.weight = var.maf() + 0.5
      else: var.weight = var_weight[var_type(var.ref, var.alt)]
      if var.pos in skip or var.is_indel(): bdata.add_var(var, select=False)
      else: bdata.add_var(var)
    l = len(bdata)
    print('{} bdata {}:{}-{} len: {} sel_len: {}, select: {}'.format(time.ctime(), chr, start, stop, l, len(bdata.select), bdata.select))
    if len(bdata.select) < 3: continue ###
    if len(bdata.select) > maxsel: 
      sel0 = bdata.select
      bdata.select = tools.downsample_even(sel0, maxsel)
      print('Sel downsample: {}'.format(bdata.select))
    bdata.all_ld2(nproc = args.numProc, sel = bdata.select) ##
    print('{} get_ld2 1'.format(time.ctime()))
    sld = bdata.get_ld2(bdata.select)
    sarr = np.array(sld)
    #print('{} low_ld 1'.format(time.ctime()))
    selskip = low_ld(sarr, th = args.ldth)
    sel = []
    for si, i in enumerate(bdata.select):
      if si in selskip: continue
      sel.append(i)
    print('{} all skip'.format(time.ctime()))
    arr = np.array(bdata.get_ld2(range(len(bdata)), keepnone=True))
    #print('{} low_ld 2'.format(time.ctime()))
    allskip = low_ld_1(arr, sel = sel, th = args.ldth)
    skip2 = {}
    for i in allskip: # remove low LD indels
      if bdata.data[i].is_indel():
        j = i - 1
        while j >= 0 and bdata.data[j].pos == bdata.data[i].pos:
          if bdata.data[j].is_indel(): skip2[j] = 1
          j -= 1
        j = i + 1
        while j < len(bdata) and bdata.data[j].pos == bdata.data[i].pos:
          if bdata.data[j].is_indel(): skip2[j] = 1
          j += 1
    for j in skip2: allskip[j] = 1
    print('{} allskip: {}'.format(time.ctime(), allskip))
    sel2 = [i for i in range(len(bdata)) if i not in allskip]
    #if args.pop is not None:
    haps = bdata.get_haps(sel2, all = True) # popi=popi, popskip=popskip) # GT_tuple: [count/weight, first hap sample, all samples list]
    weights = []
    for i in range(l):
      #print(i)
      if i in sel2:
        weights.append(bdata.data[i].weight)
        #print('{}\t{}'.format(i, bdata.data[i]))
      else:
        print('{}\t{}\tskipped'.format(i, bdata.data[i]))
    hapl, hapids = [], [] #haps.keys()
    for h in haps:
      hapl.append(''.join(map(str, h))) # h)
      hid = '{}_{}'.format(haps[h][1], haps[h][0])
      hapids.append(hid)
      #print('{}\t{}'.format(h, hid))

    print('{} Calculating distance matrix'.format(time.ctime()))
    dm =  hapgrp.getDisMatrix(hapl, weights, euclidean=args.euclidean, nproc = args.numProc)
    #print('{} hierarchy.linkage'.format(time.ctime()))
    Z = hierarchy.linkage(dm, 'average')
    #print('{} hierarchy.to_tree'.format(time.ctime()))
    tree = hierarchy.to_tree(Z, False)
    #print('{} hapgrp.get_newick'.format(time.ctime()))
    nwk = hapgrp.get_newick(tree, tree.dist, hapids, "")
    print('{} nwk: \n{}'.format(time.ctime(), nwk))
    #outree.write(io.tabjoin(chr, start, stop, str(nwk), '\n'))
    t = Tree(nwk)
    print('{} Grouping'.format(time.ctime()))
    typed_dict = hapgrp.haplogrouping(t, dynamic = args.dynamic, cutoff = args.cutoff) # , nproc = args.numProc)
    #print('typed_dict: {}'.format(typed_dict))
    for node in t.traverse():
      if node.name in typed_dict: node.name = '{} {}'.format(typed_dict[node.name], node.name)
    print('{} #grouped tree\n{}'.format(time.ctime(), t))
    outgt = ['|'] * len(var) # (len(data) - 9)
    hugs = {}
    poplist, pophugs, hapops = [], {}, {}
    for h, hd in haps.items():
      hid = '{}_{}'.format(hd[1], hd[0])
      hapops[hid] = {}
      hg = typed_dict[hid]
      if hg not in hugs: hugs[hg] = hd[0]
      else: hugs[hg] += hd[0]
      for i in hd[2]:
        i2 = i // 2
        if i % 2 == 0: outgt[i2] = hg + outgt[i2]
        else: outgt[i2] += hg
        if popi is None or popi[i2] is None: continue
        for p in popi[i2]:
          if p not in pophugs:
            pophugs[p] = {}
            poplist.append(p)
          if hg not in pophugs[p]: pophugs[p][hg] = 0
          pophugs[p][hg] += 1
          if p not in hapops[hid]: hapops[hid][p] = 0
          hapops[hid][p] += 1
    print('{} out vcf'.format(time.ctime()))
    huglst = sorted(hugs)
    infostr = 'END={}'.format(stop)
    if args.gid is not None: infostr += ';GID={}'.format(lst[args.gid])
    if args.symbol is not None: infostr += ';SYMBOL={}'.format(lst[args.symbol])
    out = [chr, start+1, '.', bdata.data[0].ref[0], '<HUG>', '.', '.', '{};HUG={}'.format(infostr, '|'.join(huglst)), 'GT:HUG'] #vcftpl[:]
    out += ['./.:' + hgt for hgt in outgt]
    outvcf.write(io.tabjoin(out, '\n'))

    print('{} out haps'.format(time.ctime()))
    outhaps.write(io.tabjoin(chr, start, stop, 'variants', ','.join(map(str, ['{}:{}'.format(bdata.data[i], round(bdata.data[i].af, 4)) for i in sel2])), '\n'))
    #outhaps.write(io.tabjoin(chr, start, stop, 'hapunit_linkage', str(Z), '\n'))
    outhaps.write(io.tabjoin(chr, start, stop, 'hapunit_tree', str(nwk), '\n'))
    outhaps.write(io.tabjoin(chr, start, stop, 'hapunit_groups', ','.join(['{}:{}'.format(hg, hugs[hg]) for hg in huglst]), '\n'))
    for p in poplist:
      outhaps.write(io.tabjoin(chr, start, stop, 'groups_{}_{}'.format(p, 2*popn[p]), ','.join(['{}:{}'.format(hg, pophugs[p][hg]) for hg in huglst if hg in pophugs[p]]), '\n'))

    hgcorr, hgcnt = {}, {}
    hapinfo = []
    for h in haps:
      hd = haps[h]
      hid = '{}_{}'.format(hd[1], hd[0])
      #n = float(sum([hapops[hid][p] for p in hapops[hid]]))
      hpl = [(hapops[hid][p]/float(2*popn[p]), hapops[hid][p], p) for p in hapops[hid]]
      hpl.sort(reverse=True)
      hps = ','.join(['{}:{}:{}'.format(ph[2], ph[1], round(ph[0], 3)) for ph in hpl])
      hg = typed_dict[hid]
      hi = [hg, 0-hd[0], hd, h, hid, hps]
      hapinfo.append(hi)

      for i in range(1, len(hg)+1):
        shg = hg[0:i]
        if shg not in hgcorr:
          hgcorr[shg] = [[0,0,0,j,0] for j in range(len(sel2))] # r^2, r, cnt11, j, D'
          hgcnt[shg] = [0,0,0] # cnt, af, pq
        for j, gt in enumerate(h):
          if gt == 1: hgcorr[shg][j][2] += hd[0]
        hgcnt[shg][0] += hd[0]

    for i in sel2:
      vaf = bdata.data[i].af
      bdata.data[i].pq = vaf * (1 - vaf)
    total = float(bdata.data[i].nh)
    for shg in sorted(hgcorr):
      af = hgcnt[shg][0] / total
      pq = af * (1 - af)
      if hgcnt[shg][0] == 0 or pq == 0:
        outhaps.write(io.tabjoin(chr, start, stop, shg, 'GAF={}'.format(af), '\n'))
        continue
      hgcnt[shg][1] = af
      hgcnt[shg][2] = pq
      for j in range(len(sel2)):
        i = sel2[j]
        d = hgcorr[shg][j][2] / total - af * bdata.data[i].af
        try: r = d / (pq * bdata.data[i].pq) ** 0.5
        except ZeroDivisionError:
          print(total, pq, bdata.data[i], bdata.data[i].pq)
          exit(1)
        r2 = r * r
        hgcorr[shg][j][0] = r2
        hgcorr[shg][j][1] = r
        if d > 0: m = min((1-af) * bdata.data[i].af, af * (1-bdata.data[i].af))
        else: m = min((1-af) * (1-bdata.data[i].af), af * bdata.data[i].af)
        hgcorr[shg][j][4] = d / m

      hgcorr[shg].sort(reverse=True)
      ic = 0
      while ic < len(sel2):
        if hgcorr[shg][ic][0] < 0.99: break
        ic += 1
      if ic < 3: ic = min(3, len(sel2))
      vars = []
      for iic in range(ic):
        i = sel2[hgcorr[shg][iic][3]]
        vars.append('{}:{}:{}'.format(bdata.data[i], round(hgcorr[shg][iic][1], 4), round(hgcorr[shg][iic][4], 4)))
      outhaps.write(io.tabjoin(chr, start, stop, shg, 'GAF={};CORVAR={}'.format(round(af, 4), ','.join(vars)), '\n')) ##


    hapinfo.sort()

    for hi in hapinfo:
      outhaps.write(io.tabjoin(chr, start, stop, '{}|{}'.format(hi[4], typed_dict[hi[4]]), ''.join(map(str, hi[3])), hi[5], '\n'))

    if args.plot: #True:
      print('Ploting dendrogram...')
      import matplotlib.pylab as plt
      typedids = ['{} {}'.format(typed_dict[name], name) for name in hapids]
      plt.figure(figsize=(6, 0.15*len(hapl)))
      hierarchy.dendrogram(Z, labels=['{} {}'.format(hap, typedids[i]) for i, hap in enumerate(hapl)], orientation="left") #''.join(map(str, hap))
      plt.savefig('{}_dendrogram_{}_{}_{}.pdf'.format(args.output, chr, start, stop), bbox_inches='tight')


    #break
    #if var.pos > 1000000: break
'''


if __name__ == '__main__':
  import sys, argparse
  p = argparse.ArgumentParser()
  set_parser(p)
  if len(sys.argv)==1:
    print(p.print_help())
    exit(0)
  run(p.parse_args())
  print('{} Completed.'.format(time.ctime()))
