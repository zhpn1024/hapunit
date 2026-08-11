#from zbio import io
import time
#import sys
#import tabix
#from scipy.cluster import hierarchy
#import matplotlib.pylab as plt
#from ete3 import Tree
#from grouping1 import haplogrouping2

ab2 = list(range(1,100))
ab1 = list(map(chr, [96+i for i in range(1, 27)]))
tsd = {'A>G': 1, 'G>A': 1, 'C>T': 1, 'T>C': 1}
w = 2.1
qcth = 4
acth = 4

def distance(hap1, hap2, weights, euclidean = False):
  d = 0
  if len(hap1) != len(hap2):
    print(hap1, len(hap1))
    print(hap2, len(hap2))
    exit(-1)
  for i in range(len(hap1)):
    if hap1[i] == hap2[i]: continue
    if euclidean: d += weights[i] ** 2
    else: d += weights[i]
    #if vs[i] in vts: d += 1
    #else: d += w
  if euclidean: return d ** 0.5
  else: return d #round(d, 2)

def dis_batch(haps, i, weights, euclidean = False):
  out = []
  l = len(haps)
  for j in range(i+1, l):
    out.append(distance(haps[i], haps[j], weights, euclidean = euclidean))
  return out

def getDisMatrix(haps, weights, euclidean = False, nproc = 1): #vs, vts):
  if nproc > 1:
    #import multiprocessing as mp
    from multiprocessing import Pool
    pool = Pool(processes = nproc)
    #manager = mp.Manager()
    #mhaps = manager.list(haps)
  dm = []
  l = len(haps)
  if nproc > 1:
    #b = max(nproc, 10)
    #for i0 in range(0, l-1, b):
    multi_res = [pool.apply_async(dis_batch, (haps, i, weights, euclidean,)) for i in range(l-1)] # range(i0, min(i0+b, l-1))
    for res in multi_res:
      dm += res.get()
    pool.close()
  else:
    for i in range(l-1):
    #if nproc > 1:
    #  multi_res = [pool.apply_async(distance, (haps[i], haps[j], weights, euclidean,)) for j in range(i+1, l)]
    #  dm += [res.get() for res in multi_res]
    #else:
      dm += dis_batch(haps, i, weights, euclidean) #[distance(haps[i], haps[j], weights, euclidean = euclidean) for j in range(i+1, l)]
    #for j in range(i+1, l):
      #dm.append(distance(haps[i], haps[j], weights, euclidean = euclidean))
  return dm

def get_newick(node, parent_dist, leaf_names, newick=''):
  if node.is_leaf():
    return "%s:%.3f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
  else:
    if len(newick) > 0:
      newick = "):%.3f%s" % (parent_dist - node.dist, newick)
    else:
      newick = ");"
    newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
    newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
    newick = "(%s" % (newick)
    return newick

def tree_group(t, gnodes, counts, j, lc, maxc):
  freq_node = []
  single = True
  for name in gnodes[j]:
    if gnodes[j][name] != '': continue
    tl = t.search_nodes(name=name)
    if len(tl) == 0: continue
    t1 = tl[0]
    if j + 1 < lc:
      s1 = tree_group(t1, gnodes, counts, j+1, lc, maxc)
      if not s1: single = False
    else:
      tn = len(t1)
      if tn > maxc[j+1]: maxc[j+1] = tn
    freq_node.append((counts[t1.name], t1.name))
  freq_node.sort(reverse = True)
  ln = len(freq_node)
  if ln >= 2: single = False
  if ln > 26: print('Too many subnodes! j: {} len: {} freq_node: {}'.format(j, len(freq_node), freq_node))
  if ln > maxc[j]: maxc[j] = ln
  for i, f in enumerate(freq_node):
    c, name = f
    if single:
      gnodes[j][name] = ''
      continue
    if j == 0:
      if i < len(ab1): gnodes[j][name] = ab1[i].upper()
      else: gnodes[j][name] = ab1[i - len(ab1)]
    else:
      if c < acth: gnodes[j][name] = ''
      elif i >= len(ab1): gnodes[j][name] = ''
      elif j % 2 == 0: gnodes[j][name] = ab1[i]
      else: gnodes[j][name] = ab2[i]
  return single

def get_dist(t, name):
  if name == 'i1': return name, 0
  d = t.get_distance('i1', name)
  return name, d

def haplogrouping(t, dynamic = False, nproc = 1, cutoff = 0.7):
  cutoffs = [cutoff, cutoff/2, cutoff/4, cutoff/8] # [0.7, 0.35, 0.7/4, 0.7/8]
  nproc = 1 ##
  if nproc > 1:
    from multiprocessing import Pool
    pool = Pool(processes = nproc)

  ### inner node names
  nodes_dict = {}
  inner_nodes, leaf_nodes = [], []
  n = 1
  for node in t.traverse():
    if not node.is_leaf():
      node.name = 'i{}'.format(n)
      n += 1
      inner_nodes.append(node.name)
    else:
      leaf_nodes.append(node.name)
    nodes_dict[node.name] = node
  print('{} haplogrouping 1'.format(time.ctime()))
  ###  node counts
  counts = {}
  for name in leaf_nodes:
    c = int(name.split('_')[-1])
    counts[name] = c
    node = nodes_dict[name]
    while node.up:
      node = node.up
      if node.name not in counts: counts[node.name] = 0
      counts[node.name] += c
  print('{} haplogrouping 2'.format(time.ctime()))
  ### node distance
  leaf_dis = []
  node_dis = {}
  #node_dis[name] = t.get_distance("i1", name)
  #print('{} haplogrouping 2.0'.format(time.ctime()))
  if nproc == 1:
    for name in nodes_dict:
      if name == "i1": node_dis[name] = 0
      else: node_dis[name] = t.get_distance("i1", name)
      #if nodes_dict[name].is_leaf(): leaf_dis.append(node_dis[name])
  else:
    #node_dis['i1'] = 0
    #names = [name for name in nodes_dict if name != 'i1']
    #print('{} haplogrouping 2.1'.format(time.ctime()))
    multi_res = [pool.apply_async(get_dist, (t, name, )) for name in nodes_dict] # (t.get_distance, ('i1', name,)) for name in names]
    print('{} haplogrouping 2.2'.format(time.ctime()))
    for res in multi_res:
      name, d = res.get() # names[i]
      node_dis[name] = d # res.get()
      #if nodes_dict[name].is_leaf(): leaf_dis.append(node_dis[name])
    print('{} haplogrouping 2.3'.format(time.ctime()))
    pool.close()

  for name in leaf_nodes: leaf_dis.append(node_dis[name])

  height = sum(leaf_dis) / len(leaf_dis)
  #print(height, leaf_dis)
  print('{} haplogrouping 3'.format(time.ctime()))
  norm_dis = {}
  for name in nodes_dict:
    norm_dis[name] = 1 - node_dis[name] / height
  ### grouping
  lc = len(cutoffs)
  gnodes = {j: {} for j in range(lc)} # nodes for grouping
  j_dict = {}
  print('{} haplogrouping 4'.format(time.ctime()))
  for node in t.traverse():
    d = norm_dis[node.name]
    if node.up is not None: j0 = j_dict[node.up.name]
    else: j0 = 0
    if not dynamic:
      j = 0
      while j < lc:
        if d >= cutoffs[j]: break
        j += 1
      j_dict[node.name] = j
      if j == j0: continue
      for j1 in range(j0, j):
        gnodes[j1][node.name] = ''
    else: # dynamic cutoffs
      if j0 == 0: co = 1 * cutoffs[0]
      elif j0 >= len(cutoffs):
        j_dict[node.name] = j0
        continue
      else:
        nu = node.up
        while nu.name not in gnodes[j0-1]: nu = nu.up
        co = norm_dis[nu.name] * cutoffs[0]
      if d >= co:
        j_dict[node.name] = j0
      else:
        j_dict[node.name] = j0 + 1
        gnodes[j0][node.name] = ''
        if node.is_leaf():
          for j in range(j0+1, lc):
            gnodes[j][node.name] = ''
  maxc = [0, 0, 0, 0, 0]
  tree_group(t, gnodes, counts, 0, lc, maxc)
  print('maxc: {}'.format(maxc))
  print('{} haplogrouping 5'.format(time.ctime()))
  #print(gnodes)
  typed_dict = {}
  for name0 in leaf_nodes:
    type = ''
    name = name0
    for j in range(lc-1, -1, -1):
      while name not in gnodes[j]:
        name = nodes_dict[name].up.name
      type = '{}{}'.format(gnodes[j][name], type)
    typed_dict[name0] = type
    #print(name0, type)
  return typed_dict

def hapgrp(chr, start, stop, logfile):
  #print('Loading variants...')
  vts = {}
  vinfo = {}
  for lst in io.splitIter('varids_merge_{}_{}_{}.txt'.format(chr, start, stop), skip=1):
    lst2 = lst[1].split(':')
    if lst2[2] in tsd: vts[lst[0]] = lst[1]
    vinfo[lst[0]] = lst[1]

  #print('Loading sample info...')
  batch = {}
  for lst in io.splitIter('/Parastor300s_G30S/lius/genealogy_tree/mega-tree/batch_test/1kg_hgdp_nyuwa_7410_popinfo.txt'):
    if lst[2] == 'China': batch[lst[0]] = 'NyuWa'
    elif len(lst[2]) == 3: batch[lst[0]] = '1KGP'
    else: batch[lst[0]] = 'HGDP'

  for lst in io.splitIter('/Parastor300s_G30S/lius/Haplotype/haplotype_merge/chr{}_merge_own_1kgp_hgdp_siteandhap.vcf.gz'.format(chr)):
    if lst[0].startswith('#CHROM'):
      samples = lst[9:]
      break

  #print('Loading haplotypes...')
  idhaps, idcnts = {}, {}
  allhaps, allids = [], []
  goodhaps = []
  qchaps = {'1KGP':[], 'HGDP':[], 'NyuWa':[]}
  qcdict = {}
  selecthaps, selectids = [], []
  for lst in io.splitIter('haps_merge_{}_{}_{}.txt'.format(chr, start, stop)):
    if lst[0].startswith('#'):
      vs = lst[0][1:].split(',')
      if len(vs) <= 3:
        logfile.write('Too few variants ({}) in the block, skip...\n'.format(len(vs)))
        return
    else:
      lst2 = lst[1].split('_')
      n = int(lst2[2])
      if n >= 10: goodhaps.append([lst[0], lst[1]])
      elif n <= 1: break ####
      if lst[3] != '':
        md = None
        for hap, hid in goodhaps:
          d = distance(lst[0], hap, vs, vts)
          if md is None or md > d:
            md = d
            mhap = hap
        dvs = {}
        qcdict[lst[1]] = dvs
        for i in range(len(mhap)):
          if mhap[i] == lst[0][i]: continue
          pos = int(vinfo[vs[i]].split(':')[1])
          dvs[pos] = [] #.append(int(pos))
        for hid2 in lst[3].split(','):
          lst3 = hid2.split('_')
          sid = samples[int(lst3[0][1:])]
          qchaps[batch[sid]].append([sid, lst[1], dvs])
      else:
        selecthaps.append(lst[0])
        selectids.append(lst[1])

      idhaps[lst[1]] = lst[0]
      idcnts[lst[1]] = lst[2]
      allhaps.append(lst[0])
      allids.append(lst[1])

  #print('Performing haplotype QC...')
  rawfiles = {'1KGP': '/Parastor300s_G30S/Resources/1KGP3_origin_vcf/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr{}.recalibrated_variants.vcf.gz'.format(chr),
              'HGDP': '/parastor300/Resource/hgdp_wgs/hgdp_wgs.20190516.full.chr{}.vcf.gz'.format(chr),
              'NyuWa': '/parastor300/zhangp/haplotype/phasing/pop5/norm/chr{}.vcf.gz'.format(chr)
             }
  for b in qchaps:
    for lst in io.splitIter(rawfiles[b]):
      if lst[0].startswith('#CHROM'):
        rawsamples = lst
        break
    sdict = {}
    pdict = {}
    for sd in qchaps[b]:
      sdict[sd[0]] = {}
      for p in sd[2]:
        pdict[p] = 1
    if len(pdict) == 0: continue
    qcstart, qcstop = min(pdict)-1, max(pdict)
    qci = {sid: i for i, sid in enumerate(rawsamples) if sid in sdict}
    tb1 = tabix.open(rawfiles[b])
    for data in tb1.query('chr'+chr, start-1, stop):
      pos = int(data[1])
      if pos not in pdict: continue
      format = data[8].split(':')
      if b != 'NyuWa': qi = format.index('GQ')
      else: qi = format.index('PL')
      for sd in qchaps[b]:
        if pos not in sd[2]: continue
        i = qci[sd[0]]
        gq = data[i].split(':')[qi]
        if gq == '.': gq = 0
        if b != 'NyuWa': sd[2][pos].append(int(gq))
        else:
          if len(data[3]) > 1 or len(data[4]) > 1: continue
          qs = eval('[{}]'.format(gq))
          qs.sort()
          if len(qs) > 1 and qs[0] == 0: sd[2][pos].append(qs[1])
          else: sd[2][pos].append(qs[0])

  fail = 0
  for hid in allids:
    if hid not in qcdict: continue
    good = True
    for pos in qcdict[hid]:
      gqs = qcdict[hid][pos]
      if len(gqs) == 0 or max(gqs) < 20:
        good = False
        logfile.write('{} {} {}\n'.format(hid, pos, gqs)) #print(hid, pos, gqs)
        break
    if good:
      selectids.append(hid)
      selecthaps.append(idhaps[hid])
    else:
      fail += 1
  logfile.write('{} haplotypes failed QC, {} selected haps.\n'.format(fail, len(selectids)))

  #print('Output data...')
  outfile = open('groups_merge_{}_{}_{}.txt'.format(chr, start, stop), 'w')
  outfile.write('#{}\n'.format(','.join(vs)))
  for i, h in enumerate(selecthaps):
    outfile.write(io.tabjoin(h, selectids[i]) + '\n')

  #print('Calculating distance matrix...')
  dm = getDisMatrix(selecthaps, vs, vts)

  #print('Calculating linkage...')
  Z = hierarchy.linkage(dm, 'average')
  tree = hierarchy.to_tree(Z, False)

  nwk = get_newick(tree, tree.dist, selectids, "")
  outfile.write('#nwk\n')
  outfile.write(nwk)
  outfile.write('\n')

  t = Tree(nwk)
  #print('Grouping...')
  typed_dict = haplogrouping(t)
  outfile.write('#grouping\n')
  outfile.write(str(typed_dict))
  outfile.write('\n')

  for node in t.traverse():
    if node.name in typed_dict: node.name = '{} {} {}'.format(typed_dict[node.name], node.name, idcnts[node.name][0:150])

  outfile.write('#grouped tree\n')
  outfile.write(str(t))
  outfile.write('\n')

