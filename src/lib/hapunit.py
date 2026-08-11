from . import block
import time

try: cmp
except:
  def cmp(a, b):
    return (a > b) - (a < b)

def changechr(chr):
  '''change between two chr versions
  '''
  if chr[0].isdigit() or chr in ('X','Y','M'): return 'chr' + chr
  elif chr == 'MT' : return 'chrM'
  elif chr == 'chrM' : return 'MT'
  elif chr[0:3] == 'chr' : return chr[3:]
  else : return chr

def export_corr(corr, filename, data = None, sel = None):
  outfile = open(filename, 'w')
  from hapunit.zbio import io
  l = len(corr)
  if data is None:
    ids = list(range(l))
  elif sel is not None:
    ids = ['{}_{}'.format(i, data[i].pos) for i in sel]
  else:
    ids = ['{}_{}'.format(i, data[i].pos) for i in range(l)]
  outfile.write(io.tabjoin('corr', ids, '\n'))
  for i in range(l):
    outfile.write(io.tabjoin(ids[i], corr[i], '\n'))

def get_ld2(haps, afs, total, i, j):
  p1, p2 = afs[i], afs[j]
  c12 = 0
  for h, c in haps.items():
    if c[0] == 0: continue
    if h[i] == 0 or h[j] == 0: continue
    c12 += c[0]
  p12 = c12 / total
  if min(p1, p2) <= 0 or max(p1, p2) >= 1 or p12 <= max(0, p1+p2-1) or p12 >= min(p1, p2):
    return 1.0 #
  d = p12 - p1 * p2
  if d > 0: m = min((1-p1)*p2, p1*(1-p2))
  else: m = min((1-p1)*(1-p2), p1*p2)
  return (d/m) ** 2

def ld2_batch(haps, afs, total, i, sel = None, half = True):
  l = len(afs)
  if half: start = i + 1
  else: start = 0
  if sel is None: return [get_ld2(haps, afs, total, i, j) if j > i else None for j in range(start, l)]
  else:
    dsel = set(sel) #{}
    #for j in sel: dsel[j] = 1
    return [get_ld2(haps, afs, total, i, j) if j > i or j not in dsel else None for j in range(start, l)]

def get_d(haps, afs, total, i, j, norm = False):
  p1, p2 = afs[i], afs[j]
  c12 = 0
  for h, c in haps.items():
    if c[0] == 0: continue
    if h[i] == 0 or h[j] == 0: continue
    c12 += c[0]
  p12 = c12 / total
  d = p12 - p1 * p2
  if not norm: return d
  if min(p1, p2) <= 0 or p12 <= max(0, p1+p2-1):
    return -1.0
  elif max(p1, p2) >= 1 or p12 >= min(p1, p2):
    return 1.0 #
  if d > 0: m = min((1-p1)*p2, p1*(1-p2))
  else: m = min((1-p1)*(1-p2), p1*p2)
  return d / m

def get_r(haps, afs, total, i, j):
  p1, p2 = afs[i], afs[j]
  v1, v2 = p1 * (1-p1), p2 * (1-p2)
  if v1 <= 0 or v2 <= 0: return 'NA'
  c12 = 0
  for h, c in haps.items():
    if c[0] == 0: continue
    if h[i] == 0 or h[j] == 0: continue
    c12 += c[0]
  p12 = c12 / total
  r = (p12 - p1 * p2) / (v1 * v2) ** 0.5
  return r



class Var:

  sample_w = None
  #popi = None
  #popskip = {}

  def __init__(self, lst = None, copy = None):
    if lst is not None:
      self.chr, self.pos = lst[0], int(lst[1])
      self.ref, self.alt = lst[3], lst[4]
      self.vid, self.info = lst[2], lst[7]
      self.gts = [None] * (len(lst) - 9)
      #print(lst[1], lst[9])
      self.ac, self.nh = 0, 0
      for i in range(9, len(lst)):
      #print(i, lst[i], )
        i2 = i-9
        self.gts[i2] = tuple(map(int, lst[i].split(':')[0].split('|')))  # ensure phased GTs
        if self.sample_w is not None:
          if self.sample_w[i2] == 0: continue
        self.ac += 2 - self.gts[i2].count(0) # sum(self.gts[i2])
        self.nh += 2
      #self.ac = sum([sum(gt) for gt in self.gts])
      self.af = float(self.ac) / self.nh #2 / len(self.gts)
    elif copy is not None:
      self.chr, self.pos, self.ref, self.alt, self.ac, self.nh, self.af = copy.chr, copy.pos, copy.ref, copy.alt, copy.ac, copy.nh, copy.af
      self.vid, self.info = copy.vid, copy.info
      self.gts = copy.gts[:]
    self._lendis = None
  
  def __len__(self):
    return len(self.gts)

  def __str__(self):
    return '{}:{}:{}>{}'.format(self.chr, self.pos, self.ref, self.alt)

  def __repr__(self):
    return  '{} {} samples'.format(str(self), len(self.gts))

  def num_haps(self):
    return self.nh
    #if self.popi is None:
    #  return len(self.gts) * 2
    #else:

  def hap_gt(self, i):
    return self.gts[i//2][i%2]

  def __cmp__(self, other):
    return cmp(self.chr, other.chr) or cmp(self.pos, other.pos) or cmp(self.alt, other.alt)

  def __lt__(self, other):
    c = self.__cmp__(other)
    return c < 0

  def __gt__(self, other):
    c = self.__cmp__(other)
    return c > 0

  def maf(self):
    if self.af <= 0.5: return self.af
    else: return 1 - self.af

  def is_indel(self):
    if len(self.ref) > 1: return True
    if len(self.alt) > 1: return True
    return False

  def is_multi(self):
    return self.alt.find(',') > 0

  def num_alts(self):
    alts = self.alt.split(',')
    return len(alts)

  def split_alts(self):
    if not self.is_multi(): return [Var(copy = self)]
    alts = self.alt.split(',')
    out = []
    for i, a in enumerate(alts):
      va = Var(copy = self)
      va.alt, va.ac = a, 0
      gi = i + 1
      for j, gt in enumerate(self.gts):
        va.gts[j] = tuple([1 if gi == g else 0 for g in gt])
        if self.sample_w is not None:
          if self.sample_w[j] == 0: continue
        va.ac += va.gts[j].count(1)
      va.af = float(va.ac) / va.nh
      out.append(va)
    return out

  def split_len(self, i, shorter = False, alt = None):
    #lendis = self.len_dis()
    alts = self.alt.split(',')
    va = Var(copy = self)
    if alt is None:
      if shorter: va.alt = '<TR:S{}>'.format(i)
      else: va.alt = '<TR:L{}>'.format(i)
    else: va.alt = alt
    va.ac = 0
    alleles = [self.ref] + alts
    if shorter: al = [1 if len(a) <= i else 0 for a in alleles]
    else: al = [1 if len(a) >= i else 0 for a in alleles]
    for j, gt in enumerate(self.gts):
      va.gts[j] = tuple(al[g] for g in gt)
      if self.sample_w is not None:
        if self.sample_w[j] == 0: continue
      va.ac += va.gts[j].count(1)
    va.af = float(va.ac) / va.nh
    return va


  def len_dis(self): 
    if self._lendis is not None: return self._lendis
    alts = self.alt.split(',')
    lens = [len(self.ref)] + [len(a) for a in alts]
    lendis = [0] * (max(lens) + 1)
    for j, gt in enumerate(self.gts):
      if self.sample_w is not None:
        if self.sample_w[j] == 0: continue
      for g in gt: lendis[lens[g]] += 1
    self._lendis = lendis
    return lendis

  def infer_period(self):
    lendis = self.len_dis()
    lm = len(lendis)
    peaks = []
    n = 0
    for i, c in enumerate(lendis):
      if c == 0: continue
      n += 1
      if i > 0 and c < lendis[i-1]: continue
      if i < lm - 1 and c < lendis[i+1]: continue #peaks.append(i)
      peaks.append(i)
    print(n, len(peaks))
    if len(peaks) <= 1: return 1
    elif len(peaks) < n / 2: return 1
    dps = {}
    for j in range(len(peaks) - 1):
      d = peaks[j+1] - peaks[j]
      if d not in dps: dps[d] = 0
      dps[d] += 1
    if len(dps) == 1: return d
      #if n >= 3: return 1
      #else: return d
    da = []
    for d in dps: da.append([dps[d], d])
    da.sort(reverse=True)
    print(da, n)
    if da[0][1] < da[1][1]: return da[0][1]
    elif da[1][0] > 2 and da[1][0] > da[0][0] / 2: return da[1][1]
    else: return da[0][1]

  def split_tr(self, outlier_th = 0.05):
    alts = self.alt.split(',')
    if len(alts) < 3: return self.split_alts()
    lendis = self.len_dis()
    lm, total = len(lendis), sum(lendis)
    oth = total * outlier_th
    out = []
    i = il = 0
    s = lendis[0]
    for i in range(lm):
      s += lendis[i]
      if s < oth: il, sl = i, s
      if s > 0.5 * total: break
    if sl > 0: # short outliers
      vl = self.split_len(il, shorter=True, alt='<TR:SO{}>'.format(il))
      out.append(vl)
    sr = 0 # lendis[lm-1]
    for ir in range(lm-1, -1, -1):
      sr += lendis[ir]
      if sr >= oth: break # irl, srl = ir, sr
    sr -= lendis[ir]
    ir += 1
    if sr > 0: # long outliers
      vr = self.split_len(ir, shorter=False, alt='<TR:LO{}>'.format(ir))
      out.append(vr)

    im1 = im2 = i
    sm1 = sm2 = s
    while im1 > 1 and lendis[im1-1] <= lendis[im1]:
      sm1 -= lendis[im1]
      im1 -= 1
      if lendis[im1] == 0: break
    im, sm = None, None
    if sm1 > sl: im, sm = im1, sm1

    while im2 < lm-1 and lendis[im2+1] <= lendis[im2]:
      im2 += 1
      sm2 +=lendis[im2]
      if lendis[im2] == 0: break
    sm2r = total - sm2
    if sm2r > sr:
      if im is None: im, sm = im2 + 1, sm2r 
      elif sm2r > sm1: im, sm = im2 + 1, sm2r
    if im is not None:
      vm = self.split_len(im, shorter=False, alt='<TR:L{}>'.format(im))
      out.append(vm)
    return out

class BlockData:

  sample_w = None

  def __init__(self, afth = 0.2, skip = [10, 5]):
    self.data = []
    self.ld2 = []
    self.select = []
    self.skipnum = 0
    self.sel_step, self.sel_lim = 50, 100
    self.afth, self.skip = afth, skip
    self.nhaps = None

  def __len__(self):
    return len(self.data)

  def last_pos(self):
    if len(self.data) > 0: return self.data[-1].pos
    else: return -1

  def _select(self):
    var = self.data[-1]
    if var.af < self.afth: return False
    self.skipnum += 1
    sel = False
    if self.skipnum <= self.skip[1]: sel = True
    if self.skipnum >= self.skip[0]: self.skipnum = 0
    return sel

  def add_var(self, var, select = None):
    self.data.append(var)
    #if len(self.data) >= 2 and self.data[-1].pos < self.data[-2].pos:
      #self.data.sort()
    if select is None:
      if self._select():
        self.select.append(len(self.data)-1)
    elif select:
      self.select.append(len(self.data)-1)

  def remove(self, n):
    if n <= 0: return
    if n == len(self.data):
      i = len(self.select)
    else:
      i = 0
      for i in range(len(self.select)):
        if self.select[i] >= n: break
    self.data[0:n] = []
    self.ld2[0:n] = []
    #for i in range(len(self.select)):
      #if self.select[i] >= n: break
    #else: i = len(self.select)
    self.select[0:i] = []
    #print('remove n', n, 'sel', i)
    for i in range(len(self.select)):
      self.select[i] -= n

  def get_haps(self, sel_iter, all = False): #popi = None, popskip = {}):
    sel = list(sel_iter)
    haps = {}
    if len(sel) == 0: return haps
    #sd = [self.data[i] for i in self.select]
    l = len(self.data[0].gts) * 2 # self.data[0].num_haps()
    for j in range(l):
      w = 1 # weight for popskip
      if self.sample_w is not None:
        j2 = j // 2
        w = self.sample_w[j2]
        #if self.popi[j2] is None: w = 0 # continue
        #elif len(popskip) > 0:
        #  for p in self.popi[j2]:
        #    if p in popskip:
        #      w = 0
        #      break
          #if ps: continue
        if not all and w == 0: continue
      h = tuple([self.data[i].hap_gt(j) for i in sel])
      if h not in haps: haps[h] = [w, j, [j]]
      else:
        haps[h][0] += w
        haps[h][2].append(j)
      #print('l: {}, w: {}, j: {}, j2: {}'.format(l, w, j, j2))
    return haps

  def check_ld2(self):
    l = len(self.data) #  the first is self corr
    if len(self.ld2) < l:
      d = l - len(self.ld2)
      self.ld2 += [[] for i in range(d)]
    for i in range(l):
      d = l - i - len(self.ld2[i]) # triangle matrix
      if d > 0:
        self.ld2[i] += [None] * d

  #def ld2_batch(self, i): pass

  def all_ld2(self, nproc, sel = None):
    l = len(self)
    alls = list(range(l))
    haps = self.get_haps(alls)
    afs = [self.data[i].af for i in alls]
    if sel is None: sel = alls
    total = 0.0 # float(self.data[0].nh)
    for h, c in haps.items():
      total += c[0]
    self.check_ld2()
    if nproc > 1:
      from multiprocessing import Pool
      pool = Pool(processes = nproc)
      multi_res = [pool.apply_async(ld2_batch, (haps, afs, total, i, sel, False, )) for i in sel] # range(l-1)]
      for j, res in enumerate(multi_res):
        i = sel[j]
        r = res.get()
        self.ld2[i][1:] = r[i+1:] # res.get() ##
        for i1 in range(i):
          if r[i1] is None: continue
          j1 = i - i1
          self.ld2[i1][j1] = r[i1]
      pool.close()
    else:
      for i in sel: # range(l-1):
        r = ld2_batch(haps, afs, total, i, sel, False) ##
        self.ld2[i][1:] = r[i+1:] # ld2_batch(haps, afs, total, i, ) # sel, )
        for i1 in range(i):
          if r[i1] is None: continue
          j1 = i - i1
          self.ld2[i1][j1] = r[i1]

  def haps_ld2(self, haps, sel, i, j, keepnone = False):
    i1, j1 = sel[i], sel[j]
    #print(i, j, i1, j1, len(self.ld2), len(self.ld2[i1]))
    if keepnone or self.ld2[i1][j1-i1] is not None: return self.ld2[i1][j1-i1]

    v1, v2 = self.data[i1], self.data[j1]
    af1 = [1-v1.af, v1.af]
    af2 = [1-v2.af, v2.af]
    t = 0 # = v1.num_haps()
    cnts = [[0,0], [0,0]]
    for h, c in haps.items():
      if c[0] == 0: continue
      g1, g2 = h[i], h[j] #sd[k1][i], sd[k2][i]
      cnts[g1][g2] += c[0] #w[i]
      t += c[0]

    d = cnts[1][1] / float(t) - af1[1] * af2[1]
    if d > 0: m = min(af1[0] * af2[1], af1[1] * af2[0])
    else: m = min(af1[0] * af2[0], af1[1] * af2[1])
    ld2 = (d/m) ** 2

    #m, mi, mj = 1, 0, 0
    #fe = [[1,1], [1,1]]
    #for i in range(2):
    #  for j in range(2):
    #    fe[i][j] = cnts[i][j] / (t*af1[i]*af2[j])
    #    if m > fe[i][j]: m, mi, mj = fe[i][j], i, j
    #ld2 = (1-m) ** 2

    #if self.ld2[i1][j1-i1] is not None and self.ld2[i1][j1-i1] != ld2:
      #print('ld2 inconsistent: {} {} {} {} {} {} {}'.format(ld2, self.ld2[i1][j1-i1], i1, j1, i, j, sel))
    self.ld2[i1][j1-i1] = ld2
    #print('ld2: {}, p1: {}, p2: {}, i1: {}, j1: {}, cnts: {}, i: {}, j: {}'.format(ld2, self.data[i1].pos, self.data[j1].pos, i1, j1, cnts, i, j))
    return ld2 #1 - f/f0 #float(cm) / (m + cm)

  def haps_dr(self, haps, sel, i, j):
    i1, j1 = sel[i], sel[j]
    v1, v2 = self.data[i1], self.data[j1]
    af1 = [1-v1.af, v1.af]
    af2 = [1-v2.af, v2.af]
    var1, var2 = af1[0]*af1[1], af2[0]*af2[1]
    t = 0 # = v1.num_haps()
    cnts = [[0,0], [0,0]]
    for h, c in haps.items():
      if c[0] == 0: continue
      g1, g2 = h[i], h[j] #sd[k1][i], sd[k2][i]
      cnts[g1][g2] += c[0] #w[i]
      t += c[0]
    p12 = float(cnts[1][1]) / t
    d = p12 - af1[1] * af2[1]
    if var1 <= 0 or var2 <= 0: return d, 'NA'
    r = d / (af1[0]*af1[1]*af2[0]*af2[1]) ** 0.5
    return d, r

  def get_ld2(self, sel_iter, keepnone = False): #, popi = None, popskip = {}):
    sel = list(sel_iter)
    haps = self.get_haps(sel) #, popi=popi, popskip=popskip)
    self.check_ld2()
    l = len(sel)
    outld2 = [[None] * l for i in range(l)] # squire matrix
    for i in range(l):
      outld2[i][i] = 1
      for j in range(i+1, l):
        outld2[i][j] = outld2[j][i] = self.haps_ld2(haps, sel, i, j, keepnone)
        #if self.ld2[sel[i]][sel[j]] is None: pass
        #else:
          #outld2[i][j] = self.ld2[sel[i]][sel[j]]
    return outld2

  def export_cluster(self, i1, i2): #, popi = None, popskip = {}):
    print('export_cluster get_ld2')
    ld2 = self.get_ld2(range(i1, i2)) #, popi=popi, popskip=popskip)
    print('export_cluster block.get_blocks')
    try: cluster, rpd = block.get_blocks(ld2)
    except:
      export_corr(ld2, 'errcorr_{}_{}_{}_{}.txt'.format(self.data[0].chr, self.data[0].pos, i1, i2))
      exit()
    print('export_cluster out')
    out = []
    la = len(self.data)
    for i, c in enumerate(cluster):
      c1 = [c[0]+i1, c[1]+i1]
      #if la - c1[1] < 10: break
      c = tuple(c)
      s = ''
      if c in rpd:
        rp = rpd[c]
        s = ','.join(['{}:{}'.format(j+i1, self.data[j+i1].pos) for j in rp])
      out.append([self.data[c1[0]].chr, self.data[c1[0]].pos, self.data[c1[1]-1].pos, c1, s]) ## c1 - 1
    return out

  def find_blocks(self, clear = False, keepflank = 10, remove = True): #, popi = None, popskip = {}):
    if clear: remove = True
    ls, la = len(self.select), len(self.data)
    print(time.ctime(), 'pos1', self.data[0].pos, 'pos2', self.data[-1].pos, 'ls', ls, 'la', la)
    sld2 = self.get_ld2(self.select) #, popi=popi, popskip=popskip)
    sel_cluster, sel_rp = block.get_blocks(sld2)

    selstep = 20
    sel_start, sel_last, bstart = 0, 0, 0
    outblocks = []
    si = 0
    while sel_last < ls:
      if si < len(sel_cluster): sel_last = sel_cluster[si][0] # next cluster
      else: sel_last = ls - 1 # no next cluster
      # the blank before sel block
      end = False
      while(sel_last - sel_start >= selstep):
        bstop = self.select[sel_start + selstep]
        print('blank bstart', bstart, 'bstop', bstop, 'sel_start', sel_start, 'sel_last', sel_last)
        out = self.export_cluster(bstart, bstop) #, popi=popi, popskip=popskip)
        for o in out:
          if bstop - o[3][1] < keepflank:
            bstart = max(o[3][0] - keepflank, bstart)
            print('blank bstart = max(o[3][0] - keepflank, bstart)', bstart)
            end = True
            break
          outblocks.append(o)
          bstart = o[3][1] ##
          print('blank bstart = o[3][1]', bstart)
          print(o)
        else:
          bstart = max(bstop - keepflank, bstart)
          print('blank bstart = max(bstop - keepflank, bstart)', bstart)
        while sel_start < len(self.select) and bstart > self.select[sel_start]:
          sel_start += 1
        if end: break

      # the block with flank
      if sel_last + keepflank >= ls: break
      if sel_start >= sel_cluster[si][1]:
        si += 1
        continue
      sel_last = sel_cluster[si][1] + keepflank
      print('sel_cluster', si, sel_cluster[si])
      if sel_last >= ls: break
      bstop = self.select[sel_last]
      print('block bstart', bstart, 'bstop', bstop, 'sel_start', sel_start, 'sel_last', sel_last)
      out = self.export_cluster(bstart, bstop) #, popi=popi, popskip=popskip)
      print('block out', out)
      for o in out:
        if bstop - o[3][1] < keepflank:
          bstart = max(o[3][0] - keepflank, bstart)
          print('block bstart = max(o[3][0] - keepflank, bstart)', bstart)
          break
        outblocks.append(o)
        bstart = o[3][1] ##
        print('block bstart = o[3][1]', bstart)
        print(o)
      else:
        bstart = max(bstop - keepflank, bstart)
        print('block bstart = max(bstop - keepflank, bstart)', bstart)
      while sel_start < len(self.select) and bstart > self.select[sel_start]:
        sel_start += 1
      si += 1

    # release analyzed data
    if remove: self.remove(bstart)
    # output all
    if clear:
      bstart, bstop = 0, len(self)
      print('clear bstart', bstart, 'bstop', bstop)
      out = self.export_cluster(bstart, bstop) #, popi=popi, popskip=popskip)
      for o in out:
        print(o)
      outblocks += out
      self.remove(bstop)

    return outblocks
