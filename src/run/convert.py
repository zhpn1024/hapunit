#!/bin/env python3
from hapunit.lib import hapunit, hapgrp, block
from hapunit.zbio import io, stat
import time
import tabix
from scipy.stats import chi2
import sys

def set_parser(parser):
  parser.add_argument("-i", type=str, dest="input", required=True, help="Input hapunit group vcf data")
  parser.add_argument("-o", type=str, dest="output", required=True, help="Output hapunit group 0/1 genotype vcf")
  #parser.add_argument("-p", type=int, dest="numproc", default=1, help="Number of processes (1)")
  parser.add_argument("--region", type=str, default=None, help="Chr region to process, default: all")
  parser.add_argument("--filt", type=str, default=None, help="Ignore the listed variant sites")
  parser.add_argument("--minaf", type=float, default=0.005, help="Minor allele frequency threshold, default: 0.005")
  parser.add_argument("--minac", type=int, default=3, help="Minimum HUG allele count, default: 3")
  parser.add_argument("--chrmap", type=str, help="Input chromosome id mapping table file if chr ids are not the same in different files")
  parser.add_argument("--pop", type=str, help="Input sample population, format: sample pop1 pop2")
  parser.add_argument("--popskip", type=str, help="Ignore certain pops in total counts, comma separated")
  parser.add_argument("--popexclude", type=str, help="Do not export samples in these populations, comma separated")
  parser.add_argument("--popinclude", type=str, help="Only export samples in these populations, comma separated")
  parser.add_argument("--exclude", type=str, help="Do not export these samples in the sample id list file")
  parser.add_argument("--include", type=str, help="Only export these samples in the sample id list file")
  #parser.add_argument("--pair", action="store_true", help="Paired end")

tsd = {'A>G': 1, 'G>A': 1, 'C>T': 1, 'T>C': 1}

vcftpl = ['', '', '', '', '<HUG>', '.', '.', '', 'GT:HUG']

def chunhe_test(n, m, k, alt = 'two.tailed'): #pass
  p = 0
  for i in range(m+1):
    p += stat.hypergeo(n*2, n, m, i) * stat.hypergeo_test(n, i, m-i, k, alt=alt)
  return p

def HWtest(AA, Aa, aa):
  n = float(AA + Aa + aa)
  if n == 0: return 1
  pA = (2*AA + Aa)/(2*n)
  pa = 1 - pA # (Aa + 2*aa)/(2*n)
  EAA = n*(pA**2)
  EAa = n*(2*pA*pa)
  Eaa = n*(pa**2)
  if (EAA>=5) and (EAa>=5) and (Eaa>=5) and n>=40:
    chisqr = ((AA-EAA)**2)/EAA + ((Aa-EAa)**2)/EAa + ((aa-Eaa)**2)/Eaa
    p = chi2.sf(chisqr, 1)
  else:
    m = 2*aa + Aa
    n = int(n)
    if aa > Eaa: p = 2 * chunhe_test(n, m, aa, alt='g')
    else: p = 2 * chunhe_test(n, m, aa, alt='l')
    if p > 1: p = 1
  return p, EAa


def run(args):

  chrmap = {}
  if args.chrmap is not None:
    for lst in io.splitIter(args.chrmap, sep=None):
      if len(lst) < 2: continue
      chrmap[lst[0]] = lst[1]
      chrmap[lst[1]] = lst[0]

  filt = {}
  if args.filt is not None:
    print('{} Loading filt data {}...'.format(time.ctime(), args.filt))
    for lst in io.splitIter(args.filt):
      if lst[0] not in filt:
        filt[lst[0]] = {} #chr
        if lst[0] in chrmap: filt[chrmap[lst[0]]] = filt[lst[0]]
      filt[lst[0]][lst[1]] = 1

  samples = []

  chr = ''
  if args.region is not None: chr = args.region.split(':')[0]

  #pops = {}
  if args.pop is not None:
    pops = {}
    for lst in io.splitIter(args.pop, sep=None):
      pops[lst[0]] = lst[1:]
    #hapunit.BlockData.pops = hapunit.Var.pops = pops

  popskip = {}
  if args.popskip is not None:
    for s in args.popskip.split(','):
      popskip[s] = 1
    #hapunit.BlockData.popskip = hapunit.Var.popskip = popskip

  vcf = tabix.open(args.input)
  outvcf = open(args.output, 'w')
  #outree = open('{}_unit.bed'.format(args.output), 'w')
  #outhaps = open('{}_haps.bed'.format(args.output), 'w')

  headeradd = {'ALT': '##ALT=<ID=HUG,Description="Hapunit Group">',
               'INFO': '##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">\n##INFO=<ID=HUG,Number=.,Type=String,Description="Hapunit group subclass">',
               'FORMAT': '##FORMAT=<ID=HUG,Number=1,Type=String,Description="Hapunit Group (phased)">',
              } ## not used
  headeradd = {
               'INFO': '##INFO=<ID=CGT,Number=G,Type=Integer,Description="Sample counts for all genotypes">\n##INFO=<ID=HetExp,Number=1,Type=Float,Description="Expect of No. of heterozygosity samples">\n##INFO=<ID=HWeq,Number=1,Type=Float,Description="Hardy-Weinberg equilibrium p-value">',
              }  ##
  if args.popskip is not None:
    headeradd['INFO'] += '##INFO=<ID=ACALL,Number=A,Type=Integer,Description="Allele count including skipped samples">\n##INFO=<ID=AFALL,Number=A,Type=Float,Description="Allele Frequency including skipped samples">' 
  vcfheader = []
  for lst in io.splitIter(args.input):
    if lst[0].startswith('##'): vcfheader.append(io.tabjoin(lst))
    else: break
    lst2 = lst[0].split('=')
    if len(lst2) <= 1: continue
    key = lst2[0][2:]
    if key in headeradd:
      vcfheader.append(headeradd[key])
      del headeradd[key]
  for key in headeradd:
    vcfheader.append(headeradd[key])
  vcfheader.append('##hapunit_convertCommand=' + ' '.join(sys.argv))
  popi = sample_w = None
  if args.pop is not None:
    popi = [None] * (len(lst) - 9)
    sample_w = [1] * (len(lst) - 9)
  popn = {}
  if lst[0].startswith('#'):
    vcfheader.append(io.tabjoin(lst))
    if args.pop is not None:
      for i in range(9, len(lst)):
        s = lst[i]
        i2 = i - 9
        if s in pops:
          popi[i-9] = set(pops[s])
          for p in set(pops[s]):
            if p in popskip: sample_w[i2] = 0
            if p not in popn: popn[p] = 0
            popn[p] += 1
        else: sample_w[i2] = 0
      hapunit.BlockData.sample_w = hapunit.Var.sample_w = sample_w
      print('len(sample_w): {}, sum(sample_w): {}'.format(len(sample_w), sum(sample_w)))
      #print('popi: {}\nsample_w: {}'.format(popi, sample_w))
  for s in vcfheader:
    outvcf.write(s + '\n')

  n = len(lst) - 9
  l = n #
  for lst in io.splitIter(args.input):
    if lst[0].startswith('#'): continue

    infos = lst[7].split(';')
    end = gid = symbol = ''
    for s in infos:
      lst2 = s.split('=')
      if lst2[0] == 'END': end = lst2[1]
      elif lst2[0] == 'GID': gid = lst2[1]
      #elif lst2[0] == 'SYMBOL': symbol = lst2[1]
      #if lst2[0] == 'HUG': break
        #hugl = lst2[1].split('|')
    format = lst[8].split(':')
    hi = format.index('HUG')
    hugs, gts = {}, {}
    for i, s in enumerate(lst):
      if i < 9: continue
      sd = s.split(':')[hi].split('|')
      for j in range(2):
        hg = sd[j]
        for k in range(len(hg)):
          g = hg[0:k+1]
          if g not in hugs:
            hugs[g] = [0, 0] # ac & acall
            if gid != '': vid = '{}_{}_{}_{}'.format(lst[0], lst[1], gid, g)
            else: vid = '{}_{}_{}'.format(lst[0], lst[1], g)
            gts[g] = [lst[0], lst[1], vid, lst[3], '<HUG>', '.', '.', 'END={};HUG={}'.format(end, g), 'GT'] + ['0|0'] * l
          if sample_w is None or sample_w[i-9] > 0:
            hugs[g][0] += 1
          hugs[g][1] += 1
          if j == 0: gts[g][i] = '1' + gts[g][i][1:]
          else: gts[g][i] = gts[g][i][:2] + '1'

    hugl = sorted(hugs)
    for g in hugl:
      if hugs[g][0] < args.minac: continue
      af = hugs[g][0] / float(l * 2)
      if af < args.minaf: continue
      afall = hugs[g][1] / float(l * 2)

      cgt = [0,0,0] # 0/0, 0/1, 1/1
      cgtp = {p: [0,0,0] for p in popn}
      for i, gt in enumerate(gts[g][9:]):
        if gt == '0|0': igt = 0 # cgt[0] += 1
        elif gt in ('0|1', '1|0'): igt = 1 # cgt[1] += 1
        else: igt = 2 # cgt[2] += 1
        cgt[igt] += 1
        if args.pop is not None:
          for p in popi[i]: cgtp[p][igt] += 1
      hw, ea = HWtest(cgt[0], cgt[1], cgt[2])
      if hw < 0.01: hw = '%.3e' % hw
      else: hw = round(hw, 4)
      gts[g][7] += ';AC={};AF={};CGT={};HetExp={};HWeq={}'.format(hugs[g][0], round(af, 4), ','.join(map(str, cgt)), round(ea, 1), hw)
      if args.pop is not None:
        for p in popn:
          hw, ea = HWtest(cgtp[p][0], cgtp[p][1], cgtp[p][2])
          if hw < 0.01: hw = '%.3e' % hw
          else: hw = round(hw, 4)
          gts[g][7] += ';CGT_{0}={1};HetExp_{0}={2};HWeq_{0}={3}'.format(p, ','.join(map(str, cgtp[p])), round(ea, 1), hw)

      if args.popskip is not None:
        afall = hugs[g][1] / float(l * 2)
        gts[g][7] += ';ACALL={};AFALL={}'.format(hugs[g][1], round(afall, 4))

      outvcf.write(io.tabjoin(gts[g], '\n'))


if __name__ == '__main__':
  import sys, argparse
  p = argparse.ArgumentParser()
  set_parser(p)
  if len(sys.argv)==1:
    print(p.print_help())
    exit(0)
  run(p.parse_args())
