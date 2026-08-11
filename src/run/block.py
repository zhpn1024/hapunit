#!/bin/env python3
from hapunit.lib import hapunit, block
from hapunit.zbio import io
import time

def set_parser(parser):
  parser.add_argument("-i", type=str, dest="input", required=True, help="Input haplotype vcf data")
  #parser.add_argument("-a", type=str, dest="adapt", required=True, help="Input adapter sequence")
  parser.add_argument("-o", type=str, dest="output", required=True, help="Output block positions")
  #parser.add_argument("-p", type=int, dest="numproc", default=1, help="Number of processes (1)")
  parser.add_argument("--minvars", type=int, default=3, help="Minimum haplotype variants filter, default: 3")
  parser.add_argument("--chr", type=str, default=None, help="Chr to process, default: auto")
  parser.add_argument("--filt", type=str, default=None, help="Ignore the listed variant sites")
  parser.add_argument("--minaf", type=float, default=0.05, help="Minor allele frequency threshold, default: 0.05")
  parser.add_argument("--ldth", type=float, default=0.95, help="average LD score (D') in a block, default: 0.95")
  parser.add_argument("--chrmap", type=str, help="Input chromosome id mapping table file if chr ids are not the same in different files")
  parser.add_argument("--pop", type=str, help="Input sample population, format: sample pop1 pop2")
  parser.add_argument("--popskip", type=str, help="Ignore certain pops in total counts, comma separated")
  parser.add_argument("-p", type=int, dest="numProc", default=1, help="Number of processes")

  #parser.add_argument("--pair", action="store_true", help="Paired end")

def mp_iter(input, process, args, p):
  from multiprocessing import Pool
  pool = Pool(p)
  jobs, res = [], []
  for i in input: # range(n):
    if len(jobs) == p:
      if len(res) != '':
        for r in res:
          if r is not None: yield r
      res = [j.get() for j in jobs]
      jobs = []
    jobs.append(pool.apply_async(process, (i, args, )))
  if len(res) != '':
    for r in res:
      if r is not None: yield r
  res = [j.get() for j in jobs]
  for r in res:
    if r is not None: yield r

def sp_iter(input, process, args):
  for i in input:
    r = process(i, args)
    if r is not None: yield r

def getvar(l, args):
  lst = l.rstrip('\n').split('\t')
  #if chr != lst[0]:
  if args.chr is not None and lst[0] != args.chr: return None ##
    #if chr != '':
    #  clear = True
    #chr = lst[0]
  #if chr not in filt:
  #  chr1 = hapunit.changechr(chr)
  #  if chr1 in filt: filt[chr] = filt[chr1]
  #  else: filt[chr] = {}
  #if lst[1] in filt[chr]: continue
  if len(lst[3]) > 1 or len(lst[4]) > 1: return None # continue
  var = hapunit.Var(lst)
  if var.maf() < args.minaf: return None # continue
  return var

def run(args):

  chrmap = {}
  if args.chrmap is not None:
    for lst in io.splitIter(args.chrmap, sep=None):
      if len(lst) < 2: continue
      chrmap[lst[0]] = lst[1]
      chrmap[lst[1]] = lst[0]

  filt = {}
  if args.filt is not None:
    print('Loading filt data {}...'.format(args.filt))
    for lst in io.splitIter(args.filt):
      if lst[0] not in filt:
        filt[lst[0]] = {} #chr
        if lst[0] in chrmap: filt[chrmap[lst[0]]] = filt[lst[0]]
      filt[lst[0]][lst[1]] = 1

  pops = {}
  if args.pop is not None:
    for lst in io.splitIter(args.pop, sep=None):
      pops[lst[0]] = lst[1:]

  popskip = {}
  if args.popskip is not None:
    for s in args.popskip.split(','):
      popskip[s] = 1

  block.block_average = args.ldth
  selbase, selstep = 50, 20
  keepflank = 10
  samples = []
  bdata = hapunit.BlockData()

  chr = ''
  if args.chr is not None: chr = args.chr
  maxdis = 1000000
  selim = selbase
  clear = False

  infile = io.lineIter(args.input)
  outfile = open(args.output, 'w')

  for l in infile: # io.splitIter(args.input):
    if l.startswith('##'): continue
    if l.startswith('#'):
      lst = l.rstrip('\n').split('\t')
      samples = lst[9:]
      popi = sample_w = None
      popn = {}
      if args.pop is not None:
        popi = [None] * (len(lst) - 9)
        sample_w = [1] * (len(lst) - 9)
        for i in range(9, len(lst)):
          s = lst[i]
          i2 = i - 9
          if s in pops:
            popi[i2] = pops[s]
            for p in set(pops[s]):
              if p in popskip: sample_w[i2] = 0
              if p not in popn: popn[p] = 0
              popn[p] += 1
          else: sample_w[i2] = 0
        hapunit.BlockData.sample_w = hapunit.Var.sample_w = sample_w
      break #continue

  if args.numProc >= 2:
    #from multiprocessing import Pool
    #pool = Pool(processes = args.numProc - 1)
    var_iter = mp_iter(infile, getvar, args, min(2, args.numProc-1))
  else:
    var_iter = sp_iter(infile, getvar, args)

  for var in var_iter:
    if chr != var.chr: # lst[0]:
      #if args.chr is not None: continue ##
      if chr != '':
        clear = True
      chr = var.chr # lst[0]
    if chr not in filt:
      chr1 = hapunit.changechr(chr)
      if chr1 in filt: filt[chr] = filt[chr1]
      else: filt[chr] = {}
    if str(var.pos) in filt[chr]: continue
    #if len(lst[3]) > 1 or len(lst[4]) > 1: continue

    #var = hapunit.Var(lst)
    #if var.maf() < args.minaf: continue
    if len(bdata) > 0 and var.pos - bdata.data[-1].pos > maxdis:
      clear = True

    if clear:
      out = bdata.find_blocks(clear = True) #, popi=popi, popskip=popskip)
      for o in out:
        outfile.write(io.tabjoin(o, '\n'))
      selim = selbase
      clear = False

    bdata.add_var(var)

    if len(bdata.select) >= selim:
      out = bdata.find_blocks(clear = False) #, popi=popi, popskip=popskip)
      for o in out:
        outfile.write(io.tabjoin(o, '\n'))

      selim = selbase
      ls = len(bdata.select)
      while selim - ls < keepflank: selim += selstep
    #if var.pos > 1000000: break


  print('last ls', len(bdata.select), 'la', len(bdata))
  out = bdata.find_blocks(clear = True) #, popi=popi, popskip=popskip)
  for o in out:
    outfile.write(io.tabjoin(o, '\n'))


if __name__ == '__main__':
  import sys, argparse
  p = argparse.ArgumentParser()
  set_parser(p)
  if len(sys.argv)==1:
    print(p.print_help())
    exit(0)
  run(p.parse_args())
