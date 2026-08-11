__all__=[]

commands = {
  'block': 'Haplotype unit block identification',
  'group': 'Haplotype unit group classification',
  'convert': 'Convert haplotype unit group vcf to HUG genotype vcf',
  'plot': 'Plot haplotypes',
  }

def load(cmd):
  if cmd == 'block':
    from . import block
    return block
  elif cmd == 'group':
    from . import group
    return group
  elif cmd == 'convert':
    from . import convert
    return convert
  elif cmd == 'plot':
    from . import plot
    return plot
  else:
    return None
