import sys
import numpy as np

#input_file = sys.argv[1]
#outfile = sys.argv[2]

#threshold
single_line_average = 0.9
block_average = 0.95
block_minimal_size = 3
gap = 2

def sum_line(arr, start, stop, i, skip = {}, remove_diagonal = True):
  s = sum(arr[start:stop,i])
  n = stop - start
  if remove_diagonal:
    if start <= i < stop:
      s -= 1
      n -= 1
  for k in skip:
    if remove_diagonal and k == i: continue
    if start <= k < stop:
      s -= arr[k, i]
      n -= 1
  return s, n

def average_line(arr, start, stop, i, skip = {}, remove_diagonal = True):
  s, n = sum_line(arr, start, stop, i, skip, remove_diagonal)
  return s / n

def average_block(arr, start, stop, skip = {}, remove_diagonal = True):
  s, n = 0, 0
  for i in range(start, stop):
    if i in skip: continue
    sl, nl = sum_line(arr, start, stop, i, skip, remove_diagonal)
    s += sl
    n += nl
  return s / n


def sum_sel(arr, i, sel, skip = {}, remove_diagonal = True):
  s, n = 0, 0
  for j in sel:
    if j in skip: continue
    if remove_diagonal and j == i: continue
    #if arr[j, i] is None: continue
    s += arr[j, i]
    n += 1
  return s, n

def average_sel(arr, i, sel, skip = {}, remove_diagonal = True):
  s, n = sum_sel(arr, i, sel, skip, remove_diagonal)
  return s / n

def remove_severals(cluster_list, gap):
  if len(cluster_list)<2:
    return cluster_list
  c = 0
  while c < len(cluster_list) - 1:
    if cluster_list[c+1][0] - cluster_list[c][1] > gap:
      c += 1
      continue
        #if cluster_list[c+1][0]-cluster_list[c][1] <= gap:
    remove_points = cluster_list[c+1][0] - cluster_list[c][1]
    start = cluster_list[c][0]
    end = cluster_list[c+1][1]
    additional_remove_sum = 0
    ignore_points_list = []
    if len(remove_points_dict) > 0: # ignore_points_dict != {}:
                # add non_removal_cluster to dict
                #for item in cleaned_cluster:
                #    tc = tuple(item)
                #    if tc not in ignore_points_dict:
                #        ignore_points_dict[tc] = []
      tc = tuple(cluster_list[c])
      if tc in remove_points_dict: ignore_points_list += remove_points_dict[tc]
      tc1 = tuple(cluster_list[c+1])
      if tc1 in remove_points_dict: ignore_points_list += remove_points_dict[tc1]
                #ignore_points_list = ignore_points_dict[tuple(cluster_list[c])]+ignore_points_dict[tuple(cluster_list[c+1])]
            #for i in ignore_points_list:
                #additional_remove_sum += (sum(arr[i,start:cluster_list[c][1]])+sum(arr[i,cluster_list[c+1][0]:end])-1)*2
    if cluster_list[c][1] != cluster_list[c+1][0]:
      ignore_points_list += list(range(cluster_list[c][1], cluster_list[c+1][0]))
    sa, na = 0, 0
    for i in range(cluster_list[c+1][0], cluster_list[c+1][1]):
      if i in ignore_points_list: continue
      s, n = sum_line(arr, cluster_list[c][0], cluster_list[c][1], i, skip = ignore_points_list)
      sa += s
      na += n
    ba1 = sa / na
    ba2 = average_block(arr, start, end, skip = ignore_points_list)
    print('r1 c: {}, c[c]: {}, c[c+1]: {}, ignore_list: {}, ba1: {}, ba2: {}'.format(c, cluster_list[c], cluster_list[c+1], ignore_points_list, ba1, ba2))
    if ba1 > block_average:
                #if (sum(sum(arr[start:end,start:end]))-additional_remove_sum-(end-start)-sum(sum(arr[cluster_list[c][1]:cluster_list[c+1][0],start:end]))*2+ sum(sum(arr[cluster_list[c][1]:cluster_list[c+1][0],cluster_list[c][1]:cluster_list[c+1][0]]))+remove_points)*1.0/((end-start-1)*(end-start)-len(ignore_points_list)*(end-start-1)*2-remove_points*(end-start)*2+remove_points*remove_points+remove_points)>block_average:
                    #if cluster_list[c+1][0] > cluster_list[c][1]:
                        #if (start,end) in remove_points_dict:
                            #remove_points_dict[(start,end)] += [m for m in range(cluster_list[c][1],cluster_list[c+1][0])]
                        #else:
      remove_points_dict[(start,end)] = ignore_points_list # [m for m in range(cluster_list[c][1],cluster_list[c+1][0])]
                    
                    #print("yes",str(cluster_list[c][1]),'-',str(cluster_list[c+1][0]))
                   # print(remove_points_dict)
      cluster_list[c] = [start, end]
      print('r1c cluster_list[c] = {}, ignore_points_list = {}, del cluster_list[c+1] {}'.format([start, end], ignore_points_list, cluster_list[c+1]))
      del cluster_list[c+1]
    #  continue
    else:
      c += 1
        #if c == len(cluster_list)-1:
        #    break
  return cluster_list


def go_up_tuning(clusters, m, n):
  if len(clusters) > 0: last = clusters[-1][0]
  else: last = 0
  for v in range(m-1, last-1, -1):
    av2 = average_line(arr, v, n, v)
    ba2 = average_block(arr, v, n)
    print('g1 v: {}, m: {}, n: {}, av2: {}, ba2: {}'.format(v, m, n, av2, ba2))
    if len(clusters) == 0 or v >= clusters[-1][1] or v - clusters[-1][0] < 2:
      if av2 < single_line_average or ba2 < block_average:
        if len(clusters) == 0 or v+1 >= clusters[-1][1]:
          clusters.append([v+1, n])
          print('g11, clusters append', [v+1, n])
        else:
          clusters[-1]=[v+1, n]
          print('g12, clusters[-1]=', [v+1, n])
        return True
    else:
      av0 = average_line(arr, clusters[-1][0], v+1, v)
      print('g2 av0: {}'.format(av0))
      if av2 < av0 or ba2 < block_average:
        if v+1 >= clusters[-1][1]:
          clusters.append([v+1, n])
          print('g21, clusters append', [v+1, n])
        else:
          clusters[-1] = [clusters[-1][0], v+1]
          clusters.append([v+1, n])
          print('g22, clusters[-1]=', clusters[-2], ', clusters.append', [v+1, n], ', clusters is', clusters)
        return True
  else:
    if len(clusters) > 0:
      clusters[-1] = [last, n]
      print('g31, go_up_tuning False, clusters[-1] = ', [last, n])
    else:
      clusters.append([last, n])
      print('g32, go_up_tuning False, clusters.append ', [last, n])
    return False

    #if clusters[-1][1] - clusters[-1][0] < 3:
    #    for v in range(clusters[-1][1]-1, clusters[-1][0]-1, -1): ## not including the first variant? Changed
    #        #av = (sum(arr[v:n-1,v:n-1][0])-1)*1.0/(n-v-2)
    #        av2 = average_line(arr, v, n, v)
    #        ba2 = average_block(arr, v, n)
    #        print('g1 m: {}, n: {}, v: {}, av2: {}, ba2: {}'.format(m, n, v, av2, ba2))
    #        if average_line(arr, v, n, v) < single_line_average or average_block(arr, v, n) < block_average: # (sum(arr[v:n-1,v:n-1][0])-1)*1.0/(n-v-2)<single_line_average:
    #        #ready to break
    #            clusters[-1]=[v+1, n] ## override the previous one?
    #            print('g1, clusters[-1]=', [v+1, n])
    #            return True
    #else:
    #    for v in range(clusters[-1][1]-1, clusters[-1][0], -1): ## Changed
    #        #if arr[v:n-1,v:n-1][0] is reder, it will be included in the next block,else break
    #        if average_line(arr, v, n, v) < average_line(arr, clusters[-1][0], v+1, v) or average_block(arr, v, n) < block_average:
    #        #if (sum(arr[v:n-1,v:n-1][0]) - 1) * 1.0 / (n-2-v) < (sum(arr[clusters[-1][0]:(v+1),clusters[-1][0]:(v+1)][-1]) - 1) * 1.0 / (v - clusters[-1][0]):
    #            #ready to break
    #            if v < clusters[-1][1] - 1:
    #                clusters[-1] = [clusters[-1][0], v+1]
    #                clusters.append([v+1, n])
    #                print('g2, clusters[-1]=', '[clusters[-1][0],v+1]', ', clusters.append', [v+1, n], ', clusters is', clusters)
    #            else:
    #                clusters.append([m, n])
    #                print('g3, clusters.append', [m, n])
    #            return True
    #print('go_up_tuning False')
    #clusters[-1] = [clusters[-1][0], n] ## move inside the function
    #return False

def remove_tiny(cluster_list):
    cluster_1=[]                            
    for item in cluster_list:
        #m=item[0]
        #n=item[1]+1
        if item[1] - item[0] >= block_minimal_size:
            cluster_1.append(item)
    #m=cluster_list[-1][0]
    #n=cluster_list[-1][1]
    #if n-m>=block_minimal_size:       
    #    cluster_1.append(cluster_list[-1])
    return cluster_1
        

def isolate_and_one_up(cluster_list, gap=2):
    for c,i in enumerate(cluster_list):
        if c==0: last = 0
        elif cluster_list[c-1][1] - cluster_list[c-1][0] >= 3: last = cluster_list[c-1][1]
        else: last = cluster_list[c-1][0]
        break_flag = False
        #print('ii c==0', c, i, cluster_list)
        if i[0] - last < gap: continue
        skip = []
        if tuple(i) in remove_points_dict: skip += remove_points_dict[tuple(i)]
        skip_pos = 1
                #may go up gap step
        for step in range(1,gap+1):
            break_flag = False
            if i[0] - step <= last:
                break_flag = True
                print('ii1 i[0]-step-1<0, step', step)
            else:
                av2 = average_line(arr, i[0]-step-1, i[1], i[0]-step-1, skip = skip + list(range(i[0]-skip_pos, i[0])))
                ba2 = average_block(arr, i[0]-step-1, i[1], skip = skip + list(range(i[0]-skip_pos, i[0])))
                print('i1 c: {}, i: {}, step: {}, av2: {}, ba2: {}'.format(c, i, step, av2, ba2))
                if av2 < block_average or ba2 < block_average: # single_line_average
            #if (sum(sum(arr[(i[0]-step-1):i[1],(i[0]-step-1):i[1]]))-(i[1]-i[0]+step+1)-sum(sum(arr[(i[0]-step):i[0],(i[0]-step-1):i[1]]))*2+ sum(sum(arr[(i[0]-step):i[0],(i[0]-step):i[0]]))+step)*1.0/((i[1]-i[0]+step)*(i[1]-i[0]+step+1)-step*(i[1]-i[0]+step+1)*2+step*step+step)<block_average: 
                    break_flag = True
            if break_flag:
                if step > skip_pos:
                    #if (i[0]-step,i[1]) in remove_points_dict:
                        #remove_points_dict[(i[0]-step, i[1])] += skip + list(range(i[0]-skip_pos, i[0])) # changed from range(i[0]-step+1,i[0])
                    #else:
                    remove_points_dict[(i[0]-step,i[1])] = skip + list(range(i[0]-skip_pos, i[0])) # changed from range(i[0]-step+1,i[0])
                            #print(step-1,"yes2",i[0]-step+1,"-",i[0])
                    print('i1, cluster_list[{}]= {}'.format(c, str([i[0]-step,i[1]])))
                    cluster_list[c]= [i[0]-step,i[1]]
                    break
                else:
                  skip_pos = step + 1
                  #break_flag = False
                  #continue
                #break
            else:
                print('ii1 >=block_average', step)
        if not break_flag and skip_pos <= gap: # not break_flag:
            #if (i[0]-gap-1,i[1]) in remove_points_dict:
            #    remove_points_dict[(i[0]-gap-1,i[1])] += [i[0]-1] ## changed from range(i[0]-gap,i[0])
            #else:
            remove_points_dict[(i[0]-gap-1,i[1])] = skip + list(range(i[0]-skip_pos, i[0])) ## [m for m in range(i[0]-1 ,i[0])]
                    #print(gap,"yes2",i[0]-step,"-",i[0])
            print('i2, cluster_list[{}]= {}'.format(c, str([i[0]-gap-1,i[1]])))
            cluster_list[c]= [i[0]-gap-1,i[1]]

    return cluster_list
                    
                
                               
def settle_remove_points(a_remove_points_dict, cleaned_cluster):
    #final_cluster=[]
    #pairs = list(a_remove_points_dict.keys())
    #for m in range(len(pairs)):
    #    break_flag = False
    #    for n in range(len(pairs)):
    #        if pairs[m][0]>=pairs[n][0] and pairs[m][1]<=pairs[n][1] and pairs[m]!=pairs[n]:
    #            break_flag = True
    #            break
    #    if not break_flag:
    #        final_cluster.append(pairs[m])
            
    final_remove_points_dict={}
    for i in cleaned_cluster: # final_cluster:
        #final_remove_points = []
        #for item in a_remove_points_dict.keys():
        #    if item[0]>=i[0] and item[1]<=i[1]:
        #        final_remove_points += a_remove_points_dict[tuple(item)]
      ti = tuple(i)
      if ti in a_remove_points_dict:
        final_remove_points_dict[ti] = a_remove_points_dict[ti] # final_remove_points
    return final_remove_points_dict 


def correct_edge(cluster_list):
    for c,i in enumerate(cluster_list):
        for m in range(i[0],i[1]-1):
            if average_line(arr, m, i[1], m) >= single_line_average: #((sum(arr[m:i[1],m:i[1]][0]) - 1) / (i[1]-m-1)) > single_line_average:
                if m > i[0]: cluster_list[c]=[m,i[1]]
                break
    return cluster_list

def get_blocks(corr):
  global arr, remove_points_dict # , cleaned_cluster
  arr  = np.array(corr)
  print('arr', len(arr))
  remove_points_dict = {}
  m=0
  clusters = []
  while m < len(arr) - 1:
    #cloud get the last+1
    for n in range(m+2, len(arr)+1):
        #not satisfied,could reach nth+1,get the last
      #av = (sum(arr[m:n,m:n][-1]) - 1) / (n-m-1)
      #av2 = average_line(arr, m, n, n-1)
      #ba = (sum(sum(arr[m:n,m:n])) - (n-m)) / (n-m-1) / (n-m)
      #ba2 = average_block(arr, m, n)
      av2 = average_line(arr, m, n, n-1)
      ba2 = average_block(arr, m, n)
      print('m: {}, n: {}, av: {}, av2: {}, ba: {}, ba2: {}'.format(m, n, 'NA', av2, 'NA', ba2))
      if av2 < single_line_average or ba2 < block_average: # (sum(arr[m:n,m:n][-1]) - 1) < (n-m-1) * single_line_average:
        #if average_block(arr, m, n) < block_average: # (sum(sum(arr[m:n,m:n])) - (n-m)) < (n-m-1) * (n-m) * block_average:
        if n - m <= block_minimal_size - 1:
          m += 1
          continue
        #if len(clusters) > 0 and m == clusters[-1][-1]:
          #if m == clusters[-1][-1]:
                    #go up,act on clusters
        break_flag = go_up_tuning(clusters, m, n-1)
          #if break_flag:
            #m = clusters[-1][-1]
            #break
          #else:
                        # the loop is over if could be concatenate,do ,else, from v+1
            #clusters[-1] = [clusters[-1][0], n-1]
            #print('m1, clusters[-1] = ', '[clusters[-1][0], n-1]', ', clusters is', clusters)
                  #if not go up        
        #else:
          #clusters.append([m, n-1])
          #print('m2, clusters.append', [m, n-1])
              #if not go up       
        #else:
          #clusters.append([m,n-1])
          #print('m3, clusters.append', [m,n-1], ', clusters is', clusters)
        m=clusters[-1][-1]
        break
 #   print(clusters,m,n)
    if m >= len(arr) - 1:
      break
    #statisfied the contdition until the last,n included
    if n >= len(arr):
      #if len(clusters) > 0 and m == clusters[-1][-1]:
      break_flag = go_up_tuning(clusters, m, n)
        #if m == clusters[-1][-1]:
                #go up 
      #else:
        #clusters.append([m,n])
        #print('m10, clusters.append', [m, n])
      break

  if len(clusters) == 0: return [], {} # cleaned_cluster, remove_points_dict
 #   cleaned_cluster =remove_tiny(correct_edge(concatenate(remove_tiny(clusters))))
  previous_cluster = remove_tiny(isolate_and_one_up(remove_severals(correct_edge(clusters), gap), gap))
  if len(previous_cluster) == 0: return [], {}
    #if len(previous_cluster) > 0:
  cleaned_cluster = []
  while cleaned_cluster != previous_cluster: # True:
    cleaned_cluster = remove_tiny(isolate_and_one_up(remove_severals(previous_cluster, gap), gap))
    if cleaned_cluster == previous_cluster:
      break
    else:
      previous_cluster = cleaned_cluster

  remove_points_dict = settle_remove_points(remove_points_dict, cleaned_cluster)
  #last_cleaned_cluster = remove_severals(cleaned_cluster, gap)
  #remove_points_dict = settle_remove_points(remove_points_dict, cleaned_cluster)
  #while True:
    #cleaned_cluster =  remove_severals(last_cleaned_cluster, gap)
    #remove_points_dict = settle_remove_points(remove_points_dict, cleaned_cluster)
    #if cleaned_cluster == last_cleaned_cluster:
    #  break
    #else:
    #  last_cleaned_cluster=cleaned_cluster

    #else:
    #  cleaned_cluster = []
    #  remove_points_dict =  {}

  #else:
  #  cleaned_cluster = []
  #  remove_points_dict= {}
        
  return cleaned_cluster, remove_points_dict 

#print(remove_points_dict)


