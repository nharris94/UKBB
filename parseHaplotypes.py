import sys, re, math, random

seqLen=0

def parse_seqs():
    global seqLen
    seq = {}  # seq[readname][2]='A'
    for read in sys.stdin:
        a = read.strip().split('\t')
        if a[0] == '@SQ':
            if seqLen > 0:
                print 'Only one contig allowed'
                sys.exit(1)
            seqLen=int(a[2][3:])
            continue
        if a[0][0] == '@': continue
            
        readLen=len(a[9])
        readName=a[0]
        flags = int(a[1])
        b=[0]*12
        for i in range(12):
            bit=flags%2
            b[i]=bit
            flags/=2
            if b[8]:
                # secondary
                continue
            if b[9]:
                # failure
                continue
            if b[11]:
                # supplemental
                continue
            if b[2]:
                # unmapped
                continue
            if b[10]:
                # duplicate
                continue
    #if readName in seq: readName+='.2'
        if readName not in seq: seq[readName]=['X']*seqLen
    # else read is other mate in pair, treat as one joint sequence
  
        #### CIGAR operators that consume sequence
        # M alignment match (can be a sequence match or mismatch)
        # = sequence match
        # X sequence mismatch
        # I insertion to the reference
        # S soft clipping (clipped sequences present in SEQ)
        #### CIGAR operators that do NOT consume sequence
        # N skipped region from the reference
        # D deletion from the reference
        # H hard clipping (clipped sequences NOT present in SEQ)
        # P padding (silent deletion from padded reference)

        cigar = ['M']*readLen
        cig = re.findall("([0-9]+)([MISX=DN])",a[5]) # just the symbols that consume sequence, plus D to count indels
        pos = 0
        rpos=int(a[3])-1 # make 0-based
        seq1 = a[9]
        insLen = 9999
        if a[6] == '=' and a[8] != '0':insLen = abs(int(a[8]))
    #print read.strip()
        for c in cig:
            if pos >insLen/2: break # only look up to the midpoint for overlapping reads...
            c0 = int(c[0])
            if c[1]=='M':
                for i in range(c0):
                    seq[readName][rpos]=seq1[pos]
                    pos+=1
                    rpos+=1
            elif c[1]=='S':
                pos+=c0
            elif c[1]=='I':
                seq[readName][rpos]='I'
                pos+=c0
            elif c[1]=='D' or c[1]=='N':
                for i in range(c0):
                    seq[readName][rpos]='D'
                    rpos+=1
            else:
                print 'problem...',c
    return seq


def compare_seqs(seq):
    info=[0]*seqLen
    select = [] # positions worth considering
    poly = {} # poly[pos][base]=freq if freq>.05
    for i in range(seqLen):
        c={}
        d={}
        n=0
        info1=0
        for r in seq:
            b=seq[r][i]
            if b == 'X': continue
            if b not in c:c[b]=1
            else: c[b]+=1
            n+=1
        for x in c:
            p = 1.0*c[x]/n
            if p<1: info1+= (-p*math.log(p,2))
            if p>.05:
                if i not in poly: poly[i]={}
                poly[i][x]=p
    #if info1>0.3: print ' info',i,info1
        if info1>0.5:
            info[i]=1
            select.append(i)

#print select
    sim = {}
    dis = {}
    for a in seq:
        for b in seq:
            sim[(a,b)]=0
            dis[(a,b)]=0
            score = 0
            for i in select: #range(max(beg[a],beg[b]),min(end[a],end[b])):
                sa=seq[a][i]
                if sa == 'X': continue
                if sa not in poly[i]: continue # also ignores errors at informative positions
                sb=seq[b][i]
                if sb == 'X': continue
                if sb not in poly[i]: continue
                if sa == sb: sim[(a,b)]+=(-math.log(poly[i][sa],2))  #+=info[i]
                if sa != sb: dis[(a,b)]+=info[i]

    return (sim,dis)


def print_cluster(seq,clust,ccnt,cmem):
    num = 1
    for c in sorted(ccnt,key=lambda x:ccnt[x],reverse=True):
        #print ccnt[c],'ccnt'
        #for x in cmem[c]:
            #print ''.join(seq[x][0:160]),'cmem',x,clust[x]
        if ccnt[c]>0: 
            (ret,cov,covX)=get_cons(seq,clust,c)
            print '>cluster'+str(num)+'_'+str(ccnt[c])+'_'+str(covX),cov
            print ret
            num+=1


def print_cluster_seqs(seq,clust,ccnt,cmem):
    num = 1
    for c in sorted(ccnt,key=lambda x:ccnt[x],reverse=True):
        #print ccnt[c],'ccnt'
        for x in cmem[c]:
            print x+'\t'+str(num)+'\t'+''.join(seq[x])
        num+=1

# write consensus
def get_cons(seq,clust,clustnum):
    ret=''
    cov=''
    covX=0
    for i in range(seqLen):        
        c={}
        d={}
        n=0
        info1=0
        for r in seq:
            if clust[r] != clustnum: continue
            #if i==1: print "".join(seq[r][0:160]),"elem",clust[r],clustnum
            b=seq[r][i]  
            if b == 'X': continue
            if b not in c:c[b]=1
            else: c[b]+=1
            n+=1
        cov10=n/10
        if cov10 > 9: cov10 = 9
        cov+=str(cov10)
        #if i>260 and i<360: covX+=n
        covX+=n
        nt='X'
        ntCnt=0
        for x in c:
            if c[x]>ntCnt:
                nt=x
                ntCnt=c[x]
        if nt == 'D': pass  # deletion
        else: ret+=nt
    return (ret,cov,covX/seqLen) #covX/100)

def init_clusters():
    clust = {}
    ccnt = {}
    cmem = {}

# initialize each to its own cluster, or with its pair if it has one
# revised code above to group the pairs more directly
# still works as is, just no sequences with .2 in names
    n = 0
    for a in sorted(seq):
        if a[-2:] == '.2':
            # add to existing mate's cluster
            mate = a[:-2]
            clust[a]=clust[mate]
            ccnt[clust[a]]=2
            cmem[clust[a]].append(a)
        else:
            # create new cluster
            clust[a]=n
            ccnt[clust[a]]=1
            cmem[clust[a]]=[]
            cmem[clust[a]].append(a)
            n+=1
            
    return (clust,ccnt,cmem)


# sort by similarity.. join most similar first, as long as they match perfectly
def cluster(clust,ccnt,cmem,sim,dis,ssmin):
    nclust = 0
    for sa,sb in sorted(sim,key=lambda x:sim[x],reverse=True):
        if sim[sa,sb] == 0: break
        a=clust[sa]
        b=clust[sb]
        if a == b: continue # already clustered
        if ccnt[a] == 0 or ccnt[b] == 0: continue #empty cluster
        if ccnt[b] > ccnt[a]: continue # merge to smaller cluster?
        dd=0 # check for any dissimilarity
        ss=0 # overall similarity score
        nn=0 # number of reads involved in comparison
        for a1 in cmem[a]:
            for b1 in cmem[b]:
                nn+=1
                dd=dis[a1,b1]
                ss+=sim[a1,b1]
                if dd > 0:  break # do not cluster
            if dd > 0: break # do not cluster
        if dd>0: continue # different, do not cluster
        if nn==0: continue # no overlap
        if ss/nn<ssmin: continue # not really similar
            # cluster b>>a
        tmpb = list(cmem[b])
        cmem[a]=cmem[a]+cmem[b]
        cmem[b]=[]
        ccnt[a]=ccnt[a]+ccnt[b]
        ccnt[b]=0
        for b1 in tmpb:
            clust[b1]=a
    return (clust,ccnt,cmem)


def sift_cluster(clust,ccnt,cmem,sim,dis,cut):
    # look at the members of each cluster, if > frac of the members could be reassigned to another cluster, do so
    for c in sorted(ccnt,key=lambda x:ccnt[x],reverse=False): # for each cluster
        if ccnt[c] == 0: continue # empty cluster
        nn = 0 # count the number of sequences in the cluster that could be reassigned
        sift = {} # remember potential assignments
        for s in cmem[c]: # for each sequence in the cluster            
            for c1 in sorted(ccnt,key=lambda x:ccnt[x],reverse=True): # test for compatibility with each other cluster
                if ccnt[c1] < 6: continue # tiny cluster
                #if ccnt[c1] < ccnt[c]/2: continue # >2x smaller cluster
                if c == c1: continue # same cluster
                dd = 0 # check for any dissimilarity
                for s1 in cmem[c1]:
                    dd=dis[s,s1]
                    if dd > 0:  break # no match
                if dd == 0: # match for seq s to cluster c1
                    nn+=1 
                    sift[s]=c1
                    break
            # close c1
        # close s

        #sys.stderr.write(str(c)+'\t'+str(ccnt[c])+'\t'+str(nn)+'\t'+str(sorted(sift.values()))+'\n' )
        if ccnt[c] - nn > cut: continue
        # reassignable -- go through the list again and modify
        for s in sift: # move sequence s to cluster sift[s]
            clust[s] = sift[s]
            cmem[sift[s]] = cmem[sift[s]] + [s]
            cmem[c].remove(s)
            ccnt[sift[s]]+=1
            ccnt[c]-=1

    return (clust,ccnt,cmem)

def balance_cluster(clust,ccnt,cmem,sim,dis):
    # randomly swap reads to compatible clusters
    sift = {}
    for c in sorted(ccnt,key=lambda x:ccnt[x],reverse=True): # for each cluster
        if ccnt[c] < 6: continue # tiny cluster
        for s in cmem[c]: # for each sequence in the cluster
            for c1 in sorted(ccnt,key=lambda x:ccnt[x],reverse=False): # test for compatibility with each other cluster
                if ccnt[c1] < 6: continue # tiny cluster
                if c == c1: continue # same cluster
                dd = 0 # check for any dissimilarity
                for s1 in cmem[c1]:
                    dd=dis[s,s1]
                    if dd > 0:  break # no match
                if dd == 0: # match for seq s to cluster c1
                    if random.random() > 0.5: sift[s]=c1
                    break
            # close c1
        # close s
    # close c
#        sys.stderr.write(str(c)+'\t'+str(ccnt[c])+'\t'+str(nn)+'\t'+str(sorted(sift.values()))+'\n' )
#        if 1.0*(nn+1)/ccnt[c] < frac: continue

    for s in sift: # move sequence s to cluster sift[s]
        c =  clust[s]
        clust[s] = sift[s]
        cmem[sift[s]] = cmem[sift[s]] + [s]
        cmem[c].remove(s)
        ccnt[sift[s]]+=1
        ccnt[c]-=1

    return (clust,ccnt,cmem)


if __name__ == '__main__':
    seq = parse_seqs()
    (sim,dis) = compare_seqs(seq)
    (clust,ccnt,cmem) = init_clusters()
    (clust,ccnt,cmem) = cluster(clust,ccnt,cmem,sim,dis,1)

    (clust,ccnt,cmem) = sift_cluster(clust,ccnt,cmem,sim,dis,5)
    #sys.stderr.write('\n')
    #(clust,ccnt,cmem) = sift_cluster(clust,ccnt,cmem,sim,dis,0)

    (clust,ccnt,cmem) = balance_cluster(clust,ccnt,cmem,sim,dis)
    (clust,ccnt,cmem) = balance_cluster(clust,ccnt,cmem,sim,dis)
    #sys.stderr.write('\n')
    #(clust,ccnt,cmem) = sift_cluster(clust,ccnt,cmem,sim,dis,0)


    print_cluster(seq,clust,ccnt,cmem)
#    print_cluster_seqs(seq,clust,ccnt,cmem)

#    for r in seq:
#        print r+'\t'+str(clust[r])+'\t'+"".join(seq[r])
