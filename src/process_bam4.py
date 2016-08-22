import pickle
import numpy as np
import pdb
from scipy.stats import lognorm
import scipy
import matplotlib.pyplot as plt
pdb.set_trace()


def read_pickle(filename):
    with open(filename, 'rb') as handle:
        molecule = pickle.load(handle)

    all_frag_len = []
    all_reads_start_diff = [ ]
 
    curr = 0

    for MI_num,start_pos_diff in molecule.iteritems():
        print curr
        curr = curr + 1
        all_reads_start_diff.extend(start_pos_diff)
        frag_len = start_pos_diff[-1] - start_pos_diff[0] + 150
        all_frag_len.append(frag_len)



    samples =  [float(i)/1000  for i in all_frag_len if i < 200000]   # 10000
    #samples =   [float(i)/1000 for i in all_reads_start_diff]
    shape, loc, scale = scipy.stats.lognorm.fit(samples, floc=0 )
    x_fit = np.linspace(np.min(samples), np.max(samples), 1000 )
    samples_fit = scipy.stats.lognorm.pdf( x_fit, shape, loc=loc, scale=scale )
    plt.subplot( 1, 2, 1 )
    N_bins = 1000
    counts, bin_edges, ignored = plt.hist(samples, N_bins, histtype='stepfilled', label='histogram' )
    area_hist = .0
    for ii in range(counts.size):
    	area_hist += (bin_edges[ii+1]-bin_edges[ii]) * counts[ii]
    # plot fit into histogram
    plt.plot( x_fit, samples_fit*area_hist, label='fitted and area-scaled PDF', linewidth=2)
    plt.legend()
    plt.show()


    """
    numBins=1000
    fig=plt.figure()
    ax=fig.add_subplot(211)
    ax.hist(all_frag_len,numBins,color='green')
    #plt.xlim(0, 300000)

    ax2=fig.add_subplot(212)
    ax2.hist(all_reads_start_diff,numBins,color='green')
    #plt.xlim(0, 300000)
    plt.show()
    return all_frag_len

    """


def read_bam(filename):
    molecule = {}
    first = {}
    curr = 0
    with open(filename,'r') as f:
        for line in f:
            print curr
            curr = curr + 1
            data = line.split('\t')
            start_pos = int(data[3])
            MI = int(data[-1].rstrip().split(':')[2])
            if molecule.has_key(MI):
                molecule[MI].append(start_pos -first[MI])
            else:
                first[MI] = start_pos
                molecule[MI] = [start_pos - first[MI]]
    
    all_frag_len = []
    all_reads_start_diff = [ ]
    curr = 0
    with open('molecule3.pickle', 'wb') as handle:
        pickle.dump(molecule, handle)

    for MI_num,start_pos_diff in molecule.iteritems():
        print curr
        curr = curr + 1
        all_reads_start_diff.extend(start_pos_diff)
        frag_len = start_pos_diff[-1] - start_pos_diff[0] + 150
        all_frag_len.append(frag_len)


    with open('all_reads_start_diff3.pickle', 'wb') as handle:
        pickle.dump(all_reads_start_diff, handle)
   
    with open('all_frag_len3.pickle', 'wb') as handle:
        pickle.dump(all_frag_len, handle)

    return all_frag_len, molecule,all_reads_start_diff



#read_bam("NA12878.chr19_allMI.sam")
#read_bam("test.sam")
molecule = read_pickle('../molecule3.pickle')
