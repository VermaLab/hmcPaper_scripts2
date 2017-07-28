import numpy as np
import matplotlib.pyplot as plt
import random
import scipy.stats


ls


fols = [
    "jd_mc_lt50_jdAtac",
    "jd_mc_gt50_jdAtac",
    "hpne_mc_lt50_hpneAtac",
    "hpne_mc_gt50_hpneAtac",
    "jd_mc_lt50_jdRna",
    "jd_mc_gt50_jdRna",
    "hpne_mc_lt50_hpneRna",
    "hpne_mc_gt50_hpneRna"
]


f = fols[0]
f


for f in fols:
    plt.close()
    outFile = "pics4/peakSums_"+f+".png"
    x = np.load(f+"/peakSums.npy")
    plt.plot(np.arange(-3000, 3000, 1), x)
    plt.title(f)
    plt.savefig(outFile)




#layer the ATAC HMC and MC profiles on top of each other
#normalized by total count


x_mc = np.load("jd_mc_gt50_jdAtac/peakSums.npy")
x_mc = x_mc/np.mean(x_mc)
x_hmc = np.load("jd_hmc_p05_jdAtac/binnedPP.npy")
x_hmc = np.sum(x_hmc, 0)
x_hmc = x_hmc/np.mean(x_hmc)

plt.close()
plt.plot(x_mc/np.mean(x_mc), color="blue", label=">50% mc sites")
plt.plot(x_hmc/np.mean(x_hmc), color="red", label = "0.05 hmc sites")
plt.title("JD: coloc OxBS, atac-seq")
plt.legend()
plt.savefig("pics4/jd_oxbs_atac.png")


x_mc = np.load("hpne_mc_gt50_hpneAtac/peakSums.npy")
x_mc = x_mc/np.mean(x_mc)
x_hmc = np.load("hpne_hmc_p05_hpneAtac/binnedPP.npy")
x_hmc = np.sum(x_hmc, 0)
x_hmc = x_hmc/np.mean(x_hmc)

plt.close()
plt.plot(x_mc/np.mean(x_mc), color="blue", label=">50% mc sites")
plt.plot(x_hmc/np.mean(x_hmc), color="red", label = "0.05 hmc sites")
plt.title("HPNE: coloc OxBS, atac-seq")
plt.legend()
plt.savefig("pics4/hpne_oxbs_atac.png")


#devolop statistical test that shows whether two profiles have similar peak
#tighter or looser peak
#
#sample from distributions
#Do F test to show variances are different


#treat the profile as a probability density distribution
#sample n values from the distribution
def sampleFromProfile(pos, prof, n=1000):
    #step one, contruct prob dist
    csum = np.cumsum(prof)
    #normalize so last element is prob 1.0
    prob = csum/csum[-1]
    #randomly sample 1000 from uniform dist
    ix = [min(np.where(prob > random.random())[0]) for i in np.arange(n)]
    samps = pos[ix]
    return samps


def kurtosisPerm(samps_mc, samps_hmc, K=10000):
    #get stat for difference in kurtosis
    kmc = scipy.stats.kurtosis(samps_mc)
    khmc = scipy.stats.kurtosis(samps_hmc)
    kstat = khmc - kmc
    #
    #combine the samples
    both = np.array([samps_mc, samps_hmc]).ravel()
    res = []
    for i in range(K):
        if i % 100 == 0:
            print i
        #shuffle the samples
        np.random.shuffle(both)
        #compute kurtosis on the two halves of the array
        k1 = scipy.stats.kurtosis(both[:len(samps_mc)])
        k2 = scipy.stats.kurtosis(both[len(samps_hmc):])
        #save the stat of difference of kurtosis
        res.append(k2-k1)
    #
    res = np.array(res) 
    hits = np.sum(abs(kstat ) < np.abs(res))
    pval = hits / res.shape[0]
    plt.figure()
    plt.hist(res, 100)
    plt.axvline(x=kstat, color="red")
    plt.xlabel("Null Distribution (K=%d)"%K)
    return [kmc, khmc, pval]



mcFile = "jd_mc_gt50_jdAtac/peakSums.npy"
hmcFile = "jd_hmc_p05_jdAtac/binnedPP.npy"

#is the hmc profile tighter than mc?
def testTightness(mcFile, hmcFile):
    #prep the profiles, divide by mean so they on same scale
    x_mc = np.load(mcFile)
    x_mc_sum = np.sum(x_mc)
    x_mc = x_mc/np.mean(x_mc)
    x_hmc = np.load(hmcFile)
    #hmc file is the binned collection, not the sum profile
    x_hmc_sum = np.sum(x_hmc)
    x_hmc = np.sum(x_hmc, 0)
    x_hmc = x_hmc/np.mean(x_hmc)
    #get the samples
    pos = np.arange(-3000,3000,1)
    samps_mc = sampleFromProfile(pos, x_mc, 10000)
    samps_hmc = sampleFromProfile(pos, x_hmc, 10000)
    #plot the figures to make sure data makes sense
    #the orginal read sum profile
    plt.figure()
    plt.plot(pos, x_mc, label="mc", color="red")
    ax = plt.plot(pos, x_hmc, label="hmc", color="blue")
    plt.ylim(ymin = 0)
    plt.legend()
    #
    #histogram of the sampled read sum profile
    plt.figure()
    plt.hist(samps_mc, 100, alpha=0.5, label="mc")
    plt.hist(samps_hmc, 100, alpha=0.5, label="hmc")
    plt.legend()
    #
    #kurtosis test
    return kurtosisPerm(samps_mc, samps_hmc)
    #
    #variance test
    #mc_var = np.var(samps_mc)
    #hmc_var = np.var(samps_hmc)
    ##perform the levene's test on variance to see which
    ##has less (tighter), and which has more (looser)
    #res = scipy.stats.levene(samps_hmc, samps_mc)
    #return [mc_var, hmc_var, res.pvalue]



jd_tight = testTightness("jd_mc_gt50_jdAtac/peakSums.npy", "jd_hmc_p05_jdAtac/binnedPP.npy")
hpne_tight = testTightness("hpne_mc_gt50_hpneAtac/peakSums.npy", "hpne_hmc_p05_hpneAtac/binnedPP.npy")





jd_tight  
hpne_tight 


kmc = scipy.stats.kurtosis(samps_mc)
khmc = scipy.stats.kurtosis(samps_hmc)
abs(kmc ) - abs(khmc)

both = np.array([samps_mc, samps_hmc]).ravel()
res = []
for i in range(10000):
    if i % 100 == 0:
        print i
    np.random.shuffle(both)
    k1 = scipy.stats.kurtosis(both[:10000])
    k2 = scipy.stats.kurtosis(both[10000:])
    res.append(abs(k1 ) - abs(k2))

plt.hist(res, 100)

#F = min(np.var(samps_hmc) / np.var(samps_mc), np.var(samps_mc) / np.var(samps_hmc) )
#2*scipy.stats.f.cdf(F, len(samps_mc)-1, len(samps_hmc)-1)
#p_value = scipy.stats.f.cdf(F, len(samps_mc)-1, len(samps_hmc)-1)

