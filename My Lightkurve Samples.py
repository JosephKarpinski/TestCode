

# %%
print("\nSection 001\n")
import lightkurve as lk
import matplotlib.pyplot as plt
%matplotlib inline

from lightkurve import search_targetpixelfile
# %%

# %%
print("\nSection 002\n")
def FoldTarget(target, targetquarter, targetexp, targetperiod):
    pixelfile = search_targetpixelfile(target, quarter=targetquarter, exptime=targetexp).download();
    #print("\nDownloaded")

    pixelfile.plot(frame=1);

    lc = pixelfile.to_lightcurve(aperture_mask=pixelfile.pipeline_mask);
    lc.plot();
    plt.title("Light Curve")
    #print("\nMask Done!")

    """lc = pixelfile.to_lightcurve(aperture_mask='all');
    lc.plot();
    print("\nAll Done!")"""

    flat_lc = lc.flatten(window_length=401);
    flat_lc.plot();
    plt.title("Flat")
    #print("\nFlat Done!")

    folded_lc = flat_lc.fold(period=targetperiod);
    folded_lc.plot();
    plt.title("Folded")
    #print("\nFold Done!")

    binned_lc = folded_lc.bin(time_bin_size=0.01);
    binned_lc.plot();
    plt.title("Binned")
    #print("\nBinn Done!")

def ProcessTarget(target, targetquarter, targetexp):
    pixelfile = search_targetpixelfile(target, quarter=targetquarter, exptime=targetexp).download();
    print("\nDownloaded")

    pixelfile.plot(frame=1);

    lc = pixelfile.to_lightcurve(aperture_mask=pixelfile.pipeline_mask);
    lc.plot();
    print("\nMask Done!")

    """lc = pixelfile.to_lightcurve(aperture_mask='all');
    lc.plot();
    print("\nAll Done!")"""

    flat_lc = lc.flatten(window_length=401);
    flat_lc.plot();
    print("\nFlat Done!")
    
# %%

# %%
print("\nSection 003\n")
search_result = lk.search_lightcurve('KIC 8462852', author='Kepler')
print(search_result)

ProcessTarget("KIC 8462852", 16, 1800, )

print("\nThe plot reveals a short-lived 20% dip in the brightness of the star. It looks like we re-discovered one of the intriguing dips in Tabby's star.\n")
# %%




# %%
print("\nSection 004\n")
# https://docs.lightkurve.org/tutorials/1-getting-started/what-are-lightcurve-objects.html

search_result = lk.search_lightcurve('KIC 6922244', author='Kepler')
print(search_result)

FoldTarget("KIC 6922244", 4, 1800, 3.5225)

# %%


# https://arxiv.org/pdf/1712.05044.pdf


# %%
print("\nSection 005\n")

t = """
KIC 11442793. This system hosts seven transit- ing planets, ranging from slightly larger than Earth
to Jupiter-sized planets, all of which have been ei- ther validated or confirmed via transit-timing variations 
(Cabrera et al. 2014; Schmitt et al. 2014; Lissauer et al. 2014). 
We detect an additional TCE with a period of 14.4 days with S/N of 8.7."""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 11442793', author='Kepler')
print(search_result)

FoldTarget("KIC 11442793", 15, 1800, 14.4)

# %%


# %%
print("\nSection 006\n")

t = """
KIC 8480285. This system hosts two planets - one previously validated Earth-sized planet 
(Morton et al. 2016) and one candidate super-Earth in periods of 16.2 days and 29.7 days respectively. 
We detect an additional TCE at S/N = 8.4 with a period of 10.9 days."""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 8480285', author='Kepler')
print(search_result)

FoldTarget("KIC 8480285", 15, 1800, 16.2)

FoldTarget("KIC 8480285", 15, 1800, 29.7)

FoldTarget("KIC 8480285", 15, 1800, 10.9)

# %%



# %%
print("\nSection 007\n")

t = """
KIC 11568987. This system hosts one validated super-Earth (Morton et al. 2016) and one candidate Earth-sized planet in 16 and 7.3 day periods, respectively. 
We detect a new TCE with a period of 27.1 days, S/N = 9.8 and a depth corresponding to a roughly Earth-sized planet.
"""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 11568987', author='Kepler')
print(search_result)

FoldTarget("KIC 11568987", 15, 1800, 16)

FoldTarget("KIC 11568987", 15, 1800, 7.3)

FoldTarget("KIC 11568987", 15, 1800, 27.1)

# %%


# %%
print("\nSection 008\n")

t = """
KIC 4852528. This system hosts five confirmed small planets, the outer four of which are in a rare dynamical configuration with the middle three planets and the outer three planets each locked in a three-body resonance (Lissauer et al. 2014; Mac- Donald et al. 2016). We find a new TCE with a period of 14.64 days at S/N = 8.6 exterior to all five previously known planets.
"""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 4852528', author='Kepler')
print(search_result)

FoldTarget("KIC 4852528", 15, 1800, 14.64)

# %%


# %%
print("\nSection 009\n")

t = """
KIC 5972334. This system has two validated planets (Morton et al. 2016) and two candidates, including a warm Jupiter, 
an ultra-short-period planet (Sanchis-Ojeda et al. 2014), and a super- Earth in between, 
giving it a distinct architecture reminiscent of the 55 Cnc 
(Nelson et al. 2014) and WASP-47 (Becker et al. 2015) planetary systems. 
We detect a strong (S/N = 10.2) TCE with a period of 6.02 days which, if real, 
would have a radius a bit larger than that of the Earth. 
We see extra scatter in the light curve during the transit of this new TCE, however, 
which needs to be better un- derstood before this TCE can be considered a good candidate.
"""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 5972334', author='Kepler')
print(search_result)

FoldTarget("KIC 5972334", 15, 1800, 6.02)

# %%


# %%
print("\nSection 010\n")

t = """
KIC 10337517. This system hosts a sub-Earth- sized validated planet in a 4.3 day orbit (Morton et al. 2016) and a super-Earth-sized planet can- didate in a 7 day orbit. We detect an additional TCE in an 11.1 day orbit, which also would likely be sub-Earth-sized.
"""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 10337517', author='Kepler')
print(search_result)

FoldTarget("KIC 10337517", 15, 1800, 4.3)

FoldTarget("KIC 10337517", 15, 1800, 7)


# %%


# %%
print("\nSection 011\n")

t = """
KIC 11030475. This system has 4 planet can- didates, none of which have been validated. We detect a new TCE at an orbital period of 4.75 days and S/N = 8.7, with the same time of transit, depth, duration, and exactly half the period of a previously detected candidate, KOI 2248.02. This “new” TCE, therefore, appears to be the true pe- riod of the candidate KOI 2248.02.
"""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 11030475', author='Kepler')
print(search_result)

FoldTarget("KIC 11030475", 15, 1800, 4.75)


# %%


# %%
print("\nSection 012\n")

t = """
KIC 4548011. This system hosts one validated planet (Morton et al. 2016) in a 6.3 day period and one additional candidate in a 9.1 day period. We detect a new TCE with a period of 13.95 days with a relatively strong S/N = 9.7, which if real, would likely be a sub-Earth-sized planet. Subse- quently, this signal was also detected by the Kepler pipeline and included in the DR25 release of the Kepler planet candidate catalog as KOI 4288.04.
"""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 4548011', author='Kepler')
print(search_result)

FoldTarget("KIC 4548011", 15, 1800, 6.3)

FoldTarget("KIC 4548011", 15, 1800, 9.1)

FoldTarget("KIC 4548011", 15, 1800, 13.95)

# %%


# %%
print("\nSection 013\n")

t = """
KIC 3832474. This system hosts two confirmed giant planets, and a smaller confirmed super-Earth- sized planet called Kepler-30 b, which orbits with a period of 29.16 days (Fabrycky et al. 2012). All of these planets show transit timing variations (TTVs), and Kepler-30 d shows extremely large TTVs with an amplitude of more than a day (Fab- rycky et al. 2012; Panichi et al. 2017). We detect a new TCE with a period of 29.45 days - very close to the mean period of Kepler-30 d. We suspect this TCE is due to residual transits of Kepler-30 d that we were unable to completely remove from the light curve due to the large TTVs.
"""

print(t + "\n")

search_result = lk.search_lightcurve('KIC 3832474', author='Kepler')
print(search_result)

FoldTarget("KIC 3832474", 15, 1800, 29.16)

FoldTarget("KIC 3832474", 15, 1800, 29.45)


# %%

