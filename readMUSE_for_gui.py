import numpy as np
from astropy.io import fits


class MUSE_spec:

    def __init__(self, folderpath, obj_id):

        spec_file = folderpath + "combined_spec_" + obj_id + ".dat"

        self.lamb, self.flux, self.errflux = np.loadtxt(
            spec_file, usecols=(1, 2, 4), unpack=True)


class MUSE_img:

    def __init__(self, folderpath, image_id):

        img_NB = fits.open(folderpath + "NB_" + image_id + ".fits")
        self.nbdata = img_NB[0].data

        img_mask = fits.open(folderpath + "Obj_Mask_" + image_id + ".fits")
        self.maskdata = img_mask[0].data


class Cubex_table:

    def __init__(self, filename):

        self.npix, self.flux, self.eflux, \
            self.centroid, self.iso_area, self.lamb_size = np.loadtxt(
                filename, usecols=(1, 18, 19, 27, 14, 15), unpack=True)
        self.npix = self.npix.astype(int)
        self.iso_area = self.iso_area.astype(int)
        self.lamb_size = self.lamb_size.astype(int)

        self.ra, self.dec, self.id = np.loadtxt(filename, usecols=(25, 26, 0),
                                                unpack=True,
                                                dtype=str)
