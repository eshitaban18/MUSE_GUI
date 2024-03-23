#! /home/eshita/anaconda3/bin/python

''' This is a GUI to check the redshift of MUSE cubes galaxies '''

import numpy as np
import tkinter as tk
from tkinter.ttk import Label, Entry
from tkinterdnd2 import DND_FILES, TkinterDnD
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from astropy.visualization import ImageNormalize, ZScaleInterval
from astropy.stats import sigma_clipped_stats
from astropy.convolution import convolve, Box2DKernel
from muse import MUSE_spec, MUSE_img, Cubex_table
import linelist as ll
import glob
import os


class My_gui:

    def __init__(self, fname):
        '''
        This part reads the data
        '''
        filename = fname.get()
        folder_path = filename.rsplit("/", 2)[0]

        img_filepath = folder_path + "/NB/"
        img_name = "NB_"\
            + filename.rsplit("/", 1)[1].split("_")[1].split(".")[0]\
            + ".fits"
        img_filename = img_filepath + img_name
        img = MUSE_img(img_filename)
        self.image = img.data

        # Smooth the image
        kernel = Box2DKernel(npix)
        self.image = convolve(self.image, kernel)

        mean, _, std = sigma_clipped_stats(self.image)
        self.vmin = mean - 2 * std
        self.vmax = mean + 2 * std

        this_id = filename.rsplit('/', 1)[1].split('.')[0].split("_")[1]
        cubex_datafile = glob.glob(folder_path + "/*.cat")
        cubex = Cubex_table(cubex_datafile[0])
        id_column = cubex.id.astype(str)
        match_id = this_id == id_column
        self.centroid = cubex.centroid[match_id][0]

        f = MUSE_spec(filename)
        self.lamb = f.lamb
        self.flux = f.flux

    def view_spectra(self):
        '''
        Plots the spectra and select the line we want to analyse
        and gives the centroid of that line
        '''
        fig1, (ax11, ax12) = plt.subplots(1, 2, facecolor=gui_color,
                                          figsize=(12., 3),
                                          width_ratios=[3, 1])
        plt.subplots_adjust(left=0.05, right=.95, wspace=.2)
        ind_plot = (self.lamb > self.centroid -
                    30) & (self.lamb < self.centroid+30)
        plot_flux = self.flux[ind_plot]
        plot_lamb = self.lamb[ind_plot]

        ax11.step(plot_lamb, plot_flux, where="mid", color="k")
        ax11.axvline(self.centroid, color="r", ls=":", lw=.6)
        ax11.set_ylim(min(plot_flux)-50, max(plot_flux)+50)

        norm = ImageNormalize(
            self.image, interval=ZScaleInterval(),
            vmin=self.vmin, vmax=self.vmax)
        ax12.imshow(self.image, interpolation="gaussian",
                    origin='lower', cmap=cmap, norm=norm)
        ax12.contour(self.image, levels=[self.vmax])
        ax11.set_facecolor(gui_color)
        ax12.set_facecolor(gui_color)

        global canvas_spec_plot
        if canvas_spec_plot:
            canvas_spec_plot.destroy()
        if canvas_veloplot:
            canvas_veloplot.destroy()
        if hidden_canvas:
            hidden_canvas.destroy()

        canvas1 = FigureCanvasTkAgg(fig1, root)
        canvas1.draw()
        canvas_spec_plot = canvas1.get_tk_widget()
        canvas_spec_plot.configure(bg=gui_color)
        canvas_spec_plot.place(x=20, y=70)
        plt.close(fig1)

    def veloplot(self, dv, line, l_list, x_place, y_place, row, column):
        '''
        Does the velocity plot corresponding to given centroid
        '''

        z = (self.centroid/line) - 1
        dz = (dv/3e5) * (1+z)

        fig, ax = plt.subplots(row, column, sharex=True,
                               figsize=(9, 5),
                               facecolor=gui_color)
        plt.subplots_adjust(left=0.05, right=.95, wspace=.3)

        for i in range(len(l_list)):

            lamb0 = float(ll.wl[l_list[i]])
            ind = (lamb0*(1+z-dz) < self.lamb) & (self.lamb < lamb0*(1+z+dz))
            velo_x = ((self.lamb[ind]/lamb0)-1-z) * (299792.46/(1+z))
            velo_y = self.flux[ind]

            if i <= (column-1):
                ind1 = 0
                ind2 = i
            else:
                ind1 = 1
                ind2 = i-column

            if (len(velo_y) != 0):
                ax[ind1, ind2].step(velo_x, velo_y, where="mid",
                                    color="k")
                ax[ind1, ind2].set_xlim(-500, 500)
                ax[ind1, ind2].axvline(
                    x=0, ymax=1, ymin=0, ls=":", color="r", lw=.6)
                ax[ind1, ind2].set_title(l_list[i],
                                         fontsize=12, color="r")

                if (i == len(l_list)-1):
                    tik = (1 + z) * ll.wl["OII_3729"]
                    second_line = ((tik/ll.wl["OII_3727"]) -
                                   1-z) * (299792.46/(1+z))
                    ax[ind1, ind2].axvline(x=second_line, ymax=1, ymin=0,
                                           ls=":", color="r", lw=0.6)

            if ll.wl[l_list[i]] == line:
                ax[ind1, ind2].set_facecolor("lightcyan")
            else:
                ax[ind1, ind2].set_facecolor(gui_color)

        global canvas_veloplot
        if canvas_veloplot:
            canvas_veloplot.destroy()
        if hidden_canvas:
            hidden_canvas.destroy()

        canvas = FigureCanvasTkAgg(fig, root)
        canvas.draw()
        canvas_veloplot = canvas.get_tk_widget()
        canvas_veloplot.configure(bg=gui_color)
        canvas_veloplot.place(x=x_place, y=y_place)
        plt.close(fig)

        return z

    def hidden_veloplot(self, dv, z, l_list, x_place, y_place, row, column):
        '''
        Does the velocity plot corresponding to given centroid
        '''
        dz = (dv/3e5) * (1+z)

        fig_hid, axh = plt.subplots(row, column, sharex=True, figsize=(9, 5),
                                    facecolor=gui_color)
        plt.subplots_adjust(left=0.05, right=0.95, wspace=.3)

        for i in range(len(l_list)):

            lamb0 = float(ll.wl[l_list[i]])
            ind = (lamb0*(1+z-dz) < self.lamb) & (self.lamb < lamb0*(1+z+dz))
            velo_x = ((self.lamb[ind]/lamb0)-1-z) * (299792.46/(1+z))
            velo_y = self.flux[ind]

            if i <= (column-1):
                ind1 = 0
                ind2 = i
            else:
                ind1 = 1
                ind2 = i-column

            if (len(velo_y) != 0):
                axh[ind1, ind2].step(velo_x, velo_y, where="mid",
                                     color="k")
                axh[ind1, ind2].set_xlim(-500, 500)
                axh[ind1, ind2].axvline(
                    x=0, ymax=1, ymin=0, ls=":", color="r", lw=.6)
                axh[ind1, ind2].set_title(l_list[i],
                                          fontsize=12, color="r")

            axh[ind1, ind2].set_facecolor(gui_color)

        global hidden_canvas
        canvas = FigureCanvasTkAgg(fig_hid, root)
        canvas.draw()
        hidden_canvas = canvas.get_tk_widget()
        hidden_canvas.configure(bg=gui_color)
        hidden_canvas.place(x=x_place, y=y_place)
        plt.close(fig_hid)


def press(fname, num):
    ''' This acts the press/click operation of button,
    executes the velocity plot function defined under My_gui class '''

    line = ll.wl[str(num)]
    line_list = ["Ha", "Hb", "Hg", "OIII_5008", "OIII_4960", "OII_3727"]
    global test_z
    test_z = np.round(My_gui(fname).veloplot(
        dv, line, line_list, 20, 370, 2, 3), 6)
    Label(text="Redshift:", font=myfont).place(x=1200, y=20)
    if test_z > 0:
        tk.Label(text=test_z, font=myfont, fg="#000000",
                 width=10,
                 borderwidth=1.2, relief="solid").place(x=1330, y=17)
    else:
        tk.Label(text=test_z, font=myfont, fg="red",
                 width=10,
                 borderwidth=1.2, relief="solid").place(x=1330, y=17)


def show_hidden_plots(fname, show=False):
    ''' This acts the press/click operation of button,
    executes the hidden velocity plot part '''
    line_list_hidden = ["Hd", "NII_6585", "NII_6549",
                        "SII_6718", "SII_6732", "OI_6365"]
    My_gui(fname).hidden_veloplot(
        dv, test_z, line_list_hidden, 910, 370, 2, 3)


def gui_body(fname):
    '''
    contains main decorating elements of the GUI
    '''

    plot_button = tk.Button(root,
                            command=lambda: My_gui(fname).view_spectra(),
                            text="Plot",  height=1, width=6,
                            font=myfont, borderwidth=3)
    plot_button.place(x=860, y=15)

    comment_label = Label(root, text="Comment:", font=myfont)
    comment_label.place(x=20, y=950)
    comment = Entry(root, width=60, font=myfont)
    comment.pack(padx=5, side=tk.BOTTOM)
    comment.place(x=140, y=950)

    save_button = tk.Button(root,
                            command=lambda: save_text(fname, comment),
                            text="Save",  height=1, width=7,
                            font=myfont, borderwidth=3)
    save_button.place(x=890, y=940)


def gui_grid_button(fname):
    '''
    contains all the buttons on side bar
    '''

    button_1 = tk.Button(root, command=lambda: press(
        fname, "Ha"), text="H-a",
        height=1, width=8, font=myfont, borderwidth=3)
    button_1.place(x=1700, y=20)

    button_2 = tk.Button(root, command=lambda: press(
        fname, "Hb"), text="H-b",
        height=1, width=8, font=myfont, borderwidth=3)
    button_2.place(x=1700, y=60)

    button_3 = tk.Button(root, command=lambda: press(
        fname, "Hg"), text="H-g",
        height=1, width=8, font=myfont, borderwidth=3)
    button_3.place(x=1700, y=100)

    button_4 = tk.Button(root, command=lambda: press(
        fname, "OIII_5008"), text="OIII_5008",
        height=1, width=8, font=myfont, borderwidth=3)
    button_4.place(x=1700, y=140)

    button_5 = tk.Button(root, command=lambda: press(
        fname, "OIII_4960"), text="OIII_4960",
        height=1, width=8, font=myfont, borderwidth=3)
    button_5.place(x=1700, y=180)

    button_6 = tk.Button(root, command=lambda: press(
        fname, "OII_3727"), text="OII_3727",
        height=1, width=8, font=myfont, borderwidth=3)
    button_6.place(x=1700, y=220)

    button_show = tk.Button(root, command=lambda: show_hidden_plots(
        fname, show=True), text="Show All",
        height=1, width=8, font=myfont, borderwidth=3)
    button_show.place(x=1700, y=300)


def save_text(fname, comment):
    ''' Saves the comment about any object in My_comments.txt file '''

    obj = fname.get().rsplit('/', 1)[1].split('.')[0].split("_")[1]
    path = fname.get().rsplit('/', 2)[0] + '/'
    completename = path + 'My_comments.txt'

    try:
        with open(completename, "x") as f:
            L = "#Object_ID, Redshift, Comment \n"
            f.writelines(L)
            f.write(obj + "   " + str(np.round(test_z, 6)) +
                    "    " + comment.get() + "\n")

    except FileExistsError:
        with open(completename, "a") as f:
            f.write(obj + "   " + str(np.round(test_z, 6)) +
                    "    " + comment.get() + "\n")


def main():

    root.title("MUSEQuBE_GUI")
    root.geometry("1550x1200")
    # root.maxsize(900, 1000)

    name_label = Label(root, text="Enter filename:", font=myfont)
    name_label.place(x=20, y=20)
    fname = Entry(root, width=50, font=myfont)

    def drop(e):
        fname.insert(0, e.data)
        num_ch = len(e.data)
        fname.delete(num_ch, tk.END)

    fname.drop_target_register(DND_FILES)
    fname.dnd_bind('<<Drop>>', drop)
    fname.place(x=220, y=20)

    gui_body(fname)
    gui_grid_button(fname)

    def change():
        foo = fname.get()
        foo_path = foo.rsplit("/", 1)[0]
        a = sorted(glob.glob(foo_path + "/spec_*.dat"),
                   key=lambda x: int(
                       os.path.basename(x).split("_")[1].split(".")[0]))
        for i in range(len(a)):
            if foo == a[i]:
                try:
                    fname.delete(0, tk.END)
                    fname.insert(0, a[i+1])
                except IndexError:
                    tk.messagebox.showerror(
                        'Error',
                        'Error: That was all !')

    tk.Button(
        root, text="Next", command=lambda: change(),
        font=myfont, borderwidth=3).place(x=1000, y=15)

    tk.Button(
        root, text="Exit", command=root.destroy,
        font=myfont, borderwidth=3).place(x=1400, y=950)

    root.mainloop()


if __name__ == "__main__":
    root = TkinterDnD.Tk()
    myfont = (("Times New Roman"), 18)
    test_z = ""
    canvas_veloplot = None
    canvas_spec_plot = None
    hidden_canvas = None
    dv = 500
    gui_color = "#D9D9D9"
    cmap = "inferno"
    npix = 2
    main()
