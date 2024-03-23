#! /home/eshita/anaconda3/bin/python

''' This is a GUI to check the redshift of MUSE cubes galaxies '''

import glob
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import gui_linelist as ll
import tkinter as tk
from tkinter.ttk import Label, Entry
from matplotlib.backends.backend_tkagg \
    import FigureCanvasTkAgg, NavigationToolbar2Tk
from astropy.visualization import ImageNormalize, ZScaleInterval
from astropy.stats import sigma_clipped_stats
from astropy.convolution import convolve, Box2DKernel
from readMUSE_for_gui import MUSE_spec, MUSE_img, Cubex_table
import warnings
plt.rcdefaults()
plt.rcParams["axes.linewidth"] = 1.2


class My_gui:

    def __init__(self, spec_list, runid):
        '''
        This part reads the data
        '''
        self.spec_list = spec_list
        self.runid = runid
        self.filename = spec_list[int(runid.get())]
        self.folder_path = self.filename.rsplit("/", 2)[0]

        # ----------------------------------------------------------
        self.this_id = os.path.basename(self.filename).split(".dat")[
            0].rsplit("_", 1)[1]
        cubex_datafile = glob.glob(self.folder_path + "/*.cat")
        self.datacube = Cubex_table(cubex_datafile[0])
        id_column = self.datacube.id.astype(str)
        self.match_id = self.this_id == id_column
        self.centroid = self.datacube.centroid[self.match_id][0]

        # ----------------------------------------------------------
        spec_filepath = self.folder_path + "/spec/"
        f = MUSE_spec(spec_filepath, self.this_id)
        self.lamb = f.lamb
        self.flux = f.flux
        self.errflux = f.errflux

        # ------------------------------------------------------------

    def view_spectra(self):
        '''
        Plots the spectra and select the line we want to analyse
        and gives the centroid of that line
        '''
        fig1, (ax12, ax13) = plt.subplots(1, 2, facecolor=gui_color,
                                          figsize=(7., 3),
                                          width_ratios=[1, 1])
        plt.subplots_adjust(left=0.05, right=.95, wspace=.2)

        # Display the image
        # ----------------------------------------------------------

        img_filepath = self.folder_path + "/NB/"
        img_name = self.this_id

        img = MUSE_img(img_filepath, img_name)
        image = img.nbdata
        mask = img.maskdata

        # Smooth the image

        kernel = Box2DKernel(npix)
        image = convolve(image, kernel)
        mean, _, std = sigma_clipped_stats(image)
        vmin = mean - 2 * std
        vmax = mean + 2 * std
        norm = ImageNormalize(
            image, interval=ZScaleInterval(),
            vmin=vmin, vmax=vmax)
        ax12.imshow(image, interpolation="gaussian",
                    origin='lower', cmap=cmap, norm=norm)
        ax12.contour(mask, levels=[0], colors="red", linewidths=2.5)
        ax12.tick_params(which="both", left=False, right=False,
                         labelleft=False,
                         labelbottom=False, bottom=False, top=False)
        ax12.set_facecolor(gui_color)

        # Display some parameters
        # ------------------------------------------------------------

        nvox = self.datacube.npix[self.match_id][0]
        iso_area = self.datacube.iso_area[self.match_id][0]
        lamb_size = self.datacube.lamb_size[self.match_id][0]
        ra = self.datacube.ra[self.match_id][0]
        dec = self.datacube.dec[self.match_id][0]

        z_LAE = (self.centroid/1215.67) - 1
        text_view = "Obj ID = %i" % (int(self.this_id)) + '\n' + \
            "if LAE, z = %.5f" % (z_LAE) + '\n' + \
            "Npix = %i" % (nvox) + '\n' + \
            "RA = " + str(ra) + '\n' + \
            "Dec = " + str(dec) + '\n' + \
            "Iso area = %i " % (iso_area) + '\n' + \
            r"$\lambda_{\rm size}$ = %i " % (lamb_size)
        ax13.text(.08, .2, text_view, fontsize=16)
        ax13.tick_params(which="both", left=False, right=False,
                         labelleft=False,
                         labelbottom=False, bottom=False, top=False)
        ax13.set_facecolor(gui_color)

        # placing the image
        # ------------------------------------------------------------

        global canvas_img_plot
        if canvas_img_plot:
            canvas_img_plot.destroy()
        if canvas_veloplot:
            canvas_veloplot.destroy()
        if canvas_tick_plot:
            canvas_tick_plot.destroy()

        canvas_img = FigureCanvasTkAgg(fig1, root)
        canvas_img.draw()
        canvas_img_plot = canvas_img.get_tk_widget()
        canvas_img_plot.configure(bg=gui_color)
        canvas_img_plot.place(x=900, y=100)

        plt.close(fig1)
        press(self.spec_list, self.runid, "OII_3727")

        # --------------------------------------------------------------

    def veloplot(self, dv, line, l_list, x_place, y_place, row, column):
        '''
        Does the velocity plot corresponding to given centroid
        '''
        z = (self.centroid/line) - 1
        dz = (dv/3e5) * (1+z)

        # Tick plot
        # --------------------------------------------------------------

        fig_tick, ax_tick = plt.subplots(facecolor=gui_color,
                                         figsize=(9., 3))
        plt.subplots_adjust(left=0.05, right=.95)

        ax_tick.step(self.lamb, self.flux, where="mid", color="k", lw=.8)
        ax_tick.step(self.lamb, self.errflux,
                     where="mid", color="c", lw=.8)
        ax_tick.axhline(y=0, color="red", lw=1, ls="dashed")
        ax_tick.axvline(self.centroid, color="m", ls="-.",
                        lw=1.6, zorder=300)

        ax_tick.set_ylim(min(self.flux)-50, max(self.flux)+50)
        ax_tick.set_xlim(min(self.lamb)-100, max(self.lamb)+100)
        ax_tick.set_facecolor("white")

        # Velocity plot
        # --------------------------------------------------------------

        fig, ax = plt.subplots(row, column, sharex=True,
                               figsize=(19, 5),
                               facecolor=gui_color)
        plt.subplots_adjust(left=0.05, right=.95, wspace=.3)

        for i in range(len(l_list)):
            lamb0 = float(ll.wl[l_list[i]])
            ind = (lamb0*(1+z-dz) < self.lamb) & (self.lamb < lamb0*(1+z+dz))
            velo_x = ((self.lamb[ind]/lamb0)-1-z) * (3e5/(1+z))
            velo_y = self.flux[ind]

            if i <= (column-1):
                ind1 = 0
                ind2 = i
            else:
                ind1 = 1
                ind2 = i-column

            if (len(velo_y) != 0):
                ax[ind1, ind2].step(velo_x, velo_y, where="mid", color="k")
                ax[ind1, ind2].set_xlim(-dv, dv)
                ax[ind1, ind2].axvline(
                    x=0, ymax=1, ymin=0, ls=":", color="r", lw=.9)
                ax[ind1, ind2].set_title(l_list[i],
                                         fontsize=12, color="r")
                if (i == 8):
                    tik = (1 + z) * ll.wl["OII_3729"]
                    second_line = ((tik/ll.wl["OII_3727"]) -
                                   1-z) * (3e5/(1+z))
                    ax[ind1, ind2].axvline(x=second_line, ymax=1, ymin=0,
                                           ls=":", color="r", lw=0.9)

            if ll.wl[l_list[i]] == line:
                ax[ind1, ind2].set_facecolor("lightcyan")
            else:
                ax[ind1, ind2].set_facecolor(gui_color)

            plot_tick(l_list[i], lamb0, z, ax_tick)

        # Placing the images
        # ------------------------------------------------------------

        global canvas_veloplot
        if canvas_veloplot:
            canvas_veloplot.destroy()

        canvas_velo = FigureCanvasTkAgg(fig, root)
        canvas_velo.draw()
        canvas_veloplot = canvas_velo.get_tk_widget()
        canvas_veloplot.configure(bg=gui_color)
        canvas_veloplot.place(x=x_place, y=y_place)

        global canvas_tick_plot
        if canvas_tick_plot:
            canvas_tick_plot.destroy()

        canvas_tick = FigureCanvasTkAgg(fig_tick, root)
        canvas_tick.draw()
        canvas_tick_plot = canvas_tick.get_tk_widget()
        canvas_tick_plot.configure(bg=gui_color)
        canvas_tick_plot.place(x=20, y=100)

        # Placing the zooming and other tools
        # ------------------------------------------------------------------

        toolbar = NavigationToolbar2Tk(canvas_tick, root, pack_toolbar=True)
        toolbar.update()
        toolbar.place(x=20, y=100)
        plt.close(fig)
        plt.close(fig_tick)

        return z


def plot_tick(line_name, lamb0, z, ax):
    '''
    Plots ticks based on the current redshift as decided by
    pressing the buttons
    '''
    x = (1+z)*lamb0
    ax.axvline(x, ls="dashed", lw=1.6, color="blue")
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.text(x+4, .6, line_name, transform=trans,
            fontsize=9, fontweight='bold', rotation=90,
            color="red")


def veloplot_button(runid, spec_list, line_id, line_text, xpos, ypos):
    '''
    Creates the buttons on sidebar
    '''
    button = tk.Button(root, command=lambda: press(
        spec_list, runid, line_id), text=line_text,
        height=1, width=8, font=myfont, borderwidth=3).place(x=xpos, y=ypos)
    return button


def press(spec_list, runid, num):
    ''' This acts the press/click operation of button,
    executes the velocity plot function defined under My_gui class '''

    line = ll.wl[str(num)]
    column = int(len(line_list)/2)

    global test_z
    test_z = np.round(My_gui(spec_list, runid).veloplot(
        dv, line, line_list, 20, 400, 2, column), 6)

    if test_z > 0:
        tk.Label(text=test_z, font=myfont, fg="#000000",
                 width=10,
                 borderwidth=1.2, relief="solid").place(x=1330, y=17)
    else:
        tk.Label(text=test_z, font=myfont, fg="red",
                 width=10,
                 borderwidth=1.2, relief="solid").place(x=1330, y=17)


def save_text(spec_list, runid, comment):
    ''' Saves the comment about any object in My_comments.txt file '''

    now_name = spec_list[int(runid.get())]
    obj = os.path.basename(now_name).split('.dat')[0].rsplit("_", 1)[1]
    path = now_name.rsplit('/', 2)[0] + '/'
    completename = path + 'My_comments.txt'

    try:
        with open(completename, "x") as f:
            L = "#Object_ID, Your_comment \n"
            f.writelines(L)
            f.write(obj + "   " + comment.get() + "\n")

    except FileExistsError:
        with open(completename, "a") as f:
            f.write(obj + "   " + comment.get() + "\n")


def gui_body(spec_list, runid):
    '''
    contains main decorating elements of the GUI
    '''

    plot_button = tk.Button(root,
                            command=lambda: My_gui(
                                spec_list, runid).view_spectra(),
                            text="Plot",  height=1, width=6,
                            font=myfont, borderwidth=3)
    plot_button.place(x=900, y=15)
    root.bind('<Return>',
              lambda event: My_gui(spec_list, runid).view_spectra())

    comment_label = Label(root, text="Comment:", font=myfont)
    comment_label.place(x=20, y=950)
    comment = Entry(root, width=60, font=myfont)
    comment.pack(padx=5, side=tk.BOTTOM)
    comment.place(x=140, y=950)

    save_button = tk.Button(root,
                            command=lambda: save_text(
                                spec_list, runid, comment),
                            text="Save",  height=1, width=7,
                            font=myfont, borderwidth=3)
    save_button.place(x=890, y=940)


def gui_grid_button(spec_list, runid):
    '''
    contains all the buttons on side bar
    '''
    veloplot_button(runid, spec_list, "Ha", "H-a", 1700, 20)
    veloplot_button(runid, spec_list, "Hb", "H-b", 1700, 60)
    veloplot_button(runid, spec_list, "Hg", "H-g", 1700, 100)
    veloplot_button(
        runid, spec_list, "OIII_5008", "OIII_5008", 1700, 140)
    veloplot_button(
        runid, spec_list, "OIII_4960", "OIII_4960", 1700, 180)
    veloplot_button(
        runid, spec_list, "OII_3727", "OII_3727", 1700, 220)
    veloplot_button(
        runid, spec_list, "SII_6732", "SII_6732", 1700, 260)
    veloplot_button(
        runid, spec_list, "NII_6585", "NII_6585", 1700, 300)
    veloplot_button(
        runid, spec_list, "NeIII_3869", "NeIII_3869", 1700, 340)


def main():

    root.title("MUSE_E-Z")
    root.geometry("1550x1200")

    # Decorate
    # -------------------------------------------------------------
    name_label = Label(root, text="Object ID file:", font=myfont)
    name_label.place(x=20, y=20)
    dfile_name = Entry(root, width=50, font=myfont)
    dfile_name.place(x=220, y=20)

    runid_label = Label(root, text="You are at:", font=myfont)
    runid_label.place(x=290, y=60)
    runid = Entry(root, width=4, font=myfont)
    runid.place(x=440, y=60)

    tot_obj_label = Label(root, text="Total object:", font=myfont)
    tot_obj_label.place(x=20, y=60)
    tot_obj = Entry(root, width=7, font=myfont)
    tot_obj.place(x=190, y=60)

    # Entering the filename
    # -------------------------------------------------------------------
    dfile_name.insert(0, str(sys.argv[1]))
    foo = dfile_name.get()
    foo_path = foo.rsplit("/", 1)[0]
    spec_list = sorted(glob.glob(foo_path + "/spec/combined_spec_*.dat"),
                       key=lambda x: int(
        os.path.basename(x).rsplit("_", 1)[1].split(".dat")[0]))
    runid.insert(0, "0")
    tot_obj.insert(0, "0 to " + str(len(spec_list)-1))
    tot_obj.config(state="disabled", foreground="black", background="white")

    Label(text="Redshift:", font=myfont).place(x=1200, y=20)

    # Execute the GUI
    # --------------------------------------------------------------------

    gui_body(spec_list, runid)
    gui_grid_button(spec_list, runid)

    # Creating "next" and "back" event
    # --------------------------------------------------------------------

    def next_button():
        run_num = int(runid.get())
        for i in range(len(spec_list)):
            if run_num == i:
                try:
                    runid.delete(0, tk.END)
                    runid.insert(0, i+1)
                    root.bind('<Return>', My_gui(
                        spec_list, runid).view_spectra())

                except IndexError:
                    runid.delete(0, tk.END)
                    runid.insert(0, i)
                    root.bind('<Return>', My_gui(
                        spec_list, runid).view_spectra())
                    tk.messagebox.showerror(
                        'Error',
                        'Error: That was all !')

    def back_button():
        run_num = int(runid.get())
        for i in range(len(spec_list)):
            if run_num == i and run_num > 0:
                runid.delete(0, tk.END)
                runid.insert(0, i-1)
                root.bind('<Return>', My_gui(
                    spec_list, runid).view_spectra())

    tk.Button(
        root, text=u'\u2193', command=lambda: back_button(),
        font=myfont, borderwidth=3).place(x=860, y=15)
    root.bind("<Down>",
              lambda event: back_button())

    tk.Button(
        root, text=u'\u2191', command=lambda: next_button(),
        font=myfont, borderwidth=3).place(x=990, y=15)
    root.bind("<Up>",
              lambda event: next_button())

    # Kill the interface
    # ----------------------------------------------------
    def exit():
        res = tk.messagebox.askquestion('Exit Tab',
                                        'Exit Interface?')
        if res == 'yes':
            root.destroy()
        else:
            pass

    tk.Button(
        root, text="Exit", command=lambda: exit(),
        font=myfont, borderwidth=3).place(x=1400, y=950)
    root.bind("<Escape>", lambda event: exit())

    root.mainloop()


if __name__ == "__main__":

    warnings.filterwarnings("ignore")
    root = tk.Tk()

    gui_color = "#D9D9D9"
    cmap = "viridis"
    myfont = (("Times New Roman"), 18)

    test_z = ""
    canvas_veloplot = None
    canvas_img_plot = None
    canvas_tick_plot = None

    dv = 800
    npix = 3
    line_list = ["Ha", "Hb", "Hg", "NII_6585", "NII_6549", "NeIII_3869",
                 "OIII_5008", "OIII_4960", "OII_3727", "SII_6718",
                 "SII_6732", "OI_6365"]

    main()
