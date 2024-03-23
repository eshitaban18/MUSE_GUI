To use this GUI on a Cubex extracted datacube:

1. There must be a folder containing: 


    (a) A folder named "NB" that will contain narrow band and object masked images of all the objects within that cube extracted by CubEx.
        The narrow band images will have a name type: "NB_{obj_id}.fits"
        The object masked images will have a name type: "Obj_Mask_{obj_id}.fits"
        (Names are case sensitive)

    (b) A folder named "spec" that will contain the optimally extracted spectra for all the objects within that cube extracted by Cubex.
        The files will have a name type: "combined_spec_{obj_id}.dat"

    (c) Finally, it must have the final object catalog extracted by CubEx that contains all the source Id, centroid, ra ,dec etc of the 
        extracted sources. The extension of the file should be "*.cat". Make sure there is no other file in that specific folder with  the same extension.

To run in terminal: python gui.py {the complete name of the .cat file described above in (c)} 


Your comments will be saved in a file named "My_comments.txt". The columns are (run_id, object_id, comment). Multiple comments on the same object saved at different instances will be saved as different entries in separate lines.

-----------------------------
Keyboard shortcut:
-----------------------------
<Enter> => Plot
<Up> => Next
<Down> => Back
<Control-a> => Save your comment
<Esc> => Exit

