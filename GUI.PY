from Tkinter import *
from tkFileDialog   import askopenfilename
from ttk import * 

#Allow the user to search for their potential file
def callback():
    global name
    name = askopenfilename()

#Outputting all userinputs to a file that can be used by the C program
def close_window():
    global entry
    a = E_1.get()
    b = E_2.get()
    c = E_3.get()
    width = E_4.get()
    height = E_5.get()
    depth = E_6.get()
    voltage = E_7.get()
    edose = E_8.get()
    abberation = E_9.get()
    unitcell = checkCmd
    filename = name
    View = checkCmd1
    tile1 = checkCmd2
    til2 = checkCmd3

    fo = open("userinput.txt", "wb")
    fo.write("%f\n%f\n%f\n%d\n%d\n%d\n%f\n%d\n%f\n%d\n%s\n%d\n%d\n%d" % (float(a),float(b),float(c),int(width),int(height),int(depth),float(voltage), int(edose),float(abberation),checkCmd.get(),filename,checkCmd1.get(),checkCmd2.get(),checkCmd3.get()));
    fo.close()
    root.destroy()

root = Tk()

#Choosing the style of the GUI
root.style = Style()
root.style.theme_use("clam")

angstrom = u'\u212B'.encode('utf-8')
metresquare = u"\u33A1".encode('utf-8')

checkCmd = IntVar()
checkCmd.set(0)

checkCmd1 = IntVar()
checkCmd1.set(0)

checkCmd2 = IntVar()
checkCmd1.set(0)

checkCmd3 = IntVar()
checkCmd1.set(0)

root.configure(background="#a1dbcd")
root.title("Multislice TEM Simulation")

label_1 = Label(root,text = "Lattice parameter a  (" + angstrom + "):" )
label_1.grid(row=0)

label_2 = Label(root,text = "Lattice parameter b  (" + angstrom + "):" )
label_2.grid(row=1)

label_3 = Label(root,text = "Lattice parameter c  (" + angstrom + "):" )
label_3.grid(row=2)

label_4 = Label(root,text = "Width (Voxels):")
label_4.grid(row=3)

label_5 = Label(root,text = "Height (Voxels): ")
label_5.grid(row=4)

label_6 = Label(root,text = "Depth (Voxels): ")
label_6.grid(row=5)

label_7 = Label(root,text = "Electron Energy  (eV):")
label_7.grid(row=6)

label_8 = Label(root,text = "Electron Dose (" + angstrom + " / " + metresquare + "):")
label_8.grid(row=7)

label_9 = Label(root,text = "Spherical Abberation (m): ")
label_9.grid(row=8)

label_10 = Label(root,text = "Potential File: ")
label_10.grid(row=9)

E_1 = Entry(root)
E_1.grid(row=0,column=1)

E_2 = Entry(root)
E_2.grid(row=1,column=1)

E_3 = Entry(root)
E_3.grid(row=2,column=1)

E_4 = Entry(root)
E_4.grid(row=3,column=1)

E_5 = Entry(root)
E_5.grid(row=4,column=1)

E_6 = Entry(root)
E_6.grid(row=5,column=1)

E_7 = Entry(root)
E_7.grid(row=6,column=1)

E_8 = Entry(root)
E_8.grid(row=7,column=1)

E_9 = Entry(root)
E_9.grid(row=8,column=1)

c1 = Checkbutton(root,variable=checkCmd, onvalue=1, offvalue=0, text = "Hexagonal Unit Cell")
c1.grid(row=0, column = 2)

c2 = Checkbutton(root,variable=checkCmd1, onvalue=1, offvalue=0, text = "Side View")
c2.grid(row=1, column = 2)

c3 = Checkbutton(root,variable=checkCmd2, onvalue=1, offvalue=0, text = "Verticle Retiling")
c3.grid(row=2, column = 2)

c4 = Checkbutton(root,variable=checkCmd3, onvalue=1, offvalue=0, text = "Quadrant Retiling")
c4.grid(row=3, column = 2)

B = Button(root, text = "OK", command = close_window)
B.grid(row=7,column=2)

B1 = Button(root,text='File Open', command=callback)
B1.grid(row=9,column=1)

root.mainloop()




