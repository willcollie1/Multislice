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
    voltage = E_1.get()
    edose = E_2.get()
    abberation = E_3.get()
    unitcell = checkCmd
    filename = name
    View = checkCmd1
    topview = checkCmd4
    tile1 = checkCmd2
    til2 = checkCmd3

    file = open(filename, 'r')
    lines = file.readlines()
    file.close()

    parts = lines[3].split() 
    a = parts[5]
    parts1 = lines[4].split() 
    b = parts1[5]
    parts2 = lines[5].split()
    c = parts2[5]
    parts3 = lines[8].split()
    width = parts3[0]
    parts4 = lines[8].split()
    height = parts4[1]
    parts5 = lines[8].split()
    depth = parts5[2]


    fo = open("userinput.txt", "wb")
    fo.write("%f\n%f\n%f\n%d\n%d\n%d\n%f\n%d\n%f\n%d\n%s\n%d\n%d\n%d\n%d" % (float(a),float(b),float(c),int(width),int(height),int(depth),float(voltage), int(edose),float(abberation),checkCmd.get(),filename,checkCmd1.get(),checkCmd2.get(),checkCmd3.get(),checkCmd4.get()));
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
checkCmd2.set(0)

checkCmd3 = IntVar()
checkCmd3.set(0)

checkCmd4 = IntVar()
checkCmd4.set(0)

root.configure(background="#a1dbcd")
root.title("Multislice TEM Simulation")

label_1 = Label(root,text = "Electron Energy  (eV):")
label_1.grid(row=0)

label_2 = Label(root,text = "Electron Dose (" + angstrom + " / " + metresquare + "):")
label_2.grid(row=1)

label_3 = Label(root,text = "Spherical Abberation (m): ")
label_3.grid(row=2)

label_4 = Label(root,text = "Potential File: ")
label_4.grid(row=3,column=0)

E_1 = Entry(root)
E_1.grid(row=0,column=1)

E_2 = Entry(root)
E_2.grid(row=1,column=1)

E_3 = Entry(root)
E_3.grid(row=2,column=1)

c1 = Checkbutton(root,variable=checkCmd, onvalue=1, offvalue=0, text = "Hexagonal Unit Cell")
c1.grid(row=0, column = 4)

c2 = Checkbutton(root,variable=checkCmd1, onvalue=1, offvalue=0, text = "Side View               ")
c2.grid(row=1, column = 4)

c3 = Checkbutton(root,variable=checkCmd4, onvalue=1, offvalue=0, text = "Top View               ")
c3.grid(row=2, column = 4)

c4 = Checkbutton(root,variable=checkCmd2, onvalue=1, offvalue=0, text = "Verticle Retiling      ")
c4.grid(row=3, column = 4)

c5 = Checkbutton(root,variable=checkCmd3, onvalue=1, offvalue=0, text = "Quadrant Retiling   ")
c5.grid(row=4, column = 4)

B = Button(root, text = "OK", command = close_window)
B.grid(row=10,column=4)

B1 = Button(root,text='File Open', command=callback)
B1.grid(row=3,column=1)

root.mainloop()




