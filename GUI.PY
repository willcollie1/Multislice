from Tkinter import *
from tkFileDialog   import askopenfilename
from ttk import * 

def callback():
    global name
    name = askopenfilename()


def close_window():
    global entry
    voltage = E_1.get()
    edose = E_2.get()
    abberation = E_3.get()
    unitcell = checkCmd
    filename = name

    tile = variable1.get()
    if(tile == "No Retiling"):
       tile = 0
    elif(tile == "Verticle Retiling"):
       tile = 1
    elif(tile == "Quadrant Retiling"):
       tile = 2

    view = variable.get()
    if(view == "Front View"):
       view = 0
    elif(view == "Side View"):
       view = 1
    elif(view == "Top View"):
       view = 2

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
    fo.write("%f\n%f\n%f\n%d\n%d\n%d\n%f\n%d\n%f\n%d\n%s\n%d\n%d" % (float(a),float(b),float(c),int(width),int(height),int(depth),float(voltage), int(edose),float(abberation),checkCmd.get(),filename,int(tile),int(view)));
    fo.close()
    root.destroy()

root = Tk()

root.style = Style()
root.style.theme_use("clam")

angstrom = u'\u212B'.encode('utf-8')
metresquare = u"\u33A1".encode('utf-8')

global variable
variable = StringVar(root)
variable.set("Front Value") 
option = OptionMenu(root, variable, "Front View", "Side View", "Top View")
option.grid(row=1, column = 4)

global variable1
variable1 = StringVar(root)
variable1.set("No Retiling")
option1 = OptionMenu(root, variable1, "No Retiling", " Vericle Retiling", "Quadrant Retiling")
option1.grid(row=2, column = 4)

checkCmd = IntVar()
checkCmd.set(0)

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

B = Button(root, text = "OK", command = close_window)
B.grid(row=10,column=4)

B1 = Button(root,text='File Open', command=callback)
B1.grid(row=3,column=1)

root.mainloop()





