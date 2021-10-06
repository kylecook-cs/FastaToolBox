#IMPORTS
    #Bio libs
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Seq import Seq
    #Tkinter for GUI libs
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import *
    #Pillow for banner
from PIL import Image, ImageTk
    #Operating System for CLustalw and files
import os
    #Numpy for Percent Identity Function
import numpy as np

#GUI 
## Create the window object
window = tk.Tk()

#Load the logo image for banner
banner = ImageTk.PhotoImage(Image.open("pics\sjsu_banner.png"))

#Create font
font = "Luminari"

fileName = None

#FUNCTIONS
## This function takes in a fasta file and align the sequences
def alignFile():
        # try block to prompt user to open aln file the was created by Clustalw
        try:
            alignFile = filedialog.askopenfilename(initialdir="C:\\Users\\cook1\\OneDrive\\Desktop\\CS123B",
                                          title="Select a file",
                                          filetypes= (("aln files","*.aln"),
                                          ("all files","*.*")))
            # Read the aln file and to preview alignments
            alignments = AlignIO.read(alignFile,"clustal")
            messagebox.showinfo('Aligned Fasta File', alignments)
        except:
            # Raise Error and Show user a error message
            messagebox.showerror('Error', "Please pick an aln file!")
## This function takes in a dnd file and creates a tree   
def drawTree():
    # Ask user to insert the dnd file that clustalw created
    treeFile = filedialog.askopenfilename(initialdir="C:\\Users\\cook1\\OneDrive\\Desktop\\CS123B",
                                          title="Select a file",
                                          filetypes= (("dnd files","*.dnd"),
                                          ("all files","*.*")))
    # try block to check it the correct file was loaded, if so draw tree
    try:
        tree = Phylo.read(treeFile, "newick")
        Phylo.draw(tree)
    except:
        # Raise error and Show user a error message
        messagebox.showerror('Error', "Please pick the correct file type (.dnd)")

## This function creates a distance matrix 
def createMatrix():
    fileName
    # try block to check if file was loaded if so create and display preview of matrix
    try:
        alignment = AlignIO.read(fileName,"fasta")
        calculator = DistanceCalculator('identity')
        matrix = calculator.get_distance(alignment)
        messagebox.showinfo('Distance Matrix', matrix)
    except:
        # Raise error and Show user a error message
        messagebox.showerror('Error', "A Fasta file must be loaded first")

## This function takes in 2 sequences as strings and computes the Identity Percentage between them
def percentIdentity():
    #Strings from entry boxes
    string1 = sequence1.get()
    string2 = sequence2.get()

    # Initialize matrix of zeros
    rows = len(string1)+1
    columns = len(string2)+1
    distance_matrix = np.zeros((rows,columns),dtype = int)

    # Populate matrix of zeros with the indeces of each character of both strings
    for i in range(1, rows):
        for j in range(1,columns):
            distance_matrix[i][0] = i
            distance_matrix[0][j] = j

    # Iterate over the matrix to compute the cost of deletions,insertions and/or substitutions    
    for c in range(1, columns):
        for r in range(1, rows):
            if string1[r-1] == string2[c-1]:
                # same position in strings, move is 0
                moves = 0 
            else:
                # either insertion, deletion, substitution - move counts as 2 due to multiple steps for 1 procedure
                moves = 2 
            distance_matrix[r][c] = min(distance_matrix[r-1][c] + 1,    # deletions moves
                                 distance_matrix[r][c-1] + 1,           # insertions moves
                                 distance_matrix[r-1][c-1] + moves)     # substitutions moves
    
    # length of 2 sequences - the moves with matrix / length of 2 sequences * 100 to get the percentage
    percent = ((len(string1)+len(string2)) - distance_matrix[r][c]) / (len(string1)+len(string2)) * 100
    #Format and show result in messagebox
    messagebox.showinfo('Levenshtein Ratio',"The identity similiarity of these sequences is {:.2f} %".format(percent))
    
## This function prompts user to select file for directory and saves as global for other functions
def openFile():
    global fileName
    fileName = filedialog.askopenfilename(initialdir="C:\\Users\\cook1\\OneDrive\\Desktop\\CS123B",
                                          title="Select a file",
                                          filetypes= (("fasta files 1","*.FASTA"),("fastsa files 2","*.fa"),("text files","*.txt"),
                                          ("all files","*.*")))
## This function uses the Clustalw.exe in order to create .aln and .dnd files  
def clusterFunction():
    # Try block to check whether clusterw is installed and checks if file is right type
    try:
        clustalw_exe = r"programs\ClustalW2\clustalw2.exe"
        clustalw_cline = ClustalwCommandline(clustalw_exe, infile = fileName)
        assert os.path.isfile(clustalw_exe), "Clustal_W executable is missing or not found"
        # this writes the results of the alignments to the dedicated files as same name as fasta file, different ext.
        stdout, stderr= clustalw_cline()
    # Raise Error and show user the error message
    except:
        messagebox.showerror('Error', "Clustalw needs fasta file loaded to create .aln and .dnd files")

def seqManipulation():
    # Strings from entry boxes
    string1 = Seq(sequence1.get())
    string2 = Seq(sequence2.get())

    # Convert sequences to upper case if not
    string1 = string1.upper()
    string2 = string2.upper()
    try:
        # Get sequence complements
        seqComp1 = string1.complement()
        seqComp2 = string2.complement()

        # Get sequence reverse complements
        seqReverseComp1 = string1.reverse_complement()
        seqReverseComp2 = string2.reverse_complement()

        # Get sequence transcriptions
        seqTranscribe1 = string1.transcribe()
        seqTranscribe2 = string2.transcribe()

        # Get sequences Translations
        seqTranslate1 = string1.translate()
        seqTranslate2 = string2.translate()

        #Create message box to display to user all sequence manipulations
        messagebox.showinfo('Sequence Manipulation', "Original Sequence 1 = " + string1 + "\n\n"
                                               + "Original Sequence 2 = " + string2 + "\n\n"
                                               + "Sequence 1 complement = " + seqComp1 + "\n\n" 
                                               + "Sequence 2 complement = " + seqComp2 + "\n\n"
                                               + "Sequence 1 reverse complement = " + seqReverseComp1 + "\n\n"
                                               + "Sequence 2 reverse complement = " + seqReverseComp2 + "\n\n"
                                               + "Sequence 1 transcribed = " + seqTranscribe1 + "\n\n"
                                               + "Sequence 2 transcribed = " + seqTranscribe2 + "\n\n"
                                               + "Sequence 1 translated = " + seqTranslate1 + "\n\n"
                                               + "Sequence 2 translated = " + seqTranslate2 + "\n"
                                               )
    except:
        messagebox.showerror('Error', "Please enter right format for DNA Sequence")
##This function opens the READme when prompted by user
def readmeFunction():
    #Use operating system to open file right away
    try:
        os.startfile("files\READme.pdf")
    except:
        messagebox.showerror('Error', "READme.pdf was not found. Make sure it was not moved or deleted")
#FRAMES
instructionFrame = tk.LabelFrame(window, text = "Instructions", font =(font, 25))
buttonFrame = tk.Frame(window)
bannerFrame = tk.Frame(window)
entryFrame = tk.Frame(window)

#LABELS
## Instruction Label
instructions = tk.Label(instructionFrame,text = ("Welcome to the SJSU Fasta ToolBox\n\n " + 
                                            "Step 1 : Open a Fasta file you wish to work with\n\n" +
                                            " Step 2 : Then click Clustalw File Creator \n\n" +  
                                            " Step 3 : The aligned sequences can be viewed with the Align File button\n\n" +
                                            " Step 4 : Distance Matrix can be viewed with Create Distance Matrix button\n\n" +
                                            " Step 5 : Click Draw Tree button and select dnd file created by Clustalw, creating a Tree \n\n" +
                                            " Step 6 : Enter sequences in text entry areas and do desired functions with either button\n"
                                            ),justify = LEFT)
## Button Label
buttonTitle = tk.Label(buttonFrame, text="Program Buttons", font=(font,14))
bannerLabel = tk.Label(bannerFrame, image=banner)

## Sequence Labels
sequenceTitle = tk.Label(entryFrame, text="Sequences to Compare", font=(font,15))
seq1Label = tk.Label(entryFrame, text = "Sequence 1 :")
seq2Label = tk.Label(entryFrame, text = "Sequence 2 :")

#Entries
## Sequence Entries
sequence1 = tk.Entry(entryFrame, width=50)
sequence2 = tk.Entry(entryFrame,width=50)
gap = tk.Label(entryFrame, text= "  ")
gap2 = tk.Label(entryFrame, text = "             ")

#BUTTONS
## Buttom Frame Buttons
fileButton = tk.Button(buttonFrame,text="Open",command=openFile)
alignButton = tk.Button(buttonFrame, text = "Align File", command=alignFile)
matrixButton = tk.Button(buttonFrame, text="Create Distance Matrix", command=createMatrix)
percentIdentityButton = tk.Button(entryFrame, text="Percent Identity", command=percentIdentity)
clusterButton = tk.Button(buttonFrame, text="Clustalw File Creator",command=clusterFunction)
readMeButton = tk.Button(buttonFrame, text = "READme",command=readmeFunction)

## Enrty Frame Buttons
drawTreeButton = tk.Button(buttonFrame, text="Draw Tree", command=drawTree)
patternMatchButton = tk.Button(entryFrame, text="DNA Sequence Manipulation", command=seqManipulation)

# This is the section of code which creates the main window
window.geometry('650x650')
window.configure(background='palegoldenrod' )
window.title('SJSU Fasta ToolBox')

# Pack Frames
bannerFrame.pack(side = TOP)
instructionFrame.pack(side =TOP )
buttonFrame.pack(side = RIGHT,padx=40)
entryFrame.pack(side = LEFT,padx =10, ipadx=10)

# Pack Widgets
instructions.pack(pady=10)
bannerLabel.pack(side=TOP,fill=X, expand=1)
buttonTitle.pack(side = TOP,pady=5)
fileButton.pack(side = TOP,pady=5)
clusterButton.pack(side = TOP,pady=5)
alignButton.pack(side = TOP,pady=5)
matrixButton.pack(side = TOP,pady=5)
drawTreeButton.pack(side = TOP,pady=5)
readMeButton.pack(side=TOP, pady=5)

#Grid for EntryFrame
sequenceTitle.grid(row=0,column=0,columnspan =3, padx=10,pady=5)
seq1Label.grid(row=1,column=0,pady = 5)
sequence1.grid(row=1,column=1,pady=5)
seq2Label.grid(row=2,column=0,pady=5)
sequence2.grid(row=2,column=1,pady=5)

#Button Grids
percentIdentityButton.grid(row=3,column=0)
gap.grid(row=3,column=4)
gap2.grid(row=3,column=3)
patternMatchButton.grid(row=3,column=1)

#Main Loop to Run GUI
window.mainloop()
