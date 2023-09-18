# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 16:16:12 2022

@author: Nico Kubetschek & Stuart Goldie
"""

import tkinter as tk

 
material = ""
methode = ""


def function_name1():
    global material, methode 
    material = "MoS2"

def function_name2():
    global material, methode 
    material = "MoSe2"
    
def function_name3():
    global material, methode 
    material = "WS2"
    
def function_name4():
    global material, methode 
    material = "WSe2"

def function_name5():
    global material, methode
    material = "NiPS3"

def function_name6():
    global material, methode
    material = "FePS3"

def function_name7():
    global material, methode
    material = "PtSe2"

def function_name8():
    global material, methode 
    methode = "Abs"
    
def function_name9():
    global material, methode 
    methode = "Ext"
    
def function_quite():
    root.destroy()
    root.quit()
        
# creating root window
root = tk.Tk()
 
# root window title and dimension
root.title("Choose the material and the method of your spectra")
root.geometry("300x550")
 

# creating button
btn1 = tk.Button(root, text="MoS2", command = function_name1, activebackground="green",bg = "lightgreen",height=2,font=12,width=18,relief="groove")
btn1.pack()

btn2 = tk.Button(root, text="MoSe2", command = function_name2,activebackground="green",bg = "lightgreen",height=2,font=12,width=18,relief="groove")
btn2.pack()

btn3 = tk.Button(root, text="WS2", command = function_name3,activebackground="green",bg = "lightgreen",height=2,font=12,width=18,relief="groove")
btn3.pack()

btn4 = tk.Button(root, text="WSe2", command = function_name4,activebackground="green",bg = "lightgreen",height=2,font=12,width=18,relief="groove")
btn4.pack()

btn5 = tk.Button(root, text="NiPS3", command = function_name5,activebackground="green",bg = "lightgreen",height=2,font=12,width=18,relief="groove")
btn5.pack()

btn6 = tk.Button(root, text="FePS3", command = function_name6,activebackground="green",bg = "lightgreen",height=2,font=12,width=18,relief="groove")
btn6.pack()

btn7 = tk.Button(root, text="PtSe2", command = function_name7,activebackground="green",bg = "lightgreen",height=2,font=12,width=18,relief="groove")
btn7.pack()

btn8 = tk.Button(root, text="Absorption", command = function_name8,activebackground="green",bg = "steelblue",height=2,font=12,width=18,relief="groove")
btn8.pack()

btn9 = tk.Button(root, text="Extiction", command = function_name9,activebackground="green",bg = "steelblue",height=2,font=12,width=18,relief="groove")
btn9.pack()

btn10 = tk.Button(root, text="Load file", command=function_quite,bg = "red",height=2,font=12,width=18,relief="groove")
btn10.pack()
root.mainloop()

print("Material:\t",material)
print("Method:\t\t",methode)



