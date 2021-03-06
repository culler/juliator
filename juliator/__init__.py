# -*- coding: utf-8 -*-

# Copyright© 2012-2017 by Marc Culler and others.
# This file is part of Juliator.
#
# Juliator is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Juliator is distributed in the hope that it will be useful
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Juliator.  If not, see <http://www.gnu.org/licenses/>.

try:
    import Tkinter as Tk_
    import ttk
    import tkFileDialog
    import tkMessageBox
except ImportError: # Python3
    import tkinter as Tk_
    from tkinter import ttk
    import tkinter.filedialog as tkFileDialog
    import tkinter.messagebox as tkMessageBox
try:
    from PIL import Image, ImageTk
except ImportError:
    import Image, ImageTk
import os, sys, math, cmath, colorsys, time
from iterator import C_Iterator, Z_Iterator, boundary, boxcount

if sys.platform == 'darwin':
    WindowBG = 'SystemDialogBackgroundActive'
    GroupBG = '#e0e0e0'
    BrowserBG = '#a8a8a8'
elif sys.platform == 'linux2':
    WindowBG = '#e0e0e0'
    GroupBG = '#e0e0e0'
    BrowserBG = '#a8a8a8'

cx_format = '{0.real:0<10.7g}{0.imag:0<+10.7g} i'
cx_format = '{0.real:<.5g}{0.imag:<+.5g} i'

class Chooser(ttk.Menubutton):
    """
    The ttk Menubutton is useless out of the box.
    """
    def __init__(self, master, choices=[''], command=lambda x:None):
        ttk.Menubutton.__init__(self, master)
        self.command=command
        self.variable = Tk_.StringVar(master)
        self.variable.set(choices[0])
        self.config(textvariable=self.variable)
        self.menu = Tk_.Menu(master)
        for choice in choices:
            self.menu.add_command(
                label=choice,
                command=Tk_._setit(self.variable, choice, command)
                )
        self['menu'] = self.menu
        self['width'] = max([len(c) for c in choices])

    def get(self):
        return self.variable.get()

    def set(self, value):
        self.variable.set(str(value))
        
def vibgyor(size=256):
    result = [0,0,0]
    M = float(size - 1)*1.125
    for n in range(size-1):
        r, g, b = colorsys.hsv_to_rgb(n/M, 0.85, 1.0)
        result += [int(r*255), int(g*255), int(b*255)]
    return result

class Viewer(ttk.LabelFrame):
    """
    Base class for panes that display colormap images.
    """
    def __init__(self, parent, W=500, H=500, width=4.0, center=0.0+0.0j,
                 palette=vibgyor(), text='Viewer'):
        ttk.LabelFrame.__init__(self, parent, text=text)
        self.parent = parent
        self.W, self.H, self.aspect = W, H, float(H)/float(W)
        self.image = None
        self.dragging = 0
        self.ctrl_dragging = 0
        self.box = None
        self.palette = palette
        self.canvas = Tk_.Canvas(self, width=W, height=H, bd=0, bg=GroupBG,
                                 highlightbackground=GroupBG)
        self.canvas.bind('<Motion>', self.motion)
        self.canvas.bind('<Leave>', self.leave)
        self.canvas.bind('<Button-1>', self.click)
        self.canvas.bind('<Control-Button-1>', self.ctrl_click)
        self.canvas.bind('<Button-2>', self.ctrl_click)
        self.canvas.bind('<ButtonRelease-1>', self.unclick)
        self.canvas.bind('<Control-ButtonRelease-1>', self.unclick)
        self.canvas.bind('<ButtonRelease-2>', self.unclick)
        self.canvas.bind('<Key>', self.keypress)
        self.canvas.pack()
        self.controlpanel = Tk_.Frame(self, bg=GroupBG, height=30)
        label = Tk_.Label(self.controlpanel, text='  Cutoff: ', bg=GroupBG)
        label.pack(side=Tk_.LEFT)
        self.max_choice = Chooser(self.controlpanel,
                                  [str(256*2**n) for n in range(6)],
                                  command=self.set_max)
        self.max_choice.pack(side=Tk_.LEFT)
        self.show_center = Tk_.StringVar(self)
        label = Tk_.Label(self.controlpanel, text='  Center:', bg=GroupBG)
        label.pack(side=Tk_.LEFT)
        label = Tk_.Label(self.controlpanel,
                          textvar=self.show_center, bg=GroupBG)
        label.pack(side=Tk_.LEFT)
        self.show_width = Tk_.StringVar(self)
        label = Tk_.Label(self.controlpanel, text='  Width:', bg=GroupBG)
        label.pack(side=Tk_.LEFT)
        label = Tk_.Label(self.controlpanel,
                          textvar=self.show_width, bg=GroupBG)
        label.pack(side=Tk_.LEFT)
        self.controlpanel.pack(expand=True, fill=Tk_.X)
        self.mouse_location = Tk_.StringVar(self)
        self.escape_time = Tk_.StringVar(self)
        self.image_name = None
        
    def keypress(self, event):
        if event.char == '-' or event.char == '_':
            self.zoom_out(event)
        if event.char == '+' or event.char == '=':
            self.zoom_in(event)

    # Subclasses should override this, and call from __init__.
    def set(self, width=None, center=None):
        pass

    # Subclasses must provide an iterator attribute.
    def display_image(self):
        #start = time.time()
        self.imagestring = self.iterator.get_image()
        #print(time.time() - start)
        self.image = Image.frombytes('P', (self.W, self.H),
                                      self.imagestring)
        self.image.putpalette(self.palette)
        self.bitmap = ImageTk.PhotoImage(self.image)
        image_name = self.canvas.create_image(1,1,
                                           anchor=Tk_.NW,
                                           image=self.bitmap)
        self.canvas.delete(self.image_name)
        self.image_name = image_name
        
    def save_image(self):
        filename = tkFileDialog.asksaveasfilename(
            initialfile='unnamed.png',
            )
        if filename:
            try:
                self.image.save(filename)
            except KeyError:
                extension = os.path.splitext(filename)[-1]
                tkMessageBox.showwarning(
                    'Save image',
                    'The image file extension %s is unknown.'%extension
                    )
            except IOError:
                tkMessageBox.showwarning(
                    'Save image',
                    'Failed to save image file.'
                    )
    
    def set_max(self, newmax):
        self.iterator.set_max(int(newmax))
        self.set(self.width, self.center)
        
    def motion(self, event):
        if self.dragging == 1:
            if self.box:
                self.canvas.delete(self.box)
            delta_x = event.x - self.x
            delta_y = int(self.aspect*delta_x)
            x0, x1 = self.x - delta_x, self.x + delta_x
            y0, y1 = self.y - delta_y, self.y + delta_y
            self.box = self.canvas.create_line(x0,y0,x1,y0,x1,y1,x0,y1,x0,y0,
                                               fill='white')
        elif self.ctrl_dragging:
            self.ctrl_motion(event)
        else:
            z = self.iterator.get_Z(event.x, event.y)
            escape = self.iterator.get_escape(event.x, event.y)
            if escape == 0:
                escape = '> cutoff'
            if z:
                self.mouse_location.set(cx_format.format(z))
                self.escape_time.set('Escape time: %s'%escape)
            else:
                self.mouse_location.set('')
                self.escape_time.set('')

    def ctrl_motion(self, event):
        pass
    
    def leave(self,event):
        self.mouse_location.set('')
        self.escape_time.set('')
        if self.box:
            self.canvas.delete(self.box)

    def click(self, event):
        self.canvas.focus_set()
        self.x, self.y = event.x, event.y
        self.dragging = 1

    def ctrl_click(self, event):
        self.canvas.focus_set()
        self.x, self.y = event.x, event.y
        if not self.dragging:
            self.ctrl_dragging = 1
            
    def unclick(self, event):
        if self.ctrl_dragging:
            self.end_ctrl_drag(event)
        if self.dragging:
            self.dragging = 0
            if self.box:
                self.canvas.delete(self.box)
            if abs(self.x - event.x) > 2:
                delta_x = self.width*(self.x - self.W/2)/self.W
                delta_y = self.aspect*self.width*(self.H/2 - self.y)/self.H
                newcenter = self.center + complex(delta_x,delta_y)
                newwidth = abs(2*self.width*float(self.x - event.x)/self.W) 
                self.set(center=newcenter, width=newwidth)

    def end_ctrl_drag(self, event):
        self.ctrl_dragging = 0
        
    def zoom_out(self, event):
        self.set(width=2*self.width)

    def zoom_in(self, event):
        self.set(width=self.width/2)

class Mandelbrot(Viewer):
    """
    Viewer for the Mandelbrot set.
    """
    def __init__(self, parent, W=500, H=500, width=2.75,
                 center=-0.75+0.0j, palette=vibgyor()):
        Viewer.__init__(self, parent, W, H, width, center, palette,
                        text='Mandelbrot Set')
        self.julia = Julia(self.parent, c=center)
        self.save_max = 256
        self.marker = ''
        self.iterator = C_Iterator(self.W, self.H, 255)
        self.set(width, center)
        
    def set(self, width=None, center=None):
        if width != None: self.width = width
        if center != None: self.center = center
        height = 1j*self.width*self.aspect
        c0 = self.center - (self.width + height)/2
        c1 = self.center + (self.width + height)/2
        self.iterator.set(c0, c1, 0+0j)
        self.display_image()
        self.show_center.set(cx_format.format(self.center))
        self.show_width.set('%.6g'%self.width)
        self.show_julia(self.W/2, self.H/2)

    def ctrl_click(self, event):
        if not self.dragging:
            self.save_max = int(self.julia.max_choice.get())
            self.julia.max_choice.set(256)
            self.julia.iterator.set_max(256)
            self.ctrl_dragging = 1
        self.show_julia(event.x, event.y)

    def ctrl_motion(self, event):
        self.show_julia(event.x, event.y)

    def end_ctrl_drag(self, event):
        self.ctrl_dragging = 0
        self.julia.max_choice.set(self.save_max)
        self.julia.iterator.set_max(self.save_max)
        self.show_julia(event.x, event.y)
        
    def show_julia(self, x, y):
        delta_x = self.width*(x - self.W/2)/self.W
        delta_y = self.aspect*self.width*(self.H/2 - y)/self.H
        center = self.center + complex(delta_x,delta_y)
        self.canvas.delete(self.marker)
        self.marker = self.canvas.create_oval(x-4, y-4, x+4, y+4,
                                              outline='white', width=2)
        self.julia.filled = True
        self.julia.set(c=center)
        
class Julia(Viewer):
    """
    Viewer for Julia sets.
    """
    def __init__(self, parent, W=500, H=500, width=4.0,
                 center=0.0+0.0j, c=0.0+0.0j, palette=vibgyor()):
        Viewer.__init__(self, parent, W, H, width, center, palette,
                 text='Julia Set')
        self.canvas.bind('<Button-2>', self.toggle_fill)
        if sys.platform == 'darwin':
            self.canvas.bind('<Control-Button-1>', self.toggle_fill)
            self.canvas.bind('<Control-ButtonRelease-1>', lambda event:None)
        self.iterator = Z_Iterator(self.W, self.H, 255)
        self.set(width, center, c)
        
    def set(self, width=None, center=None, c=None):
        if width != None: self.width = width
        if center != None: self.center = center
        if c != None:
            self.c = c
        self.filled = True
        height = 1j*self.width*float(self.H)/float(self.W)
        z0 = self.center - (self.width + height)/2
        z1 = self.center + (self.width + height)/2
        self.iterator.set(z0, z1, self.c)
        self.display_image()
        self.show_center.set(cx_format.format(self.center))
        self.show_width.set('%.6g'%self.width)

    def toggle_fill(self, event):
        if self.filled:
            self.filled = False
            bdry = boundary(self.imagestring, self.W, self.H)
            self.image = Image.frombytes('L', (self.W, self.H), bdry)
        else:
            self.filled = True
            self.image = Image.frombytes(
                'P', (self.W, self.H), self.imagestring)
            self.image.putpalette(self.palette)

        self.bitmap = ImageTk.PhotoImage(self.image)
        self.canvas.create_image(0,0, anchor=Tk_.NW, image=self.bitmap)
        if not self.filled:
            self.canvas.create_rectangle(20, self.H-14, 300, self.H,
                                         fill='white', outline = '')
            self.canvas.create_text(
                20, self.H-14, fill='red',
                anchor=Tk_.NW,  
                text= 'Box dimension is approximately %.3f'%(
                self.box_dimension(bdry) )
                )

    def box_dimension(self, boundary):
        count = boxcount(boundary, self.W, self.H)
        data = []
        for n in range(len(count)):
            try:
                data.append( (-(n+1)*math.log(2), math.log(count[n]) ) )
            except ValueError:
                pass
        X, Y, XY, XX = 0.0, 0.0, 0.0, 0.0
        if len(data) == 0:
            return 0.0
        n = float(len(data))
        for point in data:
            x, y = point
            X += x
            Y += y
            XY += x*y
            XX += x*x
        return (n*XY - X*Y)/(n*XX - X*X)

class Juliator:
    def __init__(self):
        if Tk_._default_root:
            self.window = Tk_.Toplevel(Tk_._default_root)
        else:
            self.window = Tk_.Tk()
        self.window.title('Juliator')
        self.top = Tk_.Frame(self.window, bg=WindowBG)
        self.mandelbrot = mandelbrot = Mandelbrot(self.top)
        self.julia = julia = mandelbrot.julia
        mandelbrot.grid(row=0, column=0, padx=4, pady=4)
        julia.grid(row=0, column=1, padx=4, pady=4)
        self.top.grid(row=0, column=0)
        self.separator = ttk.Separator(self.window, orient=Tk_.HORIZONTAL)
        self.separator.grid(row=1, column=0, sticky=Tk_.EW)
        self.bottom = Tk_.Frame(self.window, bg='white', height=20)
        self.bottom.columnconfigure(1, weight=1)
        self.bottom.columnconfigure(3, weight=1)
        self.mandel_where = Tk_.Label(self.bottom, width=20, anchor=Tk_.W,
                                      bg='white',
                                      textvar=mandelbrot.mouse_location)
        self.mandel_where.grid(sticky=Tk_.W, padx=10, row=0, column=0)
        self.mandel_when = Tk_.Label(self.bottom, width=20, anchor=Tk_.W,
                                      bg='white',
                                      textvar=mandelbrot.escape_time)
        self.mandel_when.grid(sticky=Tk_.W, padx=10, row=0, column=1)
        self.julia_where = Tk_.Label(self.bottom, width=20, anchor=Tk_.W,
                                     bg='white',
                                     textvar=julia.mouse_location)
        self.julia_where.grid(sticky=Tk_.W, padx=10, row=0, column=2)
        self.julia_when = Tk_.Label(self.bottom, width=20, anchor=Tk_.W,
                                     bg='white',
                                     textvar=julia.escape_time)
        self.julia_when.grid(sticky=Tk_.W, padx=10, row=0, column=3)
        self.bottom.grid(row=2, column=0, sticky=Tk_.NSEW)
        menubar = Tk_.Menu(self.window)
        file_menu = Tk_.Menu(menubar, tearoff=0)
        print_menu = Tk_.Menu(menubar, tearoff=0)
        print_menu.add_command(label='Mandelbrot',
                               command=self.mandelbrot.save_image)
        print_menu.add_command(label='Julia',
                               command=self.julia.save_image)
        file_menu.add_cascade(label='Save Image', menu=print_menu)
        file_menu.add_separator()
        menubar.add_cascade(label='File', menu=file_menu)
        self.window.config(menu=menubar)

