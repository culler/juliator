import Tkinter as Tk_
import ttk
import Image
import ImageTk
import os, sys, math, cmath, colorsys, time
from iterator import iterate, boundary, boxcount, C_Iterator, Z_Iterator
from multiprocessing import Process, Queue, cpu_count

if sys.platform == 'darwin':
    WindowBG = 'SystemDialogBackgroundActive'
    GroupBG = '#e0e0e0'
    BrowserBG = '#a8a8a8'
elif sys.platform == 'linux2':
    WindowBG = '#e0e0e0'
    GroupBG = '#e0e0e0'
    BrowserBG = '#a8a8a8'
    
class Chooser(ttk.Menubutton):
    def __init__(self, master, choices=[''], command=lambda x:None):
        ttk.Menubutton.__init__(self, master)
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
                 
def vibgyor(size=256):
    result = []
    M = float(size - 1)*1.125
    for n in range(size-1):
        r, g, b = colorsys.hsv_to_rgb(n/M, 0.85, 1.0)
        result += [int(r*255), int(g*255), int(b*255)]
    result += [0,0,0]
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
        self.box = None
        self.dot = None
        self.palette = palette
        self.canvas = Tk_.Canvas(self, width=W, height=H, bd=0, bg=GroupBG,
                                 highlightbackground=GroupBG)
        self.canvas.bind('<Motion>', self.motion)
        self.canvas.bind('<Leave>', self.leave)
        self.canvas.bind('<Button-1>', self.mouse_down)
        self.canvas.bind('<ButtonRelease-1>', self.mouse_up)
        self.canvas.bind('<Button-3>', self.zoom_out)
        self.canvas.bind('<Key>', self.keypress)
        self.canvas.pack()
        self.controlpanel = Tk_.Frame(self, bg=GroupBG, height=30)
        self.max_choice = Chooser(self.controlpanel,
                                  [str(256*2**n) for n in range(6)])
        self.max_choice.grid(row=0, column=3)
        self.controlpanel.pack(expand=True, fill=Tk_.X)
        self.mouse_location = Tk_.StringVar(self)
        self.image_name = None
        
    def keypress(self, event):
        if event.char == '-' or event.char == '_':
            self.zoom_out(event)
        if event.char == '+' or event.char == '=':
            self.zoom_in(event)

    # Subclasses should override this, and call from __init__.
    def set(self, width=None, center=None):
        if width != None: self.width = width
        if center != None: self.center = center

    # Subclasses should provide an iterator attribute
    def display_image(self):
        start = time.time()
        self.imagestring = self.iterator.get_image()
        print time.time() - start
        self.image = Image.fromstring('P', (self.W, self.H),
                                      self.imagestring)
        self.image.putpalette(self.palette)
        self.bitmap = ImageTk.PhotoImage(self.image)
        image_name = self.canvas.create_image(1,1,
                                           anchor=Tk_.NW,
                                           image=self.bitmap)
        self.canvas.delete(self.image_name)
        self.image_name = image_name
        
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
        else:
            z = self.iterator.get_Z(event.x, event.y)
            zstr = '{0.real:0<10.7g}{0.imag:0<+10.7g} i'.format(z)
            if z:
                self.mouse_location.set(zstr)
            else:
                self.mouse_location.set('')

    def leave(self,event):
        self.mouse_location.set('')
        if self.box:
            self.canvas.delete(self.box)

    def mouse_down(self, event):
        self.canvas.focus_set()
        self.x, self.y = event.x, event.y
        self.dragging = 1
        
    def mouse_up(self, event):
        self.dragging = 0
        if self.box:
            self.canvas.delete(self.box)
        if abs(self.x - event.x) > 2:
            delta_x = self.width*(self.x - self.W/2)/self.W
            delta_y = self.aspect*self.width*(self.H/2 - self.y)/self.H
            newcenter = self.center + complex(delta_x,delta_y)
            newwidth = abs(2*self.width*float(self.x - event.x)/self.W) 
            self.set(center=newcenter, width=newwidth)

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
        self.marker = ''
        self.canvas.bind('<Button-2>', self.show_julia)
        if sys.platform == 'darwin':
            self.canvas.bind('<Control-Button-1>', self.show_julia)
            self.canvas.bind('<Control-ButtonRelease-1>', lambda event:None)
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
        self.show_julia()

    def show_julia(self, event=None, center = 0+0j):
        if event:
            x, y = event.x, event.y
        else:
            x, y = self.W/2, self.H/2
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

    def toggle_fill(self, event):
        if self.filled:
            self.filled = False
            bdry = boundary(self.imagestring, self.W, self.H, 255)
            self.image = Image.fromstring('L', (self.W, self.H), bdry)
        else:
            self.filled = True
            self.image = Image.fromstring(
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
        count = boxcount(boundary, self.W, self.H, 255)
        data = []
        for i in range(len(count)):
            try:
                data.append( (-(i+1)*math.log(2), math.log(count[i]) ) )
            except OverflowError:
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
        self.bottom.columnconfigure(0, weight=1)
        self.bottom.columnconfigure(1, weight=1)
        self.mandel_where = Tk_.Label(self.bottom, width=40, anchor=Tk_.W,
                                      bg='white',
                                      textvar=mandelbrot.mouse_location)
        self.mandel_where.grid(sticky=Tk_.W, padx=10, row=0, column=0)
        self.julia_where = Tk_.Label(self.bottom, width=40, anchor=Tk_.W,
                                     bg='white',
                                     textvar=julia.mouse_location)
        self.julia_where.grid(sticky=Tk_.W, padx=10, row=0, column=1)
        self.bottom.grid(row=2, column=0, sticky=Tk_.NSEW)
