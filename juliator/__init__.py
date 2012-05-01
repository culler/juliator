import Tkinter
import Image
import ImageTk
import os, sys, math, cmath, colorsys
from iterator import iterate, boundary, boxcount, C_Iterator, Z_Iterator
from multiprocessing import Process, Queue, cpu_count

def vibgyor(size=256):
    result = []
    M = float(size - 1)*1.125
    for n in range(size-1):
        r, g, b = colorsys.hsv_to_rgb(n/M, 0.85, 1.0)
        result += [int(r*255), int(g*255), int(b*255)]
    result += [0,0,0]
    return result

class Viewer(Tkinter.Frame):
    """
    Base class for panes that display colormap images.
    """
    def __init__(self, parent, W=500, H=500, width=4.0, center=0.0+0.0j,
                 palette=vibgyor()):
        Tkinter.Frame.__init__(self, parent)
        self.parent = parent
        self.W, self.H, self.aspect = W, H, float(H)/float(W)
        self.image = None
        self.dragging = 0
        self.box = None
        self.dot = None
        self.palette = palette
        self.canvas = Tkinter.Canvas(self, width=W, height=H)
        self.canvas.bind('<Motion>', self.motion)
        self.canvas.bind('<Leave>', self.leave)
        self.canvas.bind('<Button-1>', self.mouse_down)
        self.canvas.bind('<ButtonRelease-1>', self.mouse_up)
        self.canvas.bind('<Button-3>', self.zoom_out)
        self.canvas.bind('<Key>', self.keypress)
        self.canvas.pack()
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

    # Subclasses should have an iterator attribute
    def display_image(self):
        self.imagestring = self.iterator.get_image()
        self.image = Image.fromstring('P', (self.W, self.H),
                                      self.imagestring)
        self.image.putpalette(self.palette)
        self.bitmap = ImageTk.PhotoImage(self.image)
        image_name = self.canvas.create_image(0,0,
                                           anchor=Tkinter.NW,
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

    def leave(self,event):
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

    def set_dot(self, z):
        self.clear_dot()
        x = int((0.5 + (z.real - self.center.real)/self.width)*self.W)
        y = int((0.5 + (self.center.imag - z.imag)/(self.width*self.aspect))*self.H)
        self.dot = self.canvas.create_oval(x-2, y-2, x+2, y+2, fill='white')

    def clear_dot(self):
        if self.dot:
            self.canvas.delete(self.dot)

class Mandelbrot(Viewer):
    """
    Viewer for the Mandelbrot set.
    """
    def __init__(self, parent, W=500, H=500, width=2.75,
                 center=-0.75+0.0j, palette=vibgyor()):
        Viewer.__init__(self, parent, W, H, width, center, palette)
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
#        self.window.title('Mandelbrot set -- center = %f + %fi, width = %f'%
#                          (self.center.real, self.center.imag, self.width))
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
        Viewer.__init__(self, parent, W, H, width, center, palette)
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
        self.set_title()
        height = 1j*self.width*float(self.H)/float(self.W)
        z0 = self.center - (self.width + height)/2
        z1 = self.center + (self.width + height)/2
        self.iterator.set(z0, z1, self.c)
        self.display_image()

    def set_title(self):
        pass
#        self.window.title(
#            '%s Julia set: c=%.3f+%.3fi, center=%.3f+%.3fi, width=%e'%
#            ('Filled' if self.filled else '',
#             self.c.real, self.c.imag, self.center.real,
#             self.center.imag, self.width))

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
        self.canvas.create_image(0,0, anchor=Tkinter.NW, image=self.bitmap)
        if not self.filled:
            self.canvas.create_rectangle(20, self.H-14, 300, self.H,
                                         fill='white', outline = '')
            self.canvas.create_text(
                20, self.H-14, fill='red',
                anchor=Tkinter.NW,  
                text= 'Box dimension is approximately %.3f'%(
                self.box_dimension(bdry) )
                )
        self.set_title()

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
        if Tkinter._default_root:
            self.window = Tkinter.Toplevel(Tkinter._default_root)
        else:
            self.window = Tkinter.Tk()
        self.mandelbrot = Mandelbrot(self.window)
        self.mandelbrot.grid(row=0, column=0)
        self.mandelbrot.julia.grid(row=0, column=1)
