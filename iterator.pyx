cdef extern from "stdlib.h":
    ctypedef unsigned long size_t
    void* malloc(size_t size)
    void free(void* mem)

from multiprocessing import Process, Queue, cpu_count #, heap
#from multiprocessing.sharedctypes import RawArray
from math import log

cdef class Iterator:
    cdef int W, H, max, strsize, cpus
    cdef double *real_axis, *imag_axis
    cdef unsigned char *log_scale, *image_string, *lo_count, *hi_count
    cdef BL, TR, param, queue, lo_escape, hi_escape
    
    def __cinit__(self, int W=500, int H=500, int max=255,
                  BL=-1-1j, TR=1+1j, param=0+0j):
        if W <= 0 or H <= 0:
            return
        self.log_scale = NULL
        self.max = 0
        strsize = W*H
        self.image_string = <unsigned char *>malloc(strsize*sizeof(char))
        self.lo_count = <unsigned char *>malloc(strsize*sizeof(char))
        self.hi_count = <unsigned char *>malloc(strsize*sizeof(char))
        self.real_axis = <double *>malloc(W*sizeof(double));
        self.imag_axis = <double *>malloc(H*sizeof(double));
        
    def __dealloc__(self):
        free(self.log_scale)
        free(self.image_string)
        free(self.lo_count)
        free(self.hi_count)
        free(self.real_axis)
        free(self.imag_axis)

    def __init__(self, int W=500, int H=500, int max=255,
                 BL=-1-1j, TR=1+1j, param=0+0j):
        if W <= 0 or H <= 0:
            raise ValueError('Image dimensions must be positive.')
        self.set_max(max)
        if not (self.image_string and self.lo_count and
                self.hi_count and self.real_axis and self.imag_axis):
            raise MemoryError
        self.W, self.H, self.max = W, H, max
        self.cpus = cpu_count()
        self.queue = Queue()
#        self.counts = RawArray('H', W*H)
        self.set(BL, TR, param)
        
    def __repr__(self):
        return 'Iterator: W=%s, max=%s, max=%s, BL=%s, TR=%s, param=%s'%(
            self.W, self.H, self.max, self.BL, self.TR, self.param)
    
    def build_axes(self):
        cdef double R, I, delta
        cdef int i
        R = self.BL.real
        delta = (self.TR.real - R)/self.W
        for i in range(self.W):
            self.real_axis[i] = R
            R += delta
        I = self.TR.imag # For graphics, the y axis points down
        delta = (self.BL.imag - I)/self.H
        for i in range(self.H):
            self.imag_axis[i] = I
            I += delta
        pass

    def set(self, BL, TR, param):
        self.BL, self.TR, self.param = BL, TR, param
        self.build_axes()

    def set_max(self, int max):
        cdef int n
        cdef unsigned char *log_scale
        if max < 1:
            raise ValueError('Max must be positive.')
        if max == self.max:
            return
        self.max = max
        log_scale = <unsigned char *>malloc((1+max)*sizeof(char))
        if log_scale == NULL:
            raise RuntimeError('Out of memory')
        free(self.log_scale)
        self.log_scale = log_scale
        log_scale[0] = 0
        for n in range(1, 1+max):
            log_scale[n] = min(n, int(1 + 254*log(n)/log(max)))
        
    def get_image(self):
        processes = []
        bands = [None]*self.cpus
        lows =  [None]*self.cpus
        highs = [None]*self.cpus
        # Start a worker process for each cpu.
        for n in range(self.cpus):
            P = Process(target=self.iterate, args=[n, self.cpus])
            processes.append(P)
            P.start()
        # Collect the results
        for n in range(self.cpus):
            k, image, low, high = self.queue.get()
            bands[k], lows[k], highs[k] = image, low, high
        # Join all the subprocesses - "probably good practice"
        for P in processes:
            P.join()
        # Merge the escape times
        self.lo_escape = ''.join(lows)
        self.hi_escape = ''.join(highs)
        # Merge the bands
        return ''.join(bands)

    def get_Z(self, m, n):
        try:
            return self.real_axis[m] + 1j*self.imag_axis[n]
        except IndexError:
            return None

    def get_escape(self, int m, int n):
        cdef i = n*self.W + m
        try:
            low = ord(self.lo_escape[i])
            high = ord(self.hi_escape[i])
            return low + (high << 8)
        except IndexError:
            return -1
        
    def iterate(self):
        # Subclasses override this method
        pass
    
cdef class Z_Iterator(Iterator):
    """
    Loops through all values in the specified band of the Z-rectangle,
    computing the number of steps before escape, using Z as the start.
    (Computes the Julia set for C = param.)
    """
    def iterate(self, int band=0, int num_bands=1):
        cdef int S, N, H0, H1, W, d, i, j, k, iterations, maxit=self.max+1
        cdef unsigned char* log_scale = self.log_scale
        cdef double R, I, RR, II
        cdef double Cr=self.param.real, Ci=self.param.imag
        cdef double *Zr = self.real_axis, *Zi=self.imag_axis
        cdef unsigned char* imgstr = self.image_string
        cdef unsigned char* low = self.lo_count
        cdef unsigned char* high = self.hi_count
#        cdef counts = self.counts
        d = self.H/num_bands
        H0 = band*d
        H1 = self.H if band == num_bands-1 else H0 + d
        W = self.W
        S = N = H0*W
        for j in range(H0, H1):
            for i in range(W):
                iterations = 0
                R, I = Zr[i], Zi[j]
                for k in range(maxit):
                    RR, II = R*R, I*I
                    if RR + II > 4.0:
                        iterations = k+1
                        break
                    I = 2*I*R + Ci
                    R = RR - II + Cr
                #counts[N] = iterations
                imgstr[N] = log_scale[iterations]
                low[N] = iterations & 255
                high[N] = iterations >> 8
                N += 1
        self.queue.put( (band,
                         bytes(imgstr[S:N]),
                         bytes(low[S:N]),
                         bytes(high[S:N]),
                         ) )
        
cdef class C_Iterator(Iterator):
    """
    Loops through all values in the specified band of the C-rectangle,
    computing the number of steps before the point specified as param
    escapes.  (Computes the Mandelbrot set, with param=0+0j.)
    """
    def iterate(self, int band=0, int num_bands=1):
        cdef int S, N, H0, H1, W, d, i, j, k, iterations, maxit=self.max+1
        cdef unsigned char *log_scale = self.log_scale
        cdef double R, I, RR, II
        cdef double R0=self.param.real, I0=self.param.imag
        cdef double *Cr=self.real_axis, *Ci=self.imag_axis
        cdef unsigned char* imgstr = self.image_string
        cdef unsigned char* low = self.lo_count
        cdef unsigned char* high = self.hi_count
        #cdef counts = self.counts
        # Actually, each process gets its own copy of self.image_string
        # but I am having each one write to its own part, in case I
        # decide it is worthwhile trying to share the string instead
        # of using a queue.
        d = self.H/num_bands
        H0 = band*d
        H1 = self.H if band == num_bands-1 else H0 + d
        W = self.W
        S = N = H0*W
        for j in range(H0, H1):
            for i in range(W):
                R, I = R0, I0
                iterations = 0
                for k in range(1, maxit):
                    RR, II = R*R, I*I
                    if RR + II > 4.0:
                        iterations = k
                        break
                    I = 2*I*R + Ci[j]
                    R = RR - II + Cr[i]
                #counts[N] = iterations
                imgstr[N] = self.log_scale[iterations]
                low[N] = iterations & 255
                high[N] = iterations >> 8
                N += 1
        self.queue.put( (band,
                         bytes(imgstr[S:N]),
                         bytes(low[S:N]),
                         bytes(high[S:N]),
                         ) )

def boundary(image, int W, int H):
    """
    Generates a black and white image of the boundary of a colored region.
    Arguments:
        image    the image as a string
        W        the width of the image
        H        the height of the image
    Returns an image which is the boundary of the region in which pixel
    values are 0.
    """
    cdef unsigned char *in_string = image 
    cdef unsigned char *out_string
    cdef int length, i, j, p
    length = W*H
    # Validate arguments
    if length != len(image):
        raise ValueError('Invalid image dimensions.')
    # Allocate memory
    out_string = <unsigned char *>malloc(length*sizeof(char))
    # Find the boundary
    for i in range(W):
        for j in range(H):
            p = i+j*W;
            out_string[p] = 255;
            if (in_string[p] > 0):
                if ( (i > 0 and in_string[p-1] == 0) or
                     (i < W and in_string[p+1] == 0) or
                     (j > 0 and in_string[p-W] == 0) or
                     (j < H and in_string[p+W] == 0) ):
                    out_string[p] = 0;
    result = bytes(out_string[:length])
    # Free memory
    free(out_string)
    # Return a string (as bytes)
    return result

def boxcount(image, int W, int H, int max = 255):
    """
    Inputs an image as a string. Counts how many boxes of size 2^i
    meet the image.  Returns a list of integers, where the ith integer
    is the number of boxes of size 2^i that meet the image.  The
    maximum box size is the largest power of 2 which is smaller than
    min(W,H)/4.  A box is deemed to meet the image if it contains a
    pixel of value >= max.

    Arguments:
        image    the image as a string
        W        the width of the image
        H        the height of the image
        max      threshold
    """
    cdef unsigned char *in_string = <unsigned char *>image
    cdef char *buf
    cdef int i, j, m, n, length, box_size
    cdef long count
    length = W*H
    # Validate arguments
    if length != len(image):
        raise ValueError('Invalid image dimensions.')
    if max < 0 or max > 255:
        raise ValueError('The value of max must lie between 0 and 255.')
    # Allocate memory
    buf = <char *>malloc(length*sizeof(char));
    # Compute dimension
    result = []
    for i in range(length):
        buf[i] = 1 if in_string[i] < max else 0
    box_size = (W if W < H else H) >> 2
    while box_size > 1:
        i, j = 0, 0
        for m in range(H, 0, -2):
            for n in range(W, 0, -2):
                if n > 1:
                    buf[i] = buf[j] | buf[j+1]
                    if m > 1:
                        buf[i] |= (buf[j+W] | buf[j+W+1])
                    j += 2
                else:
                    buf[i] = buf[j]
                    if m > 1:
                        buf[i] |= buf[j+W]
                    j += 1
                i += 1
            j += W;
        W = (W+1)/2
        H = (H+1)/2
        count = 0
        while i > 0:
            i -= 1
            count += buf[i]
        result.append(count)
        box_size = box_size >> 1;
    free(buf);
    return result;
