# Malloc
ctypedef unsigned long size_t
ctypedef long off_t

cdef extern from "stdlib.h":
    void* malloc(size_t size)
    void free(void* mem)

# Shared memory allocation
IF UNAME_SYSNAME == 'Windows':
    cdef void* sh_malloc(size_t size):
        cdef HANDLE fm
        fm = CreateFileMapping(INVALID_HANDLE_VALUE,NULL,PAGE_READWRITE,0,size,NULL)
        return MapViewOfFile(fm, FILE_MAP_ALL_ACCESS,0, 0, size)
        return NULL

    cdef void sh_free(void* ptr, size_t size):
        UnmapViewOfFile(ptr)
        CloseHandle(fm)

ELSE:
    cdef extern from "sys/mman.h":
        void* mmap(void *addr, size_t len, int prot, int flags, int fd,
                   off_t offset)
        int munmap(void *addr, size_t len)
        cdef enum:
            PROT_NONE
            PROT_READ
            PROT_WRITE
            PROT_EXEC
            MAP_ANON
            MAP_FILE
            MAP_SHARED
            MAP_PRIVATE
            MAP_FIXED

    cdef void* sh_malloc(size_t size):
        return mmap(NULL, size, PROT_READ|PROT_WRITE, MAP_ANON|MAP_SHARED, -1, 0)

    cdef int sh_free(void* ptr, size_t size):
        return munmap(ptr, size)

    cdef class SharedMemory:
        cdef size_t size
        cdef void* pointer
        def __cinit__(self, size_t size):
            self.size = size
            self.pointer = sh_malloc(size)
            
        def __dealloc__(self):
            sh_free(self.pointer, self.size)

from multiprocessing import Process, cpu_count
from math import log

cdef class Iterator:
    cdef int W, H, max, strsize, cpus
    cdef double* real_axis
    cdef double* imag_axis
    cdef unsigned char* log_scale
    cdef unsigned char* image_string
    cdef unsigned short *counts
    cdef BL, TR, param
    cdef SharedMemory image_memory, counts_memory
    
    def __cinit__(self, int W=500, int H=500, int max=255,
                  BL=-1-1j, TR=1+1j, param=0+0j):
        if W <= 0 or H <= 0:
            return
        self.log_scale = NULL
        self.max = 0
        self.strsize = strsize = W*H
        # Use shared memory, so subprocesses can change it.
        self.image_memory = SharedMemory(strsize*sizeof(char))
        self.image_string = <unsigned char *>self.image_memory.pointer
        self.counts_memory = SharedMemory(strsize*sizeof(short))
        self.counts = <unsigned short *>self.counts_memory.pointer
        # These are not changed by subprocesses, so malloc is fine.
        self.real_axis = <double *>malloc(W*sizeof(double));
        self.imag_axis = <double *>malloc(H*sizeof(double));
        
    def __dealloc__(self):
        free(self.log_scale)
        free(self.real_axis)
        free(self.imag_axis)

    def __init__(self, int W=500, int H=500, int max=255,
                 BL=-1-1j, TR=1+1j, param=0+0j):
        if W <= 0 or H <= 0:
            raise ValueError('Image dimensions must be positive.')
        self.set_max(max)
        if not (self.image_string and self.counts and
                self.real_axis and self.imag_axis):
            raise MemoryError
        self.W, self.H, self.max = W, H, max
        self.cpus = cpu_count()
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
        I = self.TR.imag
        # For graphics, the y axis points down
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
        if not log_scale:
            raise MemoryError
        free(self.log_scale)
        self.log_scale = log_scale
        log_scale[0] = 0
        for n in range(1, 1+max):
            log_scale[n] = min(n, int(1 + 254*log(n)/log(max)))
        
    def get_image(self):
        processes = []
        # Start a worker process for each cpu.
        for n in range(self.cpus):
            P = Process(target=self.iterate, args=[n, self.cpus])
            processes.append(P)
            P.start()
        # Join all the subprocesses - "probably good practice"
        for P in processes:
            P.join()
        return self.image_string[:self.strsize]

    def get_Z(self, m, n):
        try:
            return self.real_axis[m] + 1j*self.imag_axis[n]
        except IndexError:
            return None

    def get_escape(self, int m, int n):
        try:
            return self.counts[n*self.W + m]
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
        cdef double R, I, RR, II, RI
        cdef double Cr=self.param.real, Ci=self.param.imag
        cdef double* Zr = self.real_axis
        cdef double* Zi = self.imag_axis
        cdef unsigned char* imgstr = self.image_string
        cdef unsigned short* counts = self.counts

        d = self.H/num_bands
        H0 = band*d
        H1 = self.H if band == num_bands-1 else H0 + d
        W = self.W
        S = N = H0*W
        # For speed -- no python calls in this loop!
        for j in range(H0, H1):
            for i in range(W):
                iterations = 0
                R, I = Zr[i], Zi[j]
                for k in range(maxit):
                    RR, II = R*R, I*I
                    if RR + II > 4.0:
                        iterations = k+1
                        break
                    RI = R*I
                    I = RI + RI + Ci
                    R = RR - II + Cr
                counts[N] = iterations
                imgstr[N] = log_scale[iterations]
                N += 1
        
cdef class C_Iterator(Iterator):
    """
    Loops through all values in the specified band of the C-rectangle,
    computing the number of steps before the point specified as param
    escapes.  (Computes the Mandelbrot set, with param=0+0j.)
    """
    def iterate(self, int band=0, int num_bands=1):
        cdef int S, N, H0, H1, W, d, i, j, k, iterations, maxit=self.max+1
        cdef unsigned char *log_scale = self.log_scale
        cdef double R, I, RR, II, RI
        cdef double R0 = self.param.real, I0 = self.param.imag
        cdef double* Cr = self.real_axis
        cdef double* Ci = self.imag_axis
        cdef unsigned char* imgstr = self.image_string
        cdef unsigned short* counts = self.counts

        d = self.H/num_bands
        H0 = band*d
        H1 = self.H if band == num_bands-1 else H0 + d
        W = self.W
        S = N = H0*W
        # For speed -- no python calls in this loop!
        for j in range(H0, H1):
            for i in range(W):
                R, I = R0, I0
                iterations = 0
                for k in range(1, maxit):
                    RR, II = R*R, I*I
                    if RR + II > 4.0:
                        iterations = k
                        break
                    RI = R*I
                    I = RI + RI + Ci[j]
                    R = RR - II + Cr[i]
                counts[N] = iterations
                imgstr[N] = self.log_scale[iterations]
                N += 1

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
    pixel of value 0.

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
    # Allocate memory
    buf = <char *>malloc(length*sizeof(char));
    # Compute dimension
    result = []
    for i in range(length):
        buf[i] = 1 if in_string[i] == 0 else 0
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
